//! A borrow-guarded shared reader.
//!
//! Wraps a single underlying `R` (here a Python file-like object, which cannot be
//! duplicated) so it can be lent to iterators safely. One `GuardedReader<R>`
//! type: the handle from [`GuardedReader::new`] is the *owner*; [`Reopen::reopen`]
//! mints a *borrow* sharing the same `R`. A single type is what satisfies
//! bigtools' [`Reopen`] trait, so existing call sites (`b.reopen()`) work
//! unchanged.
//!
//! One shared `state` atom (`OWNER` | `BUSY` | a live-borrow epoch) coordinates:
//!   * At most one handle is live. `reopen` makes the new borrow live and
//!     supersedes the previous one; an owner read reclaims (`-> OWNER`).
//!   * A superseded borrow is dead forever — its acquire `compare_exchange`
//!     fails (the RMW reads the latest epoch, so this is deterministic).
//!   * Borrow op: `CAS(epoch -> BUSY)` to acquire, `Release` store to release.
//!     The owner *spins* the (rare) take-back-mid-IO wait — no condvar.
//!   * Owner op: fast path is one relaxed load (`state == OWNER`); otherwise it
//!     reclaims, spinning out any in-flight borrow.
//!
//! A read is therefore always either the live cursor's correct bytes or an error
//! — never another cursor's bytes — even without the GIL.
//!
//! Soundness rests on `reopen` (i.e. creating an iterator) and owner reads being
//! serialized against each other, which pyo3's `&mut self` borrow rules give us:
//! `records()`/`average_over_bed()` and `values()`/`info()` are all `&mut self`
//! on the reader. Driving already-created iterators on other threads is free.

use std::cell::UnsafeCell;
use std::io::{self, Read, Seek, SeekFrom};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use bigtools::utils::reopen::Reopen;

const OWNER: u64 = 0;
const BUSY: u64 = u64::MAX;

struct Shared<R> {
    inner: UnsafeCell<R>,
    state: AtomicU64,
    /// Next borrow epoch. Shared so `reopen(&self)` can mint without `&mut`.
    dispenser: AtomicU64,
}

// SAFETY: `inner` is only touched by the handle that currently holds `state`
// (OWNER for the owner, or BUSY for a borrow), so there is never aliasing `&mut R`.
unsafe impl<R: Send> Sync for Shared<R> {}

pub struct GuardedReader<R> {
    shared: Arc<Shared<R>>,
    epoch: u64,
}

impl<R> GuardedReader<R> {
    /// Create the owner handle around `inner`.
    pub fn new(inner: R) -> Self {
        GuardedReader {
            shared: Arc::new(Shared {
                inner: UnsafeCell::new(inner),
                state: AtomicU64::new(OWNER),
                dispenser: AtomicU64::new(1),
            }),
            epoch: OWNER,
        }
    }

    /// Mint a fresh borrow over the same reader, superseding any previous borrow.
    fn mint(&self) -> Self {
        let epoch = self.shared.dispenser.fetch_add(1, Ordering::Relaxed);
        let mut cur = self.shared.state.load(Ordering::Relaxed);
        loop {
            if cur == BUSY {
                std::hint::spin_loop();
                cur = self.shared.state.load(Ordering::Relaxed);
                continue;
            }
            match self.shared.state.compare_exchange(
                cur,
                epoch,
                Ordering::AcqRel,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(c) => cur = c,
            }
        }
        GuardedReader {
            shared: self.shared.clone(),
            epoch,
        }
    }

    /// (Owner only) ensure `state == OWNER`, spinning out any in-flight borrow.
    fn reclaim(&self) {
        if self.shared.state.load(Ordering::Relaxed) == OWNER {
            return; // fast path: no live borrow
        }
        let mut cur = self.shared.state.load(Ordering::Acquire);
        loop {
            if cur == OWNER {
                return;
            }
            if cur == BUSY {
                std::hint::spin_loop();
                cur = self.shared.state.load(Ordering::Acquire);
                continue;
            }
            match self.shared.state.compare_exchange(
                cur,
                OWNER,
                Ordering::AcqRel,
                Ordering::Acquire,
            ) {
                Ok(_) => return,
                Err(c) => cur = c,
            }
        }
    }

    fn with_io<T>(&self, f: impl FnOnce(&mut R) -> io::Result<T>) -> io::Result<T> {
        if self.epoch == OWNER {
            self.reclaim();
            // SAFETY: state == OWNER and (per the &mut-self contract) no borrow can
            // be minted concurrently, so no borrow holds BUSY here.
            let inner = unsafe { &mut *self.shared.inner.get() };
            f(inner)
        } else {
            if self
                .shared
                .state
                .compare_exchange(self.epoch, BUSY, Ordering::Acquire, Ordering::Relaxed)
                .is_err()
            {
                return Err(invalidated());
            }
            // SAFETY: we hold BUSY ⇒ exclusive access.
            let inner = unsafe { &mut *self.shared.inner.get() };
            let r = f(inner);
            self.shared.state.store(self.epoch, Ordering::Release);
            r
        }
    }
}

impl<R> Reopen for GuardedReader<R> {
    fn reopen(&self) -> io::Result<Self> {
        Ok(self.mint())
    }
}

impl<R: Read + Seek> Read for GuardedReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.with_io(|r| r.read(buf))
    }
}

impl<R: Read + Seek> Seek for GuardedReader<R> {
    fn seek(&mut self, pos: SeekFrom) -> io::Result<u64> {
        self.with_io(|r| r.seek(pos))
    }
}

fn invalidated() -> io::Error {
    io::Error::new(
        io::ErrorKind::Other,
        "reader handle invalidated: the file is in use by another cursor (create a \
         fresh reader, or finish/close the outstanding iterator before reusing it)",
    )
}
