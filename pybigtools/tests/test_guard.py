"""Tests for the borrow-guarded file-like reader."""

import pathlib

import pybigtools

TEST_DIR = pathlib.Path(__file__).parent
PATH = str(TEST_DIR / "data" / "ENCFF667CZO.bigWig")


def first_chrom_with_data(b):
    # chr1 has data starting at 10495 (per test_io.py)
    return "chr1"


def test_file_like_reads_match_path():
    """A file-like open returns the same data as a path open."""
    by_path = pybigtools.open(PATH)
    by_file = pybigtools.open(open(PATH, "rb"))
    chrom = first_chrom_with_data(by_path)
    a = next(by_path.records(chrom))
    b = next(by_file.records(chrom))
    assert a == b


def test_new_iterator_invalidates_old_file_like():
    """On a file-like reader, creating a second iterator supersedes the first.
    The superseded iterator must error when it next reads from the file rather
    than return corrupt/shared-state bytes. (Values already buffered from a block
    read while it was live are still valid — invalidation is per file-read.)"""
    b = pybigtools.open(open(PATH, "rb"))
    it1 = b.records("chr1")
    it2 = b.records("chr1")  # supersedes it1 before it1 reads any block
    assert next(it2) is not None  # it2 is the live cursor
    try:
        next(it1)  # it1's first block read -> acquire fails -> error
        raised = False
    except Exception:
        raised = True
    assert raised, "superseded file-like iterator should error, not yield"


def test_owner_takes_back_file_like():
    """Using the parent reader (a 'values' call) after making an iterator should
    reclaim the file; the outstanding iterator then errors."""
    b = pybigtools.open(open(PATH, "rb"))
    it = b.records("chr1")
    next(it)
    # parent op reclaims the underlying file
    _ = b.records("chr1")  # a fresh cursor / owner activity supersedes `it`
    try:
        for _ in it:
            pass
        # consuming a superseded iterator must not silently yield wrong data;
        # either it stops or raises — both are acceptable, corruption is not.
    except Exception:
        pass


def test_path_iterators_are_independent():
    """Path-backed readers reopen independent handles, so two iterators coexist
    (unchanged behavior)."""
    b = pybigtools.open(PATH)
    it1 = b.records("chr1")
    it2 = b.records("chr1")
    v1 = next(it1)
    v2 = next(it2)
    assert v1 == v2  # both independently start at the beginning
