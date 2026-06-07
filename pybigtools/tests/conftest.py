"""Shared pytest fixtures and helpers.

Some tests fetch fixtures from live servers (ENCODE, UCSC) that intermittently
return 403s or 5xx errors. The ``require_remote`` fixture lets those tests skip
on transient server problems instead of failing, while still failing on genuine
client errors (e.g. a 404) or on assertion failures once the data is fetched.
"""

import urllib.error
import urllib.request

import pytest

# HTTP statuses we treat as transient/server-side: skip rather than fail.
TRANSIENT_STATUSES = frozenset({403, 408, 425, 429, 500, 502, 503, 504})


def _require_remote(url, *, timeout=30):
    """Skip the current test unless ``url`` is reachable.

    Performs a 1-byte ranged GET (rather than HEAD, which some servers reject)
    and follows redirects. Connection failures, timeouts, and transient/
    forbidden HTTP statuses cause a skip; any other HTTP error (e.g. 404)
    propagates so real problems still fail the test.
    """
    request = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
    try:
        urllib.request.urlopen(request, timeout=timeout).close()
    except urllib.error.HTTPError as exc:
        if exc.code in TRANSIENT_STATUSES:
            pytest.skip(f"remote server returned {exc.code} for {url}")
        if exc.code == 416:  # range not satisfiable -> server is up, proceed
            return
        raise
    except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
        pytest.skip(f"remote server unreachable for {url}: {exc}")


@pytest.fixture
def require_remote():
    """Return a callable that skips the test if a URL is unavailable.

    Usage::

        def test_something(require_remote):
            require_remote(url)
            ...  # network-dependent assertions
    """
    return _require_remote
