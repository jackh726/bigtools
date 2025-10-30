import pathlib

import numpy as np
import pytest

import pybigtools as pbt
import bbi

TEST_DIR = pathlib.Path(__file__).parent
REPO_ROOT = TEST_DIR.parent.parent


@pytest.fixture
def bw_path():
    """
    The values here are identical to the coverage depth of ``bb_path``.
    The coverage breadth is 8 out of 12 bases (66.67%).

    chroms = {"chr1": 12}
    entries = [
        ("chr1", 1, 3, 1.0),
        ("chr1", 4, 8, 2.0),
        ("chr1", 8, 10, 3.0)
    ]
    """
    return str(TEST_DIR / "data/mini.bw")


@pytest.fixture
def bb_path():
    """
    The coverage depth of these intervals is identical to ``bw_path``.
    The coverage breadth is 8 out of 12 bases (66.67%).

    chroms = {"chr1": 12}
    entries = [
        ("chr1", 1, 3, ""),
        ("chr1", 4, 6, ""),
        ("chr1", 4, 10, ""),
        ("chr1", 6, 10, ""),
        ("chr1", 8, 10, ""),
    ]
    """
    return str(TEST_DIR / "data/mini.bb")


def test_values(bw_path, bb_path):
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values("chr1", 0, 12)
        assert np.allclose(x1, [0, 1, 1, 0, 2, 2, 2, 2, 3, 3, 0, 0], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values("chr1", -2, 14)
        assert np.allclose(
            x1,
            [np.nan, np.nan, 0, 1, 1, 0, 2, 2, 2, 2, 3, 3, 0, 0, np.nan, np.nan],
            equal_nan=True,
        )


@pytest.mark.parametrize("missing", [np.nan, 0, -1])
def test_binned_mean(bw_path, bb_path, missing):
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=12, exact=True, missing=missing, summary="mean"
        )
        assert np.allclose(
            x1,
            [missing, 1, 1, missing, 2, 2, 2, 2, 3, 3, missing, missing],
            equal_nan=True,
        )

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, missing=missing, summary="mean"
        )
        assert np.allclose(x1, [1, 1, 2, 2, 3, missing], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=5, exact=True, missing=missing, summary="mean"
        )
        assert np.allclose(x1, [1, 1, 2, 2.5, 3], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=2, exact=True, missing=missing, summary="mean"
        )
        assert np.allclose(x1, [1.5, 2.5], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=1, exact=True, missing=missing, summary="mean"
        )
        assert np.allclose(x1, [2], equal_nan=True)


@pytest.mark.parametrize("missing", [np.nan, 0, -1])
def test_binned_mean0(bw_path, bb_path, missing):
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=12, exact=True, missing=missing, summary="mean0"
        )
        assert np.allclose(
            x1,
            [missing, 1, 1, missing, 2, 2, 2, 2, 3, 3, missing, missing],
            equal_nan=True,
        )

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, missing=missing, summary="mean0"
        )
        assert np.allclose(x1, [0.5, 0.5, 2, 2, 3, missing], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=5, exact=True, missing=missing, summary="mean0"
        )
        assert np.allclose(
            x1, [1 / 2.4, 1 / 2.4, 6 / 2.4, 5 / 2.4, 3 / 2.4], equal_nan=True
        )

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=2, exact=True, missing=missing, summary="mean0"
        )
        assert np.allclose(x1, [1, 10 / 6], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=1, exact=True, missing=missing, summary="mean0"
        )
        assert np.allclose(x1, [16 / 12], equal_nan=True)


@pytest.mark.parametrize("missing", [np.nan, 0, -1])
def test_binned_minmax(bw_path, bb_path, missing):
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=12, exact=True, missing=missing, summary="min"
        )
        assert np.allclose(
            x1,
            [missing, 1, 1, missing, 2, 2, 2, 2, 3, 3, missing, missing],
            equal_nan=True,
        )
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=12, exact=True, missing=missing, summary="max"
        )
        assert np.allclose(
            x1,
            [missing, 1, 1, missing, 2, 2, 2, 2, 3, 3, missing, missing],
            equal_nan=True,
        )

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, missing=missing, summary="min"
        )
        assert np.allclose(x1, [1, 1, 2, 2, 3, missing], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, missing=missing, summary="max"
        )
        assert np.allclose(x1, [1, 1, 2, 2, 3, missing], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=5, exact=True, missing=missing, summary="min"
        )
        assert np.allclose(x1, [1, 1, 2, 2, 3], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=5, exact=True, missing=missing, summary="max"
        )
        assert np.allclose(x1, [1, 1, 2, 3, 3], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=2, exact=True, missing=missing, summary="min"
        )
        assert np.allclose(x1, [1, 2], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=2, exact=True, missing=missing, summary="max"
        )
        assert np.allclose(x1, [2, 3], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=1, exact=True, missing=missing, summary="min"
        )
        assert np.allclose(x1, [1], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=1, exact=True, missing=missing, summary="max"
        )
        assert np.allclose(x1, [3], equal_nan=True)


@pytest.mark.parametrize("missing", [np.nan, 0, -1])
def test_binned_minmax0(bw_path, bb_path, missing):
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=12, exact=True, missing=missing, summary="min0"
        )
        assert np.allclose(
            x1,
            [missing, 1, 1, missing, 2, 2, 2, 2, 3, 3, missing, missing],
            equal_nan=True,
        )
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=12, exact=True, missing=missing, summary="max0"
        )
        assert np.allclose(
            x1,
            [missing, 1, 1, missing, 2, 2, 2, 2, 3, 3, missing, missing],
            equal_nan=True,
        )

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, missing=missing, summary="min0"
        )
        assert np.allclose(x1, [0, 0, 2, 2, 3, missing], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, missing=missing, summary="max0"
        )
        assert np.allclose(x1, [1, 1, 2, 2, 3, missing], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=5, exact=True, missing=missing, summary="min0"
        )
        assert np.allclose(x1, [0, 0, 2, 2, 0], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=5, exact=True, missing=missing, summary="max0"
        )
        assert np.allclose(x1, [1, 1, 2, 3, 3], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=2, exact=True, missing=missing, summary="min0"
        )
        assert np.allclose(x1, [0, 0], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=2, exact=True, missing=missing, summary="max0"
        )
        assert np.allclose(x1, [2, 3], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=1, exact=True, missing=missing, summary="min0"
        )
        assert np.allclose(x1, [0], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=1, exact=True, missing=missing, summary="max0"
        )
        assert np.allclose(x1, [3], equal_nan=True)


@pytest.mark.parametrize(
    "n_bins,exact", [(n_bins, exact) for n_bins in range(1, 13) for exact in [True]]
)
def test_binned_coverage_breadth(bw_path, bb_path, n_bins, exact):
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=n_bins, exact=exact, summary="bases_covered"
        )
        assert np.nansum(x1) == 8


@pytest.mark.parametrize(
    "start,end", [(start, end) for start in [0, 2] for end in [10, 12]]
)
def test_bbi_consistent_values(bw_path, bb_path, start, end):
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values("chr1", start, end)
        x2 = bbi.open(path).fetch("chr1", start, end)
        assert np.allclose(x1, x2, equal_nan=True)


def test_bbi_consistent_values_match_coverage(bb_path, bw_path):
    bb1 = pbt.open(bb_path).values("chr1")
    bw1 = pbt.open(bw_path).values("chr1")
    assert np.allclose(bb1, bw1, equal_nan=True)

    bb2 = bbi.open(bb_path).fetch("chr1", 0, -1)
    bw2 = bbi.open(bw_path).fetch("chr1", 0, -1)
    assert np.allclose(bb2, bw2, equal_nan=True)


@pytest.mark.parametrize(
    "start,end,missing,oob",
    [
        (start, end, missing, oob)
        for start in [-2, 0, 2]
        for end in [10, 12, 14]
        for missing in [np.nan, 0, -1]
        for oob in [np.nan, 0, -1]
    ],
)
def test_bbi_consistent_values_missing_oob(bw_path, bb_path, start, end, missing, oob):
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values("chr1", start, end, missing=missing, oob=oob)
        x2 = bbi.open(path).fetch("chr1", start, end, missing=missing, oob=oob)
        assert np.allclose(x1, x2, equal_nan=True)


@pytest.mark.parametrize(
    "start,end,n_bins,exact,summary",
    [
        (start, end, n_bins, exact, summary)
        for start in [0, 2]
        for end in [10, 12]
        for n_bins in range(1, 13)
        for exact in [True]
        for summary in ["sum", "mean", "std", "min", "max", "bin_covered"]
    ],
)
def test_bbi_consistent_binned(bw_path, bb_path, start, end, n_bins, exact, summary):
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", start, end, bins=n_bins, exact=exact, summary=summary
        )
        x2 = bbi.open(path).fetch(
            "chr1",
            start,
            end,
            bins=n_bins,
            exact=exact,
            summary=summary if summary != "bin_covered" else "cov",
        )
        assert np.allclose(x1, x2, equal_nan=True)


@pytest.mark.parametrize(
    "start,end,n_bins,exact,missing,oob",
    [
        (start, end, n_bins, exact, missing, oob)
        for start in [-2, 0, 2]
        for end in [10, 12, 14]
        for n_bins in range(1, 13)
        for exact in [True]
        for missing in [np.nan, 0, -1]
        for oob in [np.nan, 0, -1]
    ],
)
def test_bbi_consistent_binned_missing_oob(
    bw_path, bb_path, start, end, n_bins, exact, missing, oob
):
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", start, end, bins=n_bins, exact=exact, missing=missing, oob=oob
        )
        x2 = bbi.open(path).fetch(
            "chr1", start, end, bins=n_bins, exact=exact, missing=missing, oob=oob
        )
        assert np.allclose(x1, x2, equal_nan=True)


@pytest.mark.parametrize(
    "start,end,n_bins,exact",
    [(0, 12, n_bins, exact) for n_bins in [12, 24, 48] for exact in [True]],
)
def test_bbi_consistent_binned_upsample(bw_path, bb_path, start, end, n_bins, exact):
    for path in [bb_path, bw_path]:
        x1 = pbt.open(path).values("chr1", start, end, bins=n_bins, exact=exact)
        x2 = bbi.open(path).fetch("chr1", start, end, bins=n_bins, exact=exact)
        assert np.allclose(x1, x2, equal_nan=True)
