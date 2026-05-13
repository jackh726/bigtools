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


@pytest.fixture
def big_bw_path():
    """Real bigWig with 10 zoom levels (chr17, 83M bases); used to exercise
    the zoom-interpolation code path (`exact=False` with non-trivial bin sizes)."""
    return str(REPO_ROOT / "bigtools/resources/test/valid.bigWig")


@pytest.fixture
def big_bb_path():
    """Real bigBed with 5 zoom levels (chr21, 48M bases)."""
    return str(TEST_DIR / "data/bigBedExample.bb")


def test_values(bw_path, bb_path):
    # Default `uncovered=None` (NaN): uncovered bases are NaN.
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values("chr1", 0, 12)
        assert np.allclose(
            x1,
            [np.nan, 1, 1, np.nan, 2, 2, 2, 2, 3, 3, np.nan, np.nan],
            equal_nan=True,
        )

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values("chr1", -2, 14)
        assert np.allclose(
            x1,
            [
                np.nan,
                np.nan,
                np.nan,
                1,
                1,
                np.nan,
                2,
                2,
                2,
                2,
                3,
                3,
                np.nan,
                np.nan,
                np.nan,
                np.nan,
            ],
            equal_nan=True,
        )

    # `uncovered=0` reproduces the previous fill-with-zero behavior.
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values("chr1", 0, 12, uncovered=0)
        assert np.allclose(x1, [0, 1, 1, 0, 2, 2, 2, 2, 3, 3, 0, 0], equal_nan=True)


def test_fillna(bw_path, bb_path):
    # `fillna` is a post-rasterization fill that replaces NaN at in-bounds
    # positions. It is INDEPENDENT of `oob`: out-of-chromosome positions are
    # filled by `oob` and not touched by `fillna`.

    # Base-level: default `uncovered=None` leaves uncovered bases as NaN;
    # `fillna=-99` replaces them after the fact.
    for path in [bw_path, bb_path]:
        x = pbt.open(path).values("chr1", 0, 12, fillna=-99)
        assert np.allclose(
            x, [-99, 1, 1, -99, 2, 2, 2, 2, 3, 3, -99, -99], equal_nan=True
        )

    # `oob` and `fillna` control different concerns: OOB positions take the
    # `oob` value; in-bounds uncovered positions take `fillna`.
    for path in [bw_path, bb_path]:
        x = pbt.open(path).values("chr1", -2, 14, oob=-1, fillna=0)
        assert np.allclose(
            x,
            [-1, -1, 0, 1, 1, 0, 2, 2, 2, 2, 3, 3, 0, 0, -1, -1],
            equal_nan=True,
        )

    # `fillna` does not "leak" into OOB positions: with `oob=NaN` (default),
    # OOB positions stay NaN even when `fillna` is set.
    for path in [bw_path, bb_path]:
        x = pbt.open(path).values("chr1", -2, 14, fillna=0)
        assert np.allclose(
            x,
            [
                np.nan,
                np.nan,
                0,
                1,
                1,
                0,
                2,
                2,
                2,
                2,
                3,
                3,
                0,
                0,
                np.nan,
                np.nan,
            ],
            equal_nan=True,
        )

    # Binned: empty bin under `uncovered=None` (default) returns NaN; `fillna=0`
    # replaces that. Partial bins are unaffected because they returned a
    # real number.
    for path in [bw_path, bb_path]:
        x = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, summary="mean", fillna=0
        )
        # Bins 0..4 partial or full → real means; bin 5 empty → fillna=0.
        assert np.allclose(x, [1, 1, 2, 2, 3, 0], equal_nan=True)


def test_binned_mean(bw_path, bb_path):
    # uncovered=None: classical mean over only covered bases; empty bins → NaN.
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=12, exact=True, uncovered=None, summary="mean"
        )
        assert np.allclose(
            x1,
            [np.nan, 1, 1, np.nan, 2, 2, 2, 2, 3, 3, np.nan, np.nan],
            equal_nan=True,
        )

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, uncovered=None, summary="mean"
        )
        assert np.allclose(x1, [1, 1, 2, 2, 3, np.nan], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=5, exact=True, uncovered=None, summary="mean"
        )
        assert np.allclose(x1, [1, 1, 2, 2.5, 3], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=2, exact=True, uncovered=None, summary="mean"
        )
        assert np.allclose(x1, [1.5, 2.5], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=1, exact=True, uncovered=None, summary="mean"
        )
        assert np.allclose(x1, [2], equal_nan=True)


def test_binned_mean_uncovered_zero(bw_path, bb_path):
    # uncovered=0: treat uncovered bases as 0 (UCSC `mean0` semantics).
    # Equivalent to (sum + 0*n_uncovered) / bin_size = sum / bin_size.
    uncovered = 0
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=12, exact=True, uncovered=uncovered, summary="mean"
        )
        assert np.allclose(
            x1,
            [uncovered, 1, 1, uncovered, 2, 2, 2, 2, 3, 3, uncovered, uncovered],
            equal_nan=True,
        )

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, uncovered=uncovered, summary="mean"
        )
        assert np.allclose(x1, [0.5, 0.5, 2, 2, 3, uncovered], equal_nan=True)

    # bins=5 has non-integer `bin_size = 2.4`. Integer-aligned bin widths are
    # [2, 2, 3, 2, 3] — see `integer_bin_bounds` in raster.rs. Mean0 uses the
    # integer bin width (not the float `bin_size`) as the denominator.
    # Bin widths: 2     2     3       2       3
    # Coverage:   1/2   1/2   3/3     2/2     1/3
    # Sums:       1     1     6       5       3
    # mean0:      0.5   0.5   2.0     2.5     1.0
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=5, exact=True, uncovered=uncovered, summary="mean"
        )
        assert np.allclose(x1, [0.5, 0.5, 2.0, 2.5, 1.0], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=2, exact=True, uncovered=uncovered, summary="mean"
        )
        assert np.allclose(x1, [1, 10 / 6], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=1, exact=True, uncovered=uncovered, summary="mean"
        )
        assert np.allclose(x1, [16 / 12], equal_nan=True)


def test_binned_minmax(bw_path, bb_path):
    # uncovered=None: classical min/max over only covered bases; empty bins → NaN.
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=12, exact=True, uncovered=None, summary="min"
        )
        assert np.allclose(
            x1,
            [np.nan, 1, 1, np.nan, 2, 2, 2, 2, 3, 3, np.nan, np.nan],
            equal_nan=True,
        )
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=12, exact=True, uncovered=None, summary="max"
        )
        assert np.allclose(
            x1,
            [np.nan, 1, 1, np.nan, 2, 2, 2, 2, 3, 3, np.nan, np.nan],
            equal_nan=True,
        )

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, uncovered=None, summary="min"
        )
        assert np.allclose(x1, [1, 1, 2, 2, 3, np.nan], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, uncovered=None, summary="max"
        )
        assert np.allclose(x1, [1, 1, 2, 2, 3, np.nan], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=5, exact=True, uncovered=None, summary="min"
        )
        assert np.allclose(x1, [1, 1, 2, 2, 3], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=5, exact=True, uncovered=None, summary="max"
        )
        assert np.allclose(x1, [1, 1, 2, 3, 3], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=2, exact=True, uncovered=None, summary="min"
        )
        assert np.allclose(x1, [1, 2], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=2, exact=True, uncovered=None, summary="max"
        )
        assert np.allclose(x1, [2, 3], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=1, exact=True, uncovered=None, summary="min"
        )
        assert np.allclose(x1, [1], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=1, exact=True, uncovered=None, summary="max"
        )
        assert np.allclose(x1, [3], equal_nan=True)


def test_binned_minmax_uncovered_zero(bw_path, bb_path):
    # uncovered=0: uncovered bases counted as 0 in min/max (UCSC `min0`/`max0` semantics).
    uncovered = 0
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=12, exact=True, uncovered=uncovered, summary="min"
        )
        assert np.allclose(
            x1,
            [uncovered, 1, 1, uncovered, 2, 2, 2, 2, 3, 3, uncovered, uncovered],
            equal_nan=True,
        )
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=12, exact=True, uncovered=uncovered, summary="max"
        )
        assert np.allclose(
            x1,
            [uncovered, 1, 1, uncovered, 2, 2, 2, 2, 3, 3, uncovered, uncovered],
            equal_nan=True,
        )

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, uncovered=uncovered, summary="min"
        )
        assert np.allclose(x1, [0, 0, 2, 2, 3, uncovered], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=6, exact=True, uncovered=uncovered, summary="max"
        )
        assert np.allclose(x1, [1, 1, 2, 2, 3, uncovered], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=5, exact=True, uncovered=uncovered, summary="min"
        )
        assert np.allclose(x1, [0, 0, 2, 2, 0], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=5, exact=True, uncovered=uncovered, summary="max"
        )
        assert np.allclose(x1, [1, 1, 2, 3, 3], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=2, exact=True, uncovered=uncovered, summary="min"
        )
        assert np.allclose(x1, [0, 0], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=2, exact=True, uncovered=uncovered, summary="max"
        )
        assert np.allclose(x1, [2, 3], equal_nan=True)

    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=1, exact=True, uncovered=uncovered, summary="min"
        )
        assert np.allclose(x1, [0], equal_nan=True)
        x1 = pbt.open(path).values(
            "chr1", 0, 12, bins=1, exact=True, uncovered=uncovered, summary="max"
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


def test_values_arr_param(bw_path, bb_path):
    # Pass a pre-allocated array. The function writes into it and returns the
    # same object.
    for path in [bw_path, bb_path]:
        arr = np.full(12, -7.0)
        ret = pbt.open(path).values("chr1", 0, 12, uncovered=0, arr=arr)
        assert ret is arr  # same object
        assert np.allclose(arr, [0, 1, 1, 0, 2, 2, 2, 2, 3, 3, 0, 0])

    # Binned: pre-allocated array of the correct length.
    for path in [bw_path, bb_path]:
        arr = np.full(6, -7.0)
        ret = pbt.open(path).values(
            "chr1",
            0,
            12,
            bins=6,
            exact=True,
            summary="mean",
            uncovered=None,
            arr=arr,
        )
        assert ret is arr
        assert np.allclose(arr, [1, 1, 2, 2, 3, np.nan], equal_nan=True)


def test_values_arr_wrong_size(bw_path):
    # Wrong-length `arr` raises ValueError.
    with pytest.raises(ValueError):
        pbt.open(bw_path).values("chr1", 0, 12, arr=np.zeros(11))  # need 12

    with pytest.raises(ValueError):
        pbt.open(bw_path).values(
            "chr1", 0, 12, bins=6, exact=True, arr=np.zeros(5)
        )  # need 6


def test_values_bad_args(bw_path):
    # Unknown chromosome.
    with pytest.raises(KeyError):
        pbt.open(bw_path).values("nope", 0, 12)

    # Unknown summary string.
    with pytest.raises(ValueError):
        pbt.open(bw_path).values("chr1", 0, 12, bins=4, summary="median")


# -- These tests check Parity with UCSC via pybbi --


@pytest.mark.parametrize(
    "start,end", [(start, end) for start in [0, 2] for end in [10, 12]]
)
def test_bbi_consistent_values(bw_path, bb_path, start, end):
    # pybigtools defaults to `uncovered=None` (NaN), pybbi to `missing=0.0`.
    # Pass matching values explicitly for the comparison.
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values("chr1", start, end, uncovered=0.0)
        x2 = bbi.open(path).fetch("chr1", start, end, missing=0.0)
        assert np.allclose(x1, x2, equal_nan=True)


def test_bbi_consistent_values_match_coverage(bb_path, bw_path):
    bb1 = pbt.open(bb_path).values("chr1", uncovered=0.0)
    bw1 = pbt.open(bw_path).values("chr1", uncovered=0.0)
    assert np.allclose(bb1, bw1, equal_nan=True)

    bb2 = bbi.open(bb_path).fetch("chr1", 0, -1)
    bw2 = bbi.open(bw_path).fetch("chr1", 0, -1)
    assert np.allclose(bb2, bw2, equal_nan=True)


@pytest.mark.parametrize(
    "uncovered,oob",
    [(uncovered, oob) for uncovered in [None, 0, -1] for oob in [np.nan, 0, -1]],
)
@pytest.mark.parametrize(
    "start,end",
    [(start, end) for start in [-2, 0, 2] for end in [10, 12, 14]],
)
def test_bbi_consistent_values_uncovered_oob(
    bw_path, bb_path, start, end, uncovered, oob
):
    # `uncovered=None` (pybigtools) ≡ `missing=NaN` (pybbi); both mean "leave
    # uncovered bases as NaN."
    pbbi_missing = np.nan if uncovered is None else uncovered
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values("chr1", start, end, uncovered=uncovered, oob=oob)
        x2 = bbi.open(path).fetch("chr1", start, end, missing=pbbi_missing, oob=oob)
        assert np.allclose(x1, x2, equal_nan=True)


@pytest.mark.parametrize(
    "start,end,n_bins,exact,summary",
    [
        (start, end, n_bins, exact, summary)
        for start in [0, 2]
        for end in [10, 12]
        for n_bins in range(1, 13)
        for exact in [True]
        # `bin_covered` (= `cov` in pybbi) is excluded: pybigtools always returns
        # the factual coverage fraction for empty bins (0.0), while pybbi fills
        # empty bins with `missing`. Intentional divergence.
        for summary in ["sum", "mean", "std", "min", "max"]
    ],
)
def test_bbi_consistent_binned(bw_path, bb_path, start, end, n_bins, exact, summary):
    # pybbi computes "regular" statistics (exclude uncovered) and fills empty bins
    # with `missing`. Under pybigtools' per-base `uncovered` semantics, the
    # equivalent comparison passes `uncovered=None` to pybigtools and
    # `missing=NaN` to pybbi: classical stats over partial bins, NaN for empty.
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1",
            start,
            end,
            bins=n_bins,
            exact=exact,
            summary=summary,
            uncovered=None,
        )
        x2 = bbi.open(path).fetch(
            "chr1",
            start,
            end,
            bins=n_bins,
            exact=exact,
            summary=summary,
            missing=np.nan,
        )
        assert np.allclose(x1, x2, equal_nan=True)


@pytest.mark.parametrize(
    "start,end,n_bins,exact,oob",
    [
        (start, end, n_bins, exact, oob)
        for start in [-2, 0, 2]
        for end in [10, 12, 14]
        for n_bins in range(1, 13)
        for exact in [True]
        for oob in [np.nan, 0, -1]
    ],
)
def test_bbi_consistent_binned_uncovered_oob(
    bw_path, bb_path, start, end, n_bins, exact, oob
):
    # `uncovered` is restricted to None for parity: pybigtools' per-base `uncovered`
    # semantics diverge from pybbi's empty-bin-fill semantics whenever the bin is
    # partially covered. With pybigtools' `uncovered=None` and pybbi's
    # `missing=NaN`, both compute classical statistics over the covered bases.
    for path in [bw_path, bb_path]:
        x1 = pbt.open(path).values(
            "chr1", start, end, bins=n_bins, exact=exact, uncovered=None, oob=oob
        )
        x2 = bbi.open(path).fetch(
            "chr1", start, end, bins=n_bins, exact=exact, missing=np.nan, oob=oob
        )
        assert np.allclose(x1, x2, equal_nan=True)


@pytest.mark.parametrize(
    "start,end,n_bins,exact",
    [(0, 12, n_bins, exact) for n_bins in [12, 24, 48] for exact in [True]],
)
def test_bbi_consistent_binned_upsample(bw_path, bb_path, start, end, n_bins, exact):
    # Upsampling produces only fully-covered or fully-empty bins (no partial bins),
    # so per-base vs empty-bin `uncovered` semantics agree.
    for path in [bb_path, bw_path]:
        x1 = pbt.open(path).values(
            "chr1", start, end, bins=n_bins, exact=exact, uncovered=None
        )
        x2 = bbi.open(path).fetch(
            "chr1", start, end, bins=n_bins, exact=exact, missing=np.nan
        )
        assert np.allclose(x1, x2, equal_nan=True)


@pytest.mark.parametrize(
    "start,end,n_bins,summary",
    [
        # Large query ranges + few bins → bin size big enough that `best_zoom`
        # selects a real zoom level, exercising the zoom-interpolation path.
        # `valid.bigWig` has zoom resolutions [10, 40, 160, 640, 2560, ...].
        # For a query of width W and N bins, max_zoom_size = W / (2*N).
        # E.g. 100kb / (2*10 bins) = 5000 ≥ resolution 2560 → uses zoom 2560.
        (100_000, 200_000, 10, "mean"),
        (100_000, 200_000, 10, "std"),
        (100_000, 200_000, 10, "min"),
        (100_000, 200_000, 10, "max"),
        (100_000, 200_000, 10, "sum"),
        (100_000, 200_000, 50, "mean"),
        (0, 1_000_000, 100, "mean"),
        (0, 1_000_000, 100, "std"),
    ],
)
def test_bbi_consistent_binned_zoom_bw(big_bw_path, start, end, n_bins, summary):
    # Parity vs pybbi when traversing the zoom-interpolation path. This is the
    # path that the UCSC ceil+rescale normalization affects, so a passing
    # parity test here validates that normalization is correctly applied.
    x1 = pbt.open(big_bw_path).values(
        "chr17",
        start,
        end,
        bins=n_bins,
        exact=False,
        summary=summary,
        uncovered=None,
    )
    x2 = bbi.open(big_bw_path).fetch(
        "chr17",
        start,
        end,
        bins=n_bins,
        exact=False,
        summary=summary,
        missing=np.nan,
    )
    assert np.allclose(x1, x2, equal_nan=True)


@pytest.mark.parametrize(
    "start,end,n_bins,summary",
    [
        (10_000_000, 20_000_000, 50, "mean"),
        (10_000_000, 20_000_000, 50, "min"),
        (10_000_000, 20_000_000, 50, "max"),
        (10_000_000, 20_000_000, 100, "sum"),
    ],
)
def test_bbi_consistent_binned_zoom_bb(big_bb_path, start, end, n_bins, summary):
    # Same parity check on the bigBed zoom path.
    x1 = pbt.open(big_bb_path).values(
        "chr21",
        start,
        end,
        bins=n_bins,
        exact=False,
        summary=summary,
        uncovered=None,
    )
    x2 = bbi.open(big_bb_path).fetch(
        "chr21",
        start,
        end,
        bins=n_bins,
        exact=False,
        summary=summary,
        missing=np.nan,
    )
    assert np.allclose(x1, x2, equal_nan=True)


def test_values_fully_out_of_bounds(bw_path, bb_path):
    # Query entirely past the end of the chromosome.
    # chr1 length is 12. Query [20, 30) is fully OOB.
    for path in [bw_path, bb_path]:
        x = pbt.open(path).values("chr1", 20, 30)
        # Default oob=NaN.
        assert len(x) == 10
        assert all(np.isnan(v) for v in x)

    # With oob=-1, the whole array is -1.
    for path in [bw_path, bb_path]:
        x = pbt.open(path).values("chr1", 20, 30, oob=-1)
        assert np.array_equal(x, np.full(10, -1.0))

    # Query entirely before the start (negative range).
    for path in [bw_path, bb_path]:
        x = pbt.open(path).values("chr1", -10, -2)
        assert len(x) == 8
        assert all(np.isnan(v) for v in x)

    # Binned version of the same — all bins OOB, all bins get oob.
    for path in [bw_path, bb_path]:
        x = pbt.open(path).values("chr1", 20, 30, bins=5, exact=True, oob=0)
        assert np.array_equal(x, np.zeros(5))
