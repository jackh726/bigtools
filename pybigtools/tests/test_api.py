import math
import pathlib
from io import BytesIO

import numpy as np
import pytest
import smart_open

import pybigtools

TEST_DIR = pathlib.Path(__file__).parent
REPO_ROOT = TEST_DIR.parent.parent


def test_open_close():
    path = str(REPO_ROOT / "bigtools/resources/test/valid.bigWig")
    b = pybigtools.open(path, "r")
    b.close()
    assert pytest.raises(pybigtools.BBIFileClosed, b.chroms)

    # Files are closed when exiting a context manager
    with pybigtools.open(path, "r") as b:
        pass
    assert pytest.raises(pybigtools.BBIFileClosed, b.chroms)

    # Invalid file-like object
    s = BytesIO()
    assert pytest.raises(pybigtools.BBIReadError, pybigtools.open, s, "r")


def test_open_pathlib_path():
    path = REPO_ROOT / "bigtools/resources/test/valid.bigWig"
    with pybigtools.open(path, "r") as b:
        assert b.chroms() == {"chr17": 83_257_441}


def test_open_raw_url():
    url = "http://genome.ucsc.edu/goldenPath/help/examples/bigLollyExample2.bb"
    with pybigtools.open(url, "r") as b:
        assert b.chroms() == {'chr21': 46_709_983}


def test_open_filelike():
    # Regular file
    with open(REPO_ROOT / "bigtools/resources/test/valid.bigWig", "rb") as f:
        with pybigtools.open(f, "r") as b:
            assert b.chroms() == {"chr17": 83_257_441}

    # BytesIO
    with open(REPO_ROOT / "bigtools/resources/test/valid.bigWig", "rb") as f:
        bw_bytes = f.read()

    with BytesIO(bw_bytes) as f:
        with pybigtools.open(f, "r") as b:
            assert b.chroms() == {"chr17": 83_257_441}

    # Other file-like objects
    url = "http://genome.ucsc.edu/goldenPath/help/examples/bigWigExample.bw"
    with smart_open.open(url, "rb") as f:
        with pybigtools.open(f, "r") as b:
            assert b.chroms() == {'chr21': 48_129_895}

def test_contextmanager_exception():
    class RaisesException(Exception):
        pass
    with pytest.raises(RaisesException):
        with pybigtools.open(REPO_ROOT / "bigtools/resources/test/valid.bigWig", "r"):
            raise RaisesException()

@pytest.fixture
def bw():
    path = str(REPO_ROOT / "bigtools/resources/test/valid.bigWig")
    return pybigtools.open(path, "r")


@pytest.fixture
def bb():
    path = str(TEST_DIR / "data/bigBedExample.bb")
    return pybigtools.open(path, "r")


def test_check_filetype(bw, bb):
    assert not bw.is_bigbed
    assert bw.is_bigwig

    assert bb.is_bigbed
    assert not bb.is_bigwig


def test_info(bw):
    assert bw.info() == {
        'version': 4,
        'isCompressed': True,
        'primaryDataSize': 603305,
        'zoomLevels': 10,
        'chromCount': 1,
        'summary': {
            'basesCovered': 137894,
            'sum': 89001714.97238067,
            'mean': 645.4357330440822,
            'min': 0.006219999864697456,
            'max': 14254.0,
            'std': 751.0146556298351
        },
    }


def test_chroms(bw, bb):
    # No args => dict
    assert bw.chroms() == {"chr17": 83_257_441}
    assert bb.chroms() == {"chr21": 48_129_895}

    # Arg with chrom name => length
    assert bw.chroms("chr17") == 83_257_441
    assert bb.chroms("chr21") == 48_129_895

    # Missing chrom => KeyError
    pytest.raises(KeyError, bw.chroms, "chr11")
    pytest.raises(KeyError, bb.chroms, "chr11")


def test_zooms(bw, bb):
    # Get a list of zooms
    assert bw.zooms() == [10, 40, 160, 640, 2560, 10240, 40960, 163840, 655360, 2621440]
    assert bb.zooms() == [3911, 39110, 391100, 3911000, 39110000]


def test_autosql(bw, bb):
    # Even bigwigs have sql (a sql representing bedGraph)
    assert "bedGraph" in bw.sql()
    # We can parse the sql
    assert bw.sql(True)['name'] == 'bedGraph'

    # Unfortunately, this test bigBed doesn't actually have autosql
    assert len(bb.sql()) == 0


def test_records(bw, bb):
    # (chrom, None, None) => all records on chrom
    assert len(list(bw.records("chr17"))) == 100_000
    assert len(list(bb.records("chr21"))) == 14_810

    # (chrom, start, None) => all records from (start, <chrom_end>)
    assert len(list(bw.records("chr17", 100_000))) == 91_360
    assert len(list(bb.records("chr21", 10_000_000))) == 14_799

    # (chrom, start, end) => all records from (start, end)
    assert len(list(bw.records("chr17", 100_000, 110_000))) == 1515
    assert len(list(bb.records("chr21", 10_000_000, 20_000_000))) == 233

    # Out of bounds start/end are truncated
    x = list(bw.records("chr17", 0, 100_000))
    assert len(x) == 8641
    assert list(bw.records("chr17", -1000, 100_000)) == x
    x = list(bb.records("chr21", 0, 10_000_000))
    assert len(x) == 11
    assert list(bb.records("chr21", -1000, 10_000_000)) == x

    y = list(bw.records("chr17", 0, bw.chroms("chr17")))
    assert len(y) == 100_000
    assert list(bw.records("chr17", 0, bw.chroms("chr17") * 2)) == y
    assert next(bw.records("chr17")) == (59898, 59900, 0.06791999936103821)
    y = list(bb.records("chr21", 0, bb.chroms("chr21")))
    assert len(y) == 14810
    assert list(bb.records("chr21", 0, bb.chroms("chr21") * 2)) == y
    assert next(bb.records("chr21")) == (9434178, 9434609)

    # Fully out of bounds ranges return no records
    assert len(list(bw.records("chr17", -1000, -500))) == 0
    assert len(list(bw.records("chr17", 83_257_441, 84_000_000))) == 0

    assert len(list(bb.records("chr21", -1000, -500))) == 0
    assert len(list(bb.records("chr21", 48_129_895, 49_000_000))) == 0

    # Unknown chrom  => exception
    assert pytest.raises(KeyError, bw.records, "chr11")
    assert pytest.raises(KeyError, bb.records, "chr11")


def test_zoom_records(bw, bb):
    # (chrom, None, None) => all records on chrom
    assert len(list(bw.zoom_records(10, "chr17"))) == 13811
    assert len(list(bb.zoom_records(3911, "chr21"))) == 1676

    # (chrom, start, None) => all records from (start, <chrom_end>)
    assert len(list(bw.zoom_records(10, "chr17", 100_000))) == 10872
    assert len(list(bb.zoom_records(3911, "chr21", 10_000_000))) == 1670

    # (chrom, start, end) => all records from (start, end)
    assert len(list(bw.zoom_records(10, "chr17", 100_000, 110_000))) == 766
    assert len(list(bb.zoom_records(3911, "chr21", 10_000_000, 20_000_000))) == 154

    # Out of bounds start/end are truncated
    x = list(bw.zoom_records(10, "chr17", 0, 100_000))
    assert len(x) == 2940
    assert list(bw.zoom_records(10, "chr17", -1000, 100_000)) == x
    x = list(bb.zoom_records(3911, "chr21", 0, 10_000_000))
    assert len(x) == 6
    assert list(bb.zoom_records(3911, "chr21", -1000, 10_000_000)) == x

    y = list(bw.zoom_records(10, "chr17", 0, bw.chroms("chr17")))
    assert len(y) == 13811
    assert list(bw.zoom_records(10, "chr17", 0, bw.chroms("chr17") * 2)) == y
    y = list(bb.zoom_records(3911, "chr21", 0, bb.chroms("chr21")))
    assert len(y) == 1676
    assert list(bb.zoom_records(3911, "chr21", 0, bb.chroms("chr21") * 2)) == y

    # Fully out of bounds ranges return no records
    assert len(list(bw.zoom_records(10, "chr17", -1000, -500))) == 0
    assert len(list(bw.zoom_records(10, "chr17", 83_257_441, 84_000_000))) == 0
    assert len(list(bb.zoom_records(3911, "chr21", -1000, -500))) == 0
    assert len(list(bb.zoom_records(3911, "chr21", 48_129_895, 49_000_000))) == 0

    assert next(bw.zoom_records(10, "chr17", 0, 100_000)) == (
        59898,
        59908,
        {
            "total_items": 0,
            "bases_covered": 10,
            "min_val": 0.06791999936103821,
            "max_val": 0.16627000272274017,
            "sum": 1.4660000801086426,
            "sum_squares": 0.2303919643163681,
        },
    )
    assert next(bw.zoom_records(160, "chr17", 0, 100000)) == (
        59898,
        60058,
        {
            "total_items": 0,
            "bases_covered": 160,
            "min_val": 0.06791999936103821,
            "max_val": 0.8688300251960754,
            "sum": 101.3516616821289,
            "sum_squares": 80.17473602294922,
        },
    )

    # Unknown zoom  => exception
    assert pytest.raises(KeyError, bw.zoom_records, 0, "chr17")
    assert pytest.raises(KeyError, bb.zoom_records, 0, "chr21")

    # Unknown chrom  => exception
    assert pytest.raises(KeyError, bw.zoom_records, 10, "chr11")
    assert pytest.raises(KeyError, bb.zoom_records, 3911, "chr11")


def test_values(bw, bb):
    assert len(bw.values("chr17", 100_000, 110_000)) == 10_000
    assert len(bb.values("chr21", 10_148_000, 10_158_000)) == 10_000

    assert len(bw.values("chr17", 100000, 110000, 10)) == 10
    assert len(bb.values("chr21", 10_148_000, 10_158_000, 10)) == 10

    assert bw.values("chr17", 100000, 110000, 10)[0] == 0.44886381671868925
    assert bb.values("chr21", 10_148_000, 10_158_000, 10)[0] == 1.0

    assert bw.values("chr17", 100000, 110000, 10, "max")[0] == 1.1978399753570557
    assert bb.values("chr21", 10_148_000, 10_158_000, 10, "max")[0] == 1.0

    assert bw.values("chr17", 100000, 110000, 10, "min")[0] == 0.05403999984264374
    assert bb.values("chr21", 10_148_000, 10_158_000, 10, "min")[0] == 0.0

    assert (
        bw.values("chr17", 100000, 110000, 10, "mean", exact=True)[0]
        == 0.4542629980980206
    )
    assert (
        bb.values("chr21", 10_148_000, 10_158_000, 10, "mean", exact=True)[0]
        == 1.0
    )

    assert list(bw.values("chr17", 59890, 59900, 10, "mean", exact=True)) == [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.06791999936103821,
        0.06791999936103821,
    ]

    assert list(
        bw.values("chr17", 59890, 59900, 10, "mean", exact=True, missing=-1.0)
    ) == [
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        -1.0,
        0.06791999936103821,
        0.06791999936103821,
    ]

    x = bw.values("chr17", -10, 10, 20, "mean", exact=True, missing=0.0)
    assert math.isnan(x[0])
    assert not math.isnan(x[19])
    x = bb.values("chr21", -10, 10, 20, "mean", exact=True, missing=0.0)
    assert math.isnan(x[0])
    assert not math.isnan(x[19])

    x = bw.values("chr17", -10, 10, 20, "mean", exact=True, missing=0.0, oob=0.0)
    assert x[0] == 0.0
    x = bb.values("chr21", -10, 10, 20, "mean", exact=True, missing=0.0, oob=0.0)
    assert x[0] == 0.0

    # The returned array is the same as the one passed, so both show the same values
    arr = np.zeros(20)
    ret_arr = bw.values(
        "chr17", -10, 10, 20, "mean", exact=True, missing=0.0, oob=np.nan, arr=arr
    )
    assert math.isnan(arr[0])
    assert arr[19] == 0.0
    assert math.isnan(ret_arr[0])
    assert ret_arr[19] == 0.0
    assert np.array_equal(arr, ret_arr, equal_nan=True)

    ret_arr = bb.values(
        "chr21", -10, 10, 20, "mean", exact=True, missing=0.0, oob=np.nan, arr=arr
    )
    assert math.isnan(arr[0])
    assert arr[19] == 0.0
    assert math.isnan(ret_arr[0])
    assert ret_arr[19] == 0.0
    assert np.array_equal(arr, ret_arr, equal_nan=True)

    # Some differences in estimates between pybigtools to other libs
    # Namely, bigtools calculates estimates by taking the
    # sum of nanmeans over covered bases (summary.sum/summary.covered_bases) and dividing by covered bases (overlap between zooom and bin)
    # So, including these as cases where the calculated value is different
    vals = bw.values("chr17", 85525, 85730, bins=2, exact=False)
    assert list(vals) == [0.15392776003070907, 2.728891665264241]
    vals = bw.values("chr17", 85525, 85730, bins=2, exact=True)
    assert list(vals) == [0.06770934917680595, 2.4864424403431347]
    vals = bw.values("chr17", 59900, 60105, bins=2, exact=False)
    assert list(vals) == [0.5358060553962108, 0.5513471488813751]
    vals = bw.values("chr17", 59900, 60105, bins=2, exact=True)
    assert list(vals) == [0.5362001863472602, 0.5527710799959679]

    vals = bb.values("chr21", 14_760_000, 14_800_000, bins=1, exact=False)
    assert list(vals) == [1.2572170068028603]
    vals = bb.values("chr21", 14_760_000, 14_800_000, bins=1, exact=True)
    assert list(vals) == [1.3408662900188324]
