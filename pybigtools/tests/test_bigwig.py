import math
import os
import pathlib
import urllib.request

import pybigtools

TEST_DIR = pathlib.Path(__file__).parent


def retrieve_encode_file(accession, filetype):
    file_url = f"https://www.encodeproject.org/files/{accession}/@@download/{accession}.{filetype}"
    path = TEST_DIR / f"data/{accession}.{filetype}"
    if not path.exists():
        urllib.request.urlretrieve(file_url, path)
    return path


def test_bigwig_read():
    path = retrieve_encode_file("ENCFF667CZO", "bigWig")
    b = pybigtools.open(str(path))

    i = b.records("chr1")
    n = next(i)
    assert n[0] == 10495
    assert n[1] == 10545
    assert math.isclose(n[2], 0.01591, abs_tol=0.00001)

    i = b.records("chr12")
    c = 0
    for _ in i:
        c += 1
    assert c == 1028747


def test_bigwig_write(tmpdir):
    chroms = ["chr1", "chr2", "chr3"]
    clengths = {"chr1": 10000, "chr2": 8000, "chr3": 6000}

    def genintervals():
        import random

        for chrom in chroms:
            clength = clengths[chrom]
            current = random.randint(0, 300)
            start = current
            while True:
                length = random.randint(1, 200)
                end = start + length
                if end > clength:
                    break
                value = round(random.random(), 5)
                yield (chrom, start, end, value)
                start = end + random.randint(20, 50)

    intervals = list(genintervals())
    b = pybigtools.open(os.path.join(tmpdir, "test.bigWig"), "w")
    b.write(clengths, iter(intervals))
    # If didn't need to test, could also be done like
    # b.write(clengths, genintervals())

    b = pybigtools.open(os.path.join(tmpdir, "test.bigWig"))
    i = []
    for chrom in chroms:
        i.extend(list(b.records(chrom)))
    c = 0
    for _ in i:
        c += 1
    assert c == len(intervals)
    for a, b in zip(i, intervals):
        assert a[0] == b[1]
        assert a[1] == b[2]
        assert math.isclose(a[2], b[3], abs_tol=0.00001)

def test_bigbed_write(tmpdir):
    f = pybigtools.open(os.path.join(tmpdir, "test.bigBed"), "w")
    f.write(
        {"chr1": 1000, "chr2": 1000}, 
        [
            ("chr1", 0, 100, "foo\tbar\t3"), 
            ("chr2", 100, 200, "xxx\tyyy\t1.0"),
            ("chr2", 200, 300, "zzz\twww\t1.0")
        ],
        autosql="""\
table bed3
"Simple bed"
(
    string chrom;        "Reference sequence chromosome or scaffold"
    uint   chromStart;   "Start position in chromosome"
    uint   chromEnd;     "End position in chromosome"
)
        """
    )
    f.close()

    f = pybigtools.open(os.path.join(tmpdir, "test.bigBed"))
    records = list(f.records("chr2"))
    assert records[0][2] == 'xxx'
    sql = f.sql()
    assert sql.startswith("table")

# TODO: bigWigAverageOverBed
# TODO: bigWigMerge
