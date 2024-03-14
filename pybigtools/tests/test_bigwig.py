import os

def install_dependencies():
    import urllib.request
    if not os.path.exists('./ENCFF667CZO.bigWig'):
        print("Downloading test bigWig")
        urllib.request.urlretrieve('https://www.encodeproject.org/files/ENCFF667CZO/@@download/ENCFF667CZO.bigWig', './ENCFF667CZO.bigWig')

print("Testing")
# Imports for test
import math

# Test
import pybigtools

def test_bigwig_write_read(tmp_path):
    install_dependencies()
    # bigWig read
    b = pybigtools.open('./ENCFF667CZO.bigWig')

    i = b.intervals("chr1")
    n = next(i)
    assert n[0] == 10495
    assert n[1] == 10545
    assert math.isclose(n[2], 0.01591, abs_tol=0.00001)

    i = b.intervals("chr12")
    c = 0
    for _ in i:
        c += 1
    assert c == 1028747

    # bigWig write
    chroms = ['chr1', 'chr2', 'chr3']
    clengths = { 'chr1': 10000, 'chr2': 8000, 'chr3': 6000 }
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
    b = pybigtools.open('./test.bigWig', 'w')
    b.write(clengths, iter(intervals))
    # If didn't need to test, could also be done like
    #b.write(clengths, genintervals())

    b = pybigtools.open('./test.bigWig')
    i = []
    for chrom in chroms:
        i.extend(list(b.intervals(chrom)))
    c = 0
    for _ in i:
        c += 1
    assert c == len(intervals)
    for a,b in zip(i, intervals):
        assert a[0] == b[1]
        assert a[1] == b[2]
        assert math.isclose(a[2], b[3], abs_tol=0.00001)


# TODO: bigWigAverageOverBed
# TODO: bigWigMerge
