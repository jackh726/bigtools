import pybigtools
import math
import numpy as np

def test_api():
    b = pybigtools.open('../bigtools/resources/test/valid.bigWig', 'r')

    # Check type of file
    assert b.is_bigbed == False
    assert b.is_bigwig == True

    # No args => dict
    assert b.chroms() == {'chr17': 83257441}
    # Arg with chrom name => length
    assert b.chroms('chr17') == 83257441
    # Missing chrom => None
    assert b.chroms('chr11') == None

    # Get a list of zooms
    assert b.zooms() == [10, 40, 160, 640, 2560, 10240, 40960, 163840, 655360, 2621440]

    # Even bigwigs have sql (a sql representing bedGraph)
    assert "bedGraph" in b.sql()
    # We can parse the sql
    assert b.sql(True)['name'] == 'bedGraph'

    # (chrom, None, None) => all records on chrom
    assert len(list(b.records('chr17'))) == 100000
    # (chrom, start, None) => all records from (start, <chrom_end>)
    assert len(list(b.records('chr17', 100000))) == 91360
    # (chrom, start, end) => all records from (start, end)
    assert len(list(b.records('chr17', 100000, 110000))) == 1515
    # Out of bounds start/end are truncated
    assert len(list(b.records('chr17', -1000, 100000))) == 8641
    assert len(list(b.records('chr17', -1000, -500))) == 0
    assert len(list(b.records('chr17', 0, 84000000))) == 100000
    assert next(b.records('chr17')) == (59898, 59900, 0.06791999936103821)
    # Unknown chrom  => exception
    try:
        b.zoom_records('chr11')
    except:
        pass
    else:
        assert False

    # (chrom, None, None) => all records on chrom
    assert len(list(b.zoom_records(10, 'chr17'))) == 13811
    # (chrom, start, None) => all records from (start, <chrom_end>)
    assert len(list(b.zoom_records(10, 'chr17', 100000))) == 10872
    # (chrom, start, end) => all records from (start, end)
    assert len(list(b.zoom_records(10, 'chr17', 100000, 110000))) == 766
    # Out of bounds start/end are truncated
    assert len(list(b.zoom_records(10, 'chr17', -1000, 100000))) == 2940
    assert len(list(b.zoom_records(10, 'chr17', -1000, -500))) == 0
    assert len(list(b.zoom_records(10, 'chr17', 0, 84000000))) == 13811
    assert next(b.zoom_records(10, 'chr17', 0, 100000)) == (59898, 59908, {'total_items': 0, 'bases_covered': 10, 'min_val': 0.06791999936103821, 'max_val': 0.16627000272274017, 'sum': 1.4660000801086426, 'sum_squares': 0.2303919643163681})
    assert next(b.zoom_records(160, 'chr17', 0, 100000)) == (59898, 60058, {'total_items': 0, 'bases_covered': 160, 'min_val': 0.06791999936103821, 'max_val': 0.8688300251960754, 'sum': 101.3516616821289, 'sum_squares': 80.17473602294922})
    # Unknown zoom  => exception
    try:
        b.zoom_records(0, 'chr17')
    except:
        pass
    else:
        assert False
    # Unknown chrom  => exception
    try:
        b.zoom_records(10, 'chr11')
    except:
        pass
    else:
        assert False

    assert len(list(b.values('chr17', 100000, 110000))) == 10000
    assert len(list(b.values('chr17', 100000, 110000, 10))) == 10
    assert b.values('chr17', 100000, 110000, 10)[0] == 0.37435242314338685
    assert b.values('chr17', 100000, 110000, 10, 'max')[0] == 1.1978399753570557
    assert b.values('chr17', 100000, 110000, 10, 'min')[0] == 0.05403999984264374
    assert b.values('chr17', 100000, 110000, 10, 'mean',  True)[0] == 0.37885534041374924
    assert list(b.values('chr17', 59890, 59900, 10, 'mean', True)) == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.06791999936103821, 0.06791999936103821]
    assert list(b.values('chr17', 59890, 59900, 10, 'mean', True, -1.0)) == [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.06791999936103821, 0.06791999936103821]
    assert math.isnan(b.values('chr17', -10, 10, 20, 'mean', True, 0.0)[0])
    assert not math.isnan(b.values('chr17', -10, 10, 20, 'mean', True, 0.0)[19])
    assert b.values('chr17', -10, 10, 20, 'mean', True, 0.0, 0.0)[0] == 0.0
    arr = np.zeros(20)
    ret_arr = b.values('chr17', -10, 10, 20, 'mean', True, 0.0, np.nan, arr)
    # The returned array is the same as the one passed, so both show the same values
    assert math.isnan(arr[0])
    assert arr[19] == 0.0
    assert math.isnan(ret_arr[0])
    assert ret_arr[19] == 0.0

    b.close()
    # Closing means file is not usable
    try:
        b.chroms()
    except:
        pass
    else:
        assert False

    b = pybigtools.open('../bigtools/resources/test/valid.bigWig', 'r')
    # Files are closed when exiting a context manager
    with b:
      pass
    try:
        b.chroms()
    except:
        pass
    else:
        assert False

