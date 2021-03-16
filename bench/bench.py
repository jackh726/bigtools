import os
import re
import resource
import struct
import subprocess
import sys
import urllib.request

from monitorps import monitor

timeregex = re.compile(r'(\d.*)user (\d.*)system (\d.*)elapsed')
ucsctoolspath = None
bigtoolspath = None
reps = 3

def download_test_data():
    print("Downloading test data")
    if not os.path.exists("./workdir"):
        os.makedirs("./workdir")
    if not os.path.exists("./workdir/output"):
        os.makedirs("./workdir/output")
    def download(url, name):
        filename = name.replace('.gz', '')
        if not os.path.exists(f'./workdir/{filename}'):
            print(f"Downloading {name}")
            urllib.request.urlretrieve(url, f'./workdir/{name}')
            if name.endswith(".gz"):
                subprocess.check_call(['gunzip', './workdir/ENCFF646AZP.bed.gz'])

    download('https://www.encodeproject.org/files/ENCFF937MNZ/@@download/ENCFF937MNZ.bigWig', 'ENCFF937MNZ.bigWig')
    download('https://www.encodeproject.org/files/ENCFF447DHW/@@download/ENCFF447DHW.bigWig', 'ENCFF447DHW.bigWig')
    download('https://www.encodeproject.org/files/ENCFF841DHZ/@@download/ENCFF841DHZ.bigWig', 'ENCFF841DHZ.bigWig')
    download('https://www.encodeproject.org/files/ENCFF518WII/@@download/ENCFF518WII.bigWig', 'ENCFF518WII.bigWig')
    download('https://www.encodeproject.org/files/ENCFF646AZP/@@download/ENCFF646AZP.bed.gz', 'ENCFF646AZP.bed.gz')
    download('https://www.encodeproject.org/files/ENCFF076CIO/@@download/ENCFF076CIO.bed.gz', 'ENCFF076CIO.bed.gz')
    download('https://raw.githubusercontent.com/ENCODE-DCC/encValData/562ab5bf03deff9bb5340991fd5c844162b82914/GRCh38/GRCh38_EBV.chrom.sizes', 'hg38.chrom.sizes')
    print("Done downloading test data")


def time(exeargs_all, bench, program, rep):
    total_seconds = 0
    for i, exeargs in enumerate(exeargs_all):
        exeargs.insert(0, '/usr/bin/time')
        exeargs = " ".join(exeargs)
        print(exeargs)
        logfile = f'./workdir/output/{bench}_{program}_{rep}_{i}.log'
        with subprocess.Popen(exeargs, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True) as process:
            monitor(process.pid, logfile, 0.1)
            out, err = process.communicate()
            err = err.decode('utf-8')
            m = timeregex.search(err)
            elapsed_split = m.group(3).split(':')
            # Seconds
            total_seconds += float(elapsed_split.pop())
            if len(elapsed_split) > 0:
                # Minutes
                total_seconds += int(elapsed_split.pop()) * 60
            if len(elapsed_split) > 0:
                # Hours
                total_seconds += int(elapsed_split.pop()) * 60 * 60

    return total_seconds

def compare(comp, bench, ucsc, bigtools_mt=None, bigtools_st=None):
    global reps
    print('Benchmarking {}'.format(bench))
    for i in range(0, reps):
        print(f"Rep {i+1}...")
        ucsctime = time(ucsc, bench, 'ucsc', i)
        print("ucsc: {}".format(round(ucsctime,3)))
        comp.write(f"{bench}\tucsc\t{round(ucsctime,3)}\n")
        if bigtools_mt is not None:
            bigtoolstime = time(bigtools_mt, bench, 'bigtools-mt', i)
            print("bigtools-mt: {}".format(round(bigtoolstime,3)))
            comp.write(f"{bench}\tbigtools-mt\t{round(bigtoolstime,3)}\n")
        if bigtools_st is not None:
            bigtoolssttime = time(bigtools_st, bench, 'bigtools-st', i)
            print("bigtools-st: {}".format(round(bigtoolssttime,3)))
            comp.write(f"{bench}\tbigtools-st\t{round(bigtoolssttime,3)}\n")

def bigwigaverageoverbed(comp):
    global ucsctoolspath
    global bigtoolspath
    # For ucsc, we have to convert narrowPeak to bed first, including adding a unique name
    ucsc = [
        ['cat ./workdir/ENCFF646AZP.bed | cut -f1-3 | awk -v OFS=\'\\t\' \'{print $1,$2,$3, NR}\' > ./workdir/ENCFF646AZP_cut.bed'],
        ['{}/bigWigAverageOverBed'.format(ucsctoolspath), './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF646AZP_cut.bed', './workdir/test_out_ucsc.bed']
        ]
    bigtools_st = [['{}/bigwigaverageoverbed'.format(bigtoolspath), './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF646AZP.bed', './workdir/test_out_bigtools.bed']]
    compare(comp, 'bigwigaverageoverbed', ucsc, None, bigtools_st)

def bigwigmerge(comp):
    global ucsctoolspath
    global bigtoolspath
    ucsc = [['{}/bigWigMerge'.format(ucsctoolspath), './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF447DHW.bigWig', './workdir/test_out_ucsc.bedGraph']]
    bigtools = [['{}/bigwigmerge'.format(bigtoolspath), './workdir/test_out_bigtools.bedGraph', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig']]
    compare(comp, 'bigwigmerge_bedgraph', ucsc, bigtools)
    ucsc = [
        ['{}/bigWigMerge'.format(ucsctoolspath), './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF447DHW.bigWig', './workdir/test_out_ucsc.bedGraph'],
        ['{}/bedGraphToBigWig'.format(ucsctoolspath), './workdir/test_out_ucsc.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_ucsc.bigWig']
        ]
    bigtools = [['{}/bigwigmerge'.format(bigtoolspath), './workdir/test_out_bigtools.bigWig', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig']]
    bigtools_st = [['{}/bigwigmerge'.format(bigtoolspath), './workdir/test_out_bigtools.bigWig', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig', '-t 1']]
    compare(comp, 'bigwigmerge_bigwig', ucsc, bigtools, bigtools_st)

def bedgraphtobigwig(comp):
    global ucsctoolspath
    global bigtoolspath
    # Need to generate bedGraph first, just use ucsc since it's what we're comparing against
    if not os.path.exists('./workdir/ENCFF518WII.bedGraph'):
        process = subprocess.check_call('{}/bigWigToBedGraph ./workdir/ENCFF518WII.bigWig ./workdir/ENCFF518WII.bedGraph'.format(ucsctoolspath), shell=True)
    if not os.path.exists('./workdir/ENCFF841DHZ.bedGraph'):
        process = subprocess.check_call('{}/bigWigToBedGraph ./workdir/ENCFF841DHZ.bigWig ./workdir/ENCFF841DHZ.bedGraph'.format(ucsctoolspath), shell=True)
    ucsc = [['{}/bedGraphToBigWig'.format(ucsctoolspath), './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_ucsc.bigWig']]
    bigtools = [['{}/bedgraphtobigwig'.format(bigtoolspath), './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig']]
    bigtools_st = [['{}/bedgraphtobigwig'.format(bigtoolspath), './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 1']]
    compare(comp, 'bedgraphtobigwig_small', ucsc, bigtools, bigtools_st)
    ucsc = [['{}/bedGraphToBigWig'.format(ucsctoolspath), './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_ucsc.bigWig']]
    bigtools = [['{}/bedgraphtobigwig'.format(bigtoolspath), './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig']]
    bigtools_st = [['{}/bedgraphtobigwig'.format(bigtoolspath), './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 1']]
    compare(comp, 'bedgraphtobigwig_medium', ucsc, bigtools, bigtools_st)

def bedtobigbed(comp):
    global ucsctoolspath
    global bigtoolspath
    if not os.path.exists('./workdir/ENCFF076CIO.sorted.bed'):
        process = subprocess.check_call('cat ./workdir/ENCFF076CIO.bed | cut -f1-3 | sort -k1,1 -k2,2n > ./workdir/ENCFF076CIO.sorted.bed', shell=True)
    ucsc = [['{}/bedGraphToBigWig'.format(ucsctoolspath), './workdir/ENCFF076CIO.sorted.bed', './workdir/hg38.chrom.sizes', './workdir/test_out_ucsc.bigBed']]
    bigtools = [['{}/bedtobigbed'.format(bigtoolspath), './workdir/ENCFF076CIO.sorted.bed', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigBed']]
    bigtools_st = [['{}/bedtobigbed'.format(bigtoolspath), './workdir/ENCFF076CIO.sorted.bed', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigBed', '-t 1']]
    compare(comp, 'bedtobigbed', ucsc, bigtools, bigtools_st)

def bigwigtobedgraph(comp):
    global ucsctoolspath
    global bigtoolspath
    ucsc = [['{}/bigWigToBedGraph'.format(ucsctoolspath), './workdir/ENCFF841DHZ.bigWig', './workdir/test_out_ucsc.bedGraph']]
    bigtools = [['{}/bigwigtobedgraph'.format(bigtoolspath), './workdir/ENCFF841DHZ.bigWig', './workdir/test_out_bigtools.bedGraph']]
    bigtools_st = [['{}/bigwigtobedgraph'.format(bigtoolspath), './workdir/ENCFF841DHZ.bigWig', './workdir/test_out_bigtools.bedGraph', '-t 1']]
    compare(comp, 'bigwigtobedgraph', ucsc, bigtools, bigtools_st)


def main():
    global ucsctoolspath
    global bigtoolspath
    global reps
    ucsctoolspath = sys.argv[1] if len(sys.argv) > 1 else "/bin"
    bigtoolspath = sys.argv[2] if len(sys.argv) > 2 else "/usr/local/bin"
    bench = sys.argv[3] if len(sys.argv) > 3 else None
    reps = 3
    download_test_data()
    with open("./workdir/output/comparison.txt", "w") as comp:
        if bench is None or bench == 'bigwigaverageoverbed':
            bigwigaverageoverbed(comp)
        if bench is None or bench == 'bigwigmerge':
            bigwigmerge(comp)
        if bench is None or bench == 'bedgraphtobigwig':
            bedgraphtobigwig(comp)
        if bench is None or bench == 'bedtobigbed':
            bedtobigbed(comp)
        if bench is None or bench == 'bigwigtobedgraph':
            bigwigtobedgraph(comp)

if __name__ == '__main__':
    main()

# To run (example): python3 bench/bench.py ~/temp /mnt/c/Users/jackh/Documents/Git/rust/bigtools/target/release
