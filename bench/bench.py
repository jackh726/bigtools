import os
import re
import subprocess
import sys
import urllib.request

from monitorps import monitor

# To run (example): python3 bench/bench.py ./workdir/kent ./target/release
# Or: time python3 bench/bench.py ./workdir/kent ./target/release > workdir/bench_out.log 2>&1

timeregex = re.compile(r'(\d.*)user (\d.*)system (\d.*)elapsed (\d.*)%CPU \((\d*)avgtext\+(\d*)avgdata (\d*)maxresident\)k')
ucsctoolspath = None
bigtoolspath = None
reps = 3

def download(url, name):
    filename = name.replace('.gz', '')
    if not os.path.exists(f'./workdir/{filename}'):
        print(f"Downloading {name}")
        urllib.request.urlretrieve(url, f'./workdir/{name}')
        if name.endswith(".gz"):
            print(f"Ungzipping {name}")
            subprocess.check_call(['gunzip', f'./workdir/{name}'])

def download_test_data():
    print("Downloading test data")
    os.makedirs("./workdir/output", exist_ok=True)

    download('https://www.encodeproject.org/files/ENCFF937MNZ/@@download/ENCFF937MNZ.bigWig', 'ENCFF937MNZ.bigWig')
    download('https://www.encodeproject.org/files/ENCFF447DHW/@@download/ENCFF447DHW.bigWig', 'ENCFF447DHW.bigWig')
    download('https://www.encodeproject.org/files/ENCFF841DHZ/@@download/ENCFF841DHZ.bigWig', 'ENCFF841DHZ.bigWig')
    download('https://www.encodeproject.org/files/ENCFF518WII/@@download/ENCFF518WII.bigWig', 'ENCFF518WII.bigWig')
    download('https://www.encodeproject.org/files/ENCFF646AZP/@@download/ENCFF646AZP.bed.gz', 'ENCFF646AZP.bed.gz')
    download('https://www.encodeproject.org/files/ENCFF076CIO/@@download/ENCFF076CIO.bed.gz', 'ENCFF076CIO.bed.gz')
    download('https://raw.githubusercontent.com/ENCODE-DCC/encValData/562ab5bf03deff9bb5340991fd5c844162b82914/GRCh38/GRCh38_EBV.chrom.sizes', 'hg38.chrom.sizes')

def download_kent_tools():
    print("Downloading kent tools")
    os.makedirs('./workdir/kent', exist_ok=True)

    download('http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed', 'kent/bedToBigBed')
    download('http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig', 'kent/bedGraphToBigWig')
    download('http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigAverageOverBed', 'kent/bigWigAverageOverBed')
    download('http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigMerge', 'kent/bigWigMerge')
    download('http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph', 'kent/bigWigToBedGraph')

def time(exeargs_all, bench, program, rep):
    total_seconds = 0
    for i, exeargs in enumerate(exeargs_all):
        exeargs = '/usr/bin/time ' + " ".join(exeargs)
        print(exeargs)
        logfile = f'./workdir/output/{bench}_{program}_{rep}_{i}.log'
        with subprocess.Popen(exeargs, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True) as process:
            monitor(process.pid, logfile, 0.25)
            out, err = process.communicate()
            out = out.decode('utf-8')
            err = err.decode('utf-8')
            #print(out, err)
            matched_regex = timeregex.search(err)
            cpu = int(matched_regex.group(4))
            max_resident = int(matched_regex.group(7))
            elapsed_split = matched_regex.group(3).split(':')
            # Seconds
            total_seconds += float(elapsed_split.pop())
            if len(elapsed_split) > 0:
                # Minutes
                total_seconds += int(elapsed_split.pop()) * 60
            if len(elapsed_split) > 0:
                # Hours
                total_seconds += int(elapsed_split.pop()) * 60 * 60

    return total_seconds,cpu,max_resident

def compare(comp, bench, benchmarks):
    global reps
    print('Benchmarking {}'.format(bench))
    def split_round_format(time_res):
        total_seconds,cpu,max_resident = time_res
        total_seconds = round(total_seconds, 3)
        return f"{total_seconds}\t{cpu}\t{max_resident}"
    for i in range(0, reps):
        print(f"Rep {i+1}...")
        for benchmark,command in benchmarks.items():
            bench_time = split_round_format(time(command, bench, benchmark, i))
            print(f"{benchmark}: {bench_time}", flush=True)
            comp.write(f"{bench}\t{benchmark}\t{bench_time}\n")

def bigwigaverageoverbed(comp):
    global ucsctoolspath
    global bigtoolspath
    # For ucsc, we have to convert narrowPeak to bed first, including adding a unique name
    if not os.path.exists('./workdir/ENCFF646AZP_cut.bed'):
        process = subprocess.check_call('cat ./workdir/ENCFF646AZP.bed | cut -f1-3 | awk -v OFS=\'\\t\' \'{print $1,$2,$3, NR}\' > ./workdir/ENCFF646AZP_cut.bed', shell=True)
    benchmarks = {
        'ucsc':             [['{}/bigWigAverageOverBed'.format(ucsctoolspath), './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF646AZP_cut.bed', '/dev/null']],
        'bigtools_1thread': [['{}/bigwigaverageoverbed'.format(bigtoolspath),  './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF646AZP_cut.bed', '/dev/null', '-t 1']],
        'bigtools_2thread': [['{}/bigwigaverageoverbed'.format(bigtoolspath),  './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF646AZP_cut.bed', '/dev/null', '-t 2']],
        'bigtools_4thread': [['{}/bigwigaverageoverbed'.format(bigtoolspath),  './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF646AZP_cut.bed', '/dev/null', '-t 4']],
        'bigtools_6thread': [['{}/bigwigaverageoverbed'.format(bigtoolspath),  './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF646AZP_cut.bed', '/dev/null', '-t 6']],
        'bigtools_8thread': [['{}/bigwigaverageoverbed'.format(bigtoolspath),  './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF646AZP_cut.bed', '/dev/null', '-t 8']],
    }
    compare(comp, 'bigwigaverageoverbed', benchmarks)

def bigwigaverageoverbed_long(comp):
    global ucsctoolspath
    global bigtoolspath
    # For ucsc, we have to convert narrowPeak to bed first, including adding a unique name
    if not os.path.exists('./workdir/ENCFF076CIO_cut_sample.bed'):
        process = subprocess.check_call(f'{bigtoolspath}/bigtools chromintersect -a ./workdir/ENCFF076CIO.bed -b ./workdir/ENCFF937MNZ.bigWig -o -' + ' | cut -f1-3 | awk -v OFS=\'\\t\' \'{print $1,$2,$3, NR}\' | shuf --random-source=./workdir/ENCFF076CIO.bed | head -1000000 | sort -k1,1 -k2,2n > ./workdir/ENCFF076CIO_cut_sample.bed', shell=True)
    benchmarks = {
        'ucsc':             [['{}/bigWigAverageOverBed'.format(ucsctoolspath), './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF076CIO_cut_sample.bed', '/dev/null']],
        'bigtools_1thread': [['{}/bigwigaverageoverbed'.format(bigtoolspath),  './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF076CIO_cut_sample.bed', '/dev/null', '-t 1']],
        'bigtools_2thread': [['{}/bigwigaverageoverbed'.format(bigtoolspath),  './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF076CIO_cut_sample.bed', '/dev/null', '-t 2']],
        'bigtools_4thread': [['{}/bigwigaverageoverbed'.format(bigtoolspath),  './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF076CIO_cut_sample.bed', '/dev/null', '-t 4']],
        'bigtools_6thread': [['{}/bigwigaverageoverbed'.format(bigtoolspath),  './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF076CIO_cut_sample.bed', '/dev/null', '-t 6']],
        'bigtools_8thread': [['{}/bigwigaverageoverbed'.format(bigtoolspath),  './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF076CIO_cut_sample.bed', '/dev/null', '-t 8']],
    }
    compare(comp, 'bigwigaverageoverbed_long', benchmarks)

def bigwigmerge_bigwig(comp):
    global ucsctoolspath
    global bigtoolspath
    benchmarks = {
        'ucsc': [
            ['{}/bigWigMerge'.format(ucsctoolspath), './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF447DHW.bigWig', './workdir/test_out_ucsc.bedGraph'],
            ['{}/bedGraphToBigWig'.format(ucsctoolspath), './workdir/test_out_ucsc.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_ucsc.bigWig'],
        ],
        'bigtools_1thread': [['{}/bigwigmerge'.format(bigtoolspath), './workdir/test_out_bigtools.bigWig', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig', '--output-type', 'bigwig', '-t 1']],
        'bigtools_2thread': [['{}/bigwigmerge'.format(bigtoolspath), './workdir/test_out_bigtools.bigWig', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig', '--output-type', 'bigwig', '-t 2']],
        'bigtools_4thread': [['{}/bigwigmerge'.format(bigtoolspath), './workdir/test_out_bigtools.bigWig', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig', '--output-type', 'bigwig', '-t 4']],
        'bigtools_6thread': [['{}/bigwigmerge'.format(bigtoolspath), './workdir/test_out_bigtools.bigWig', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig', '--output-type', 'bigwig', '-t 6']],
        'bigtools_8thread': [['{}/bigwigmerge'.format(bigtoolspath), './workdir/test_out_bigtools.bigWig', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig', '--output-type', 'bigwig', '-t 8']],
    }
    compare(comp, 'bigwigmerge_bigwig', benchmarks)
    
def bigwigmerge_bedgraph(comp):
    global ucsctoolspath
    global bigtoolspath
    benchmarks = {
        'ucsc': [['{}/bigWigMerge'.format(ucsctoolspath), './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF447DHW.bigWig', '/dev/null']],
        'bigtools_1thread': [['{}/bigwigmerge'.format(bigtoolspath), '/dev/null', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig', '--output-type', 'bedgraph', '-t 1']],
        #'bigtools_2thread': [['{}/bigwigmerge'.format(bigtoolspath), '/dev/null', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig', '--output-type', 'bedgraph', '-t 2']],
        #'bigtools_4thread': [['{}/bigwigmerge'.format(bigtoolspath), '/dev/null', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig', '--output-type', 'bedgraph', '-t 4']],
        #'bigtools_6thread': [['{}/bigwigmerge'.format(bigtoolspath), '/dev/null', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig', '--output-type', 'bedgraph', '-t 6']],
        #'bigtools_8thread': [['{}/bigwigmerge'.format(bigtoolspath), '/dev/null', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig', '--output-type', 'bedgraph', '-t 8']],
    }
    compare(comp, 'bigwigmerge_bedgraph', benchmarks)

def bigwigmerge(comp):
    bigwigmerge_bedgraph(comp)
    #bigwigmerge_bigwig(comp)


def bedgraphtobigwig(comp):
    global ucsctoolspath
    global bigtoolspath
    # Need to generate bedGraph first, just use ucsc since it's what we're comparing against
    if not os.path.exists('./workdir/ENCFF518WII.bedGraph'):
        process = subprocess.check_call('{}/bigWigToBedGraph ./workdir/ENCFF518WII.bigWig ./workdir/ENCFF518WII.bedGraph'.format(ucsctoolspath), shell=True)
    if not os.path.exists('./workdir/ENCFF841DHZ.bedGraph'):
        process = subprocess.check_call('{}/bigWigToBedGraph ./workdir/ENCFF841DHZ.bigWig ./workdir/ENCFF841DHZ.bedGraph'.format(ucsctoolspath), shell=True)
    benchmarks = {
        'ucsc':                                 [['{}/bedGraphToBigWig'.format(ucsctoolspath), './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_ucsc.bigWig']],
        'bigtools_1thread_multipass_serial':    [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 1', '-p no']],
        'bigtools_2thread_multipass_serial':    [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 2', '-p no']],
        'bigtools_4thread_multipass_serial':    [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 4', '-p no']],
        'bigtools_6thread_multipass_serial':    [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 6', '-p no']],
        'bigtools_8thread_multipass_serial':    [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 8', '-p no']],
        'bigtools_2thread_multipass_parallel':  [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 2', '-p yes']],
        'bigtools_4thread_multipass_parallel':  [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 4', '-p yes']],
        'bigtools_6thread_multipass_parallel':  [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 6', '-p yes']],
        'bigtools_8thread_multipass_parallel':  [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 8', '-p yes']],
        'bigtools_1thread_singlepass_serial':   [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 1', '-p no', '--single-pass']],
        'bigtools_2thread_singlepass_serial':   [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 2', '-p no', '--single-pass']],
        'bigtools_4thread_singlepass_serial':   [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 4', '-p no', '--single-pass']],
        'bigtools_6thread_singlepass_serial':   [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 6', '-p no', '--single-pass']],
        'bigtools_8thread_singlepass_serial':   [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 8', '-p no', '--single-pass']],
        'bigtools_1thread_singlepass_parallel': [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 1', '-p yes', '--single-pass']],
        'bigtools_2thread_singlepass_parallel': [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 2', '-p yes', '--single-pass']],
        'bigtools_4thread_singlepass_parallel': [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 4', '-p yes', '--single-pass']],
        'bigtools_6thread_singlepass_parallel': [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 6', '-p yes', '--single-pass']],
        'bigtools_8thread_singlepass_parallel': [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 8', '-p yes', '--single-pass']],
    }
    compare(comp, 'bedgraphtobigwig_small', benchmarks)
    benchmarks = {
        'ucsc':                                 [['{}/bedGraphToBigWig'.format(ucsctoolspath), './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_ucsc.bigWig']],
        'bigtools_1thread_multipass_serial':    [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 1', '-p no']],
        'bigtools_2thread_multipass_serial':    [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 2', '-p no']],
        'bigtools_4thread_multipass_serial':    [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 4', '-p no']],
        'bigtools_6thread_multipass_serial':    [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 6', '-p no']],
        'bigtools_8thread_multipass_serial':    [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 8', '-p no']],
        'bigtools_2thread_multipass_parallel':  [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 2', '-p yes']],
        'bigtools_4thread_multipass_parallel':  [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 4', '-p yes']],
        'bigtools_6thread_multipass_parallel':  [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 6', '-p yes']],
        'bigtools_8thread_multipass_parallel':  [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 8', '-p yes']],
        'bigtools_1thread_singlepass_serial':   [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 1', '-p no', '--single-pass']],
        'bigtools_2thread_singlepass_serial':   [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 2', '-p no', '--single-pass']],
        'bigtools_4thread_singlepass_serial':   [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 4', '-p no', '--single-pass']],
        'bigtools_6thread_singlepass_serial':   [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 6', '-p no', '--single-pass']],
        'bigtools_8thread_singlepass_serial':   [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 8', '-p no', '--single-pass']],
        'bigtools_2thread_singlepass_parallel': [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 2', '-p yes', '--single-pass']],
        'bigtools_4thread_singlepass_parallel': [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 4', '-p yes', '--single-pass']],
        'bigtools_6thread_singlepass_parallel': [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 6', '-p yes', '--single-pass']],
        'bigtools_8thread_singlepass_parallel': [['{}/bedgraphtobigwig'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigWig', '-t 8', '-p yes', '--single-pass']],
    }
    compare(comp, 'bedgraphtobigwig_medium', benchmarks)

def bedtobigbed(comp):
    global ucsctoolspath
    global bigtoolspath
    if not os.path.exists('./workdir/ENCFF076CIO.sorted.bed'):
        process = subprocess.check_call('cat ./workdir/ENCFF076CIO.bed | cut -f1-3 | sort -k1,1 -k2,2n > ./workdir/ENCFF076CIO.sorted.bed', shell=True)
    benchmarks = {
        'ucsc':             [['{}/bedToBigBed'.format(ucsctoolspath), './workdir/ENCFF076CIO.sorted.bed', './workdir/hg38.chrom.sizes', './workdir/test_out_ucsc.bigBed']],
        'bigtools_1thread': [['{}/bedtobigbed'.format(bigtoolspath),  './workdir/ENCFF076CIO.sorted.bed', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigBed', '-t 1']],
        'bigtools_2thread': [['{}/bedtobigbed'.format(bigtoolspath),  './workdir/ENCFF076CIO.sorted.bed', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigBed', '-t 2']],
        'bigtools_4thread': [['{}/bedtobigbed'.format(bigtoolspath),  './workdir/ENCFF076CIO.sorted.bed', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigBed', '-t 4']],
        'bigtools_6thread': [['{}/bedtobigbed'.format(bigtoolspath),  './workdir/ENCFF076CIO.sorted.bed', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigBed', '-t 6']],
        'bigtools_8thread': [['{}/bedtobigbed'.format(bigtoolspath),  './workdir/ENCFF076CIO.sorted.bed', './workdir/hg38.chrom.sizes', './workdir/test_out_bigtools.bigBed', '-t 8']],
    }
    compare(comp, 'bedtobigbed', benchmarks)

def bigwigtobedgraph(comp):
    global ucsctoolspath
    global bigtoolspath
    benchmarks = {
        'ucsc':             [['{}/bigWigToBedGraph'.format(ucsctoolspath), './workdir/ENCFF841DHZ.bigWig', '/dev/null']],
        'bigtools_1thread': [['{}/bigwigtobedgraph'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bigWig', '/dev/null', '-t 1']],
        'bigtools_2thread': [['{}/bigwigtobedgraph'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bigWig', '/dev/null', '-t 2']],
        'bigtools_4thread': [['{}/bigwigtobedgraph'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bigWig', '/dev/null', '-t 4']],
        'bigtools_6thread': [['{}/bigwigtobedgraph'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bigWig', '/dev/null', '-t 6']],
        'bigtools_8thread': [['{}/bigwigtobedgraph'.format(bigtoolspath),  './workdir/ENCFF841DHZ.bigWig', '/dev/null', '-t 8']],
    }
    compare(comp, 'bigwigtobedgraph', benchmarks)


def main():
    global ucsctoolspath
    global bigtoolspath
    ucsctoolspath = sys.argv[1] if len(sys.argv) > 1 else "/bin"
    bigtoolspath = sys.argv[2] if len(sys.argv) > 2 else "/usr/local/bin"
    bench = sys.argv[3] if len(sys.argv) > 3 else None
    download_test_data()
    download_kent_tools()
    print("Running benchmarks.")
    with open("./workdir/output/comparison.txt", "w") as comp:
        if bench is None or bench == 'bigwigaverageoverbed':
            bigwigaverageoverbed(comp)
        if bench is None or bench == 'bigwigaverageoverbed_long':
            bigwigaverageoverbed_long(comp)
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
