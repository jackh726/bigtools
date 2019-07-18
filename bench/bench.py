import os
import re
import resource
import struct
import subprocess
import sys
import urllib.request

timeregex = re.compile(r'(\d.*)user (\d.*)system (\d.*)elapsed')

def download_test_data():
    print("Downloading test data")
    if not os.path.exists("./workdir"):
        os.makedirs("./workdir")
    if not os.path.exists('./workdir/ENCFF937MNZ.bigWig'):
        urllib.request.urlretrieve('https://www.encodeproject.org/files/ENCFF937MNZ/@@download/ENCFF937MNZ.bigWig', './workdir/ENCFF937MNZ.bigWig')
    if not os.path.exists('./workdir/ENCFF447DHW.bigWig'):
        urllib.request.urlretrieve('https://www.encodeproject.org/files/ENCFF447DHW/@@download/ENCFF447DHW.bigWig', './workdir/ENCFF447DHW.bigWig')
    if not os.path.exists('./workdir/ENCFF841DHZ.bigWig'):
        urllib.request.urlretrieve('https://www.encodeproject.org/files/ENCFF841DHZ/@@download/ENCFF841DHZ.bigWig', './workdir/ENCFF841DHZ.bigWig')
    if not os.path.exists('./workdir/ENCFF518WII.bigWig'):
        urllib.request.urlretrieve('https://www.encodeproject.org/files/ENCFF518WII/@@download/ENCFF518WII.bigWig', './workdir/ENCFF841DHZ.bigWig')
    if not os.path.exists('./workdir/ENCFF646AZP.bed'):
        urllib.request.urlretrieve('https://www.encodeproject.org/files/ENCFF646AZP/@@download/ENCFF646AZP.bed.gz', './workdir/ENCFF646AZP.bed.gz')
        subprocess.check_call(['gunzip', './workdir/ENCFF646AZP.bed.gz'])
    if not os.path.exists('./workdir/hg38.chrom.sizes'):
        urllib.request.urlretrieve('https://raw.githubusercontent.com/ENCODE-DCC/encValData/562ab5bf03deff9bb5340991fd5c844162b82914/GRCh38/GRCh38_EBV.chrom.sizes', './workdir/hg38.chrom.sizes')
    print("Done downlaoding test data")



# from https://stackoverflow.com/a/13889698
def start_running(command):
    time_read_pipe, time_write_pipe = os.pipe()
    want_read_pipe, want_write_pipe = os.pipe()
    runner_pid = os.fork()
    if runner_pid != 0:
        os.close(time_write_pipe)
        os.close(want_read_pipe)
        def finish_running():
            os.write(want_write_pipe, 'x')
            os.close(want_write_pipe)
            time = os.read(time_read_pipe, struct.calcsize('f'))
            os.close(time_read_pipe)
            time = struct.unpack('f', time)[0]
            return time
        return finish_running
    os.close(time_read_pipe)
    os.close(want_write_pipe)
    sub_pid = os.fork()
    if sub_pid == 0:
        os.close(time_write_pipe)
        os.close(want_read_pipe)
        os.execvp(command[0], command)
    os.wait()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    os.read(want_read_pipe, 1)
    os.write(time_write_pipe, struct.pack('f', usage.ru_utime))
    sys.exit(0)

def time(exeargs_all):
    total_seconds = 0
    for exeargs in exeargs_all:
        exeargs.insert(0, 'time')
        exeargs = " ".join(exeargs)
        print(exeargs)
        process = subprocess.Popen(exeargs, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
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

def compare(bench, ucsc, bigwig2):
    print('Benchmarking {}'.format(bench))
    ucsctime = time(ucsc)
    bigwig2time = time(bigwig2)
    print("ucsc: {}".format(round(ucsctime,3)))
    print("bigwig2: {}".format(round(bigwig2time,3)))

def bigwigaverageoverbed():
    # For ucsc, we have to convert narrowPeak to bed first, including adding a unique name
    ucsc = [
        ['cat ./workdir/ENCFF646AZP.bed | cut -f1-3 | awk -v OFS=\'\\t\' \'{print $1,$2,$3, NR}\' > ./workdir/ENCFF646AZP_cut.bed'],
        ['{}/bigWigAverageOverBed'.format(ucsctoolspath), './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF646AZP_cut.bed', './workdir/test_out_ucsc.bed']
        ]
    bigwig2 = [['{}/bigwigaverageoverbed'.format(bigwig2path), '-s' , './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF646AZP.bed', './workdir/test_out_bigwig2.bed']]
    compare('bigwigaverageoverbed', ucsc, bigwig2)

def bigwigmerge():
    ucsc = [['{}/bigWigMerge'.format(ucsctoolspath), './workdir/ENCFF937MNZ.bigWig', './workdir/ENCFF447DHW.bigWig', './workdir/test_out_ucsc.bedGraph']]
    bigwig2 = [['{}/bigwigmerge'.format(bigwig2path), './workdir/test_out_bigwig2.bigWig', '-b ./workdir/ENCFF937MNZ.bigWig', '-b ./workdir/ENCFF447DHW.bigWig']]
    compare('bigwigmerge', ucsc, bigwig2)

def bedgraphtobigwig():
    # Need to generate bedGraph first, just use ucsc since it's what we're comparing against
    if not os.path.exists('./workdir/ENCFF518WII.bedGraph'):
        process = subprocess.check_call('{}/bigWigToBedGraph ./workdir/ENCFF518WII.bigWig ./workdir/ENCFF518WII.bedGraph'.format(ucsctoolspath), shell=True)
    if not os.path.exists('./workdir/ENCFF841DHZ.bedGraph'):
        process = subprocess.check_call('{}/bigWigToBedGraph ./workdir/ENCFF841DHZ.bigWig ./workdir/ENCFF841DHZ.bedGraph'.format(ucsctoolspath), shell=True)
    ucsc = [['{}/bedGraphToBigWig'.format(ucsctoolspath), './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_ucsc.bigWig']]
    bigwig2 = [['{}/bedgraphtobigwig'.format(bigwig2path), './workdir/ENCFF518WII.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigwig2.bigWig']]
    compare('bedgraphtobigwig', ucsc, bigwig2)
    ucsc = [['{}/bedGraphToBigWig'.format(ucsctoolspath), './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_ucsc.bigWig']]
    bigwig2 = [['{}/bedgraphtobigwig'.format(bigwig2path), './workdir/ENCFF841DHZ.bedGraph', './workdir/hg38.chrom.sizes', './workdir/test_out_bigwig2.bigWig']]
    compare('bedgraphtobigwig', ucsc, bigwig2)

def bigwigtobedgraph():
    ucsc = [['{}/bigWigToBedGraph'.format(ucsctoolspath), './workdir/ENCFF841DHZ.bigWig', './workdir/test_out_ucsc.bedGraph']]
    bigwig2 = [['{}/bigwigtobedgraph'.format(bigwig2path), './workdir/ENCFF841DHZ.bigWig', './workdir/test_out_bigwig2.bedGraph']]
    compare('bigwigtobedgraph', ucsc, bigwig2)


def main():
    global ucsctoolspath
    global bigwig2path
    ucsctoolspath = sys.argv[1]
    bigwig2path = sys.argv[2]
    download_test_data()
    bigwigaverageoverbed()
    bigwigmerge()
    bedgraphtobigwig()
    bigwigtobedgraph()

if __name__ == '__main__':
    main()

# To run (example): python3 bench/bench.py ~/temp /mnt/c/Users/jackh/Documents/Git/rust/bigwig2/target/release