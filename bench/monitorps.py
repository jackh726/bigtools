# Based on psrecord https://github.com/astrofrog/psrecord

# Copyright (c) 2013, Thomas P. Robitaille
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import subprocess
import time
import argparse
import psutil


def main():
    parser = argparse.ArgumentParser(description='Record CPU and memory usage for a process')
    parser.add_argument('command', type=str, help='the process command')
    parser.add_argument('--log', type=str, help='output the statistics to a file')
    parser.add_argument('--interval', type=float, help='how long to wait between each sample (in seconds). By default the process is sampled as often as possible.')
    args = parser.parse_args()
    start_and_monitor(args.command, args.log, args.interval)

def start_and_monitor(command, logfile, interval):
    print("Starting up command '{0}' and attaching to process".format(command))
    sprocess = subprocess.Popen(command, shell=True)
    pid = sprocess.pid
    monitor(pid, logfile, interval)
    sprocess.kill()

def monitor(pid, logfile, interval):
    pr = psutil.Process(pid)
    with open(logfile, "w") as f:
        f.write("# {0:12s} {1:12s} {2:12s} {3:12s}\n".format(
            'Elapsed time'.center(12),
            'CPU (%)'.center(12),
            'Real (MB)'.center(12),
            'Virtual (MB)'.center(12))
        )

        # Record start time
        start_time = time.time()

        children = {}

        # Start main event loop
        while True:
            # Find current time
            current_time = time.time()

            current_cpu = 0
            current_mem_real = 0
            current_mem_virt = 0

            # Process may finish at any point
            try:
                with pr.oneshot():
                    pr_status = pr.status()
                    # Check if process status indicates we should exit
                    if pr_status in [psutil.STATUS_ZOMBIE, psutil.STATUS_DEAD]:
                        break

                    current_cpu += pr.cpu_percent()
                    mem = pr.memory_info()
                    current_mem_real += mem.rss / 1024. ** 2
                    current_mem_virt += mem.vms / 1024. ** 2
                
                for child in pr.children(recursive=True):
                    if child.pid in children:
                        child_p = children[child.pid]
                    else:
                        children[child.pid] = child
                        child_p = child

                    with child_p.oneshot():
                        current_cpu += child_p.cpu_percent()
                        mem = child_p.memory_info()
                        current_mem_real += mem.rss / 1024. ** 2
                        current_mem_virt += mem.vms / 1024. ** 2

                f.write("{0:12.3f} {1:12.3f} {2:12.3f} {3:12.3f}\n".format(
                    current_time - start_time,
                    current_cpu,
                    current_mem_real,
                    current_mem_virt))
                f.flush()
            except:
                break

            if interval is not None:
                time.sleep(interval)

if __name__ == "__main__":
    main()