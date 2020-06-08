#!/usr/bin/env python

import yaml
import subprocess
import sys
import argparse
import io
import time
import re
from datetime import datetime


bcp_indexes = {
    'ub': 6,
    'root lb': 7,
    'lb': 8
}

stats_indexes = {
    'ub': 3,
    'root lb': 4,
    'lb': 5
}

output_keys = ['root-lb', 'lb', 'ub', 'status', 'stats', 'time']


def find_bounds(line, indexes):
    l_values = line.split()
    print(line[:-1])
    try:
        r_lb = float(l_values[indexes['root lb']])
        lb = float(l_values[indexes['lb']])
        ub = int(float(l_values[indexes['ub']]))
        print('root-lb={}, lb={}, ub={}{}'.format(
            r_lb, lb, ub, '*' if ub - lb < 5 - 1e-5 else ''))
        return r_lb, lb, ub, line
    except ValueError:
        return


def find_bcp_bounds(stream):
    # search for last occurrence of 'BCP'
    i = stream.rfind('BCP:')

    # search root LB, best LB, UB
    b = io.StringIO(stream)
    b.seek(i-1000)
    line = b.readline()  # return 1000 bytes before and go to next line
    # search for last occurrence of UB, ROOT LB, LB
    while b.tell() <= i:
        line = b.readline()
        if line[0:4] == 'BCP:':
            return find_bounds(line, bcp_indexes)


def find_stats_bounds(stream):
    b = io.StringIO(stream)

    # search for last occurrence of 'final statistics='
    stats = []
    i = stream.find('final statistics=')
    while i >= 0:
        # search root LB, best LB, UB
        b.seek(i)
        line = b.readline()
        r_lb, lb, ub, l = find_bounds(line, stats_indexes)
        a = re.compile(r'\s').split(line.split('=', 1)[-1])
        stats.append(' '.join(s for s in a if s))
        i = stream.find('final statistics=', i + 1)

    # find status
    i = stream.rfind('# Solution status =')
    b.seek(i)
    status = b.readline().split()[-1]

    # do not return last stat if solving by components
    return r_lb, lb, ub, status, stats[:-1] if len(stats) > 1 else stats


def run_cmd(cmd, i=0, pipe=subprocess.PIPE):
    st = datetime.now().strftime('%H:%M:%S')
    print("[{}] Running test {}: {}".format(st, i, cmd))
    start = time.time()
    p = subprocess.Popen(cmd.split(), stdout=pipe, stderr=pipe)
    stdout, stderr = p.communicate()
    end = int(time.time() - start)
    print("run time: {} s.".format(end))
    try:
        if p.returncode > 0:
            raise "ended with a positive return code: {}".format(p.returncode)
        if stdout:
            return (*find_stats_bounds(stdout.decode("utf-8")), end)
    except:
        if stdout:
            print(stdout.decode("utf-8"))
        if stderr:
            sys.stderr.write(stderr.decode("utf-8"))
        sys.stderr.write("Fail test {}: {}".format(i, cmd))
        raise


def run_benchmark(benchmark_file, exe, pipe=subprocess.PIPE, status='.'):
    with open(benchmark_file, 'r') as stream:
        try:
            d = yaml.safe_load(stream)
            ns_args = ''
            for k, v in d['benchmark'].items():
                if k == 'instances':
                    instances = v
                elif k == 'exe':
                    if not exe:
                        exe = v
                elif k == 'args':
                    if v:
                        ns_args += ' {}'.format(v)
                else:
                    ns_args += ' --{} {}'.format(k, v)

            o_file = benchmark_file.split('.')[-2]
            o_file += '_{}.yml'.format(int(time.time()))
            for i, t in enumerate(instances):
                if 'status' in t and not re.match(status, t['status'], re.IGNORECASE):
                    continue
                i_args = t['name'].split('_') + [ns_args]
                print(i_args)
                cmd = '--instance {a[0]}  --his {a[1]} ' \
                      '--weeks {a[2]}{a[3]}'.format(a=i_args)
                bds = run_cmd("{} {}".format(exe, cmd), i+1, pipe)
                if bds:
                    for j, v in enumerate(bds):
                        t[output_keys[j]] = v
                    with open(o_file, 'w') as wstream:
                        yaml.dump(d, wstream, sort_keys=False)

        except yaml.YAMLError as exc:
            print(exc)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--benchmark', dest='benchmark',
                        default='benchmark/w4.yml',
                        help='Benchmark to run. Default: benchmark/w4.yml')
    parser.add_argument('-d', '--docker', dest='docker', action='store_true',
                        help='Set if using docker. Default: False')
    parser.add_argument('-s', '--status', dest='status', default='.',
                        help='Status of the instances to solve. '
                             'Default: match anything, thus "."')
    args = parser.parse_args()

    if args.docker:
        run_cmd('docker build -t {} .'.format('ns'))
        run = 'docker run --rm {}'.format('ns')
    else:
        run = ''

    run_benchmark(args.benchmark, run, status=args.status)
