#!/usr/bin/env python

import yaml
import subprocess
import sys
import argparse
import io
import time
import re
import glob
import shutil
import os
import math
from datetime import datetime
from copy import deepcopy

stats_field = '# Final statistics ='
status_field = '# Solution status ='
ub_field = '# Upper bounds ='
lb_field = '# Lower bounds ='

bcp_indexes = {
    'ub': 6,
    'root lb': 7,
    'lb': 8
}

stats_indexes = {
    'ub': 1,
    'root lb': 2,
    'lb': 3
}

output_keys = ['root-lb', 'lb', 'ub', 'status', 'stats', 'lbs', 'ubs', 'time']

markers_plt = {'H': '.', 'D': 'v', 'B': 'x'}
markers_name = {'H': 'Heuristic', 'D': 'Dive', 'B': 'Branch'}


def is_INRC2_instance(name):
    inst = name.split('_')
    return len(inst) == 3 and len(inst[1]) == 1


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


def find_fields(stream, field, func=None):
    b = io.StringIO(stream)

    # search for last occurrence of stats_field
    lines = []
    i = stream.find(field)
    while i >= 0:
        # search root LB, best LB, UB
        b.seek(i)
        line = b.readline()
        line2 = line.split('=', 1)[-1]
        new_line = ' '.join(s for s in line2.split() if s)
        # if found same line, break
        if len(lines) > 0 and new_line == lines[-1]:
            break
        if func is not None:
            func(line2)
        lines.append(new_line)
        i = stream.find(field, i + 1)

    return lines


def find_stats_bounds(stream):
    b = io.StringIO(stream)

    # search for last occurrence of stats_field
    bounds = []

    def extract_bounds(line):
        r = find_bounds(line, stats_indexes)
        bounds.append(r)
    stats = find_fields(stream, stats_field, extract_bounds)
    if len(bounds) == 0:
        raise ValueError("No bounds found")
    # if len(bounds) >= 2:
    #     sum_bounds = [sum(bds[i] for bds in bounds) for i in range(3)]
    #     sum_bounds.append('')
    #     bounds.append(sum_bounds)
    #     stats.append('N/A')
    r_lb, lb, ub, l = bounds[-1]

    # bounds history
    lbs = find_fields(stream, lb_field)
    ubs = find_fields(stream, ub_field)

    # find status
    i = stream.rfind(status_field)
    if i >= 0:
        b.seek(i)
        status = b.readline().split()[-1]
    else:
        status = 'N/A'

    # do not return last stat if solving by components
    return r_lb, lb, ub, status, stats[:-1] if len(stats) > 1 else stats, lbs, ubs


def run_cmd(cmd, i=0, n=1, pipe=True, logPath=None, r=None):
    st = datetime.now().strftime('%H:%M:%S')
    print("[{}] Running test {:2d}/{}{}: {}".format(
        st, i, n, "" if r is None else " (%d)" % r, cmd))
    start = time.time()
    out = subprocess.PIPE if pipe else sys.stdout
    if logPath is not None:
        out = open("%s/log.txt" % logPath, 'w')
    p = subprocess.Popen(cmd.split(), stdout=out, stderr=out)
    stdout, stderr = p.communicate()
    end = int(time.time() - start)
    print("run time: {} s.".format(end))
    if logPath is not None:
        with open("%s/log.txt" % logPath, 'rb') as f:
            stdout = f.read()
    try:
        if p.returncode > 0:
            raise "ended with a positive return code: {}".format(p.returncode)
        return (*find_stats_bounds(stdout.decode("utf-8")), end) if stdout else {}
    except:
        sys.stderr.write("Fail test {}: {}".format(i, cmd))
        raise


def get_args(d_yaml, ns_args='', param=None):
    for k, v in d_yaml.items():
        if k in ['instances', 'exe', 'steps', 'replications', 'addparams', 'param']:
            continue
        elif k == 'args':
            if v:
                ns_args += ' {}'.format(v)
        else:
            ns_args += ' --{} {}'.format(k, v)
    if param:
        ns_args += ' --{} {}'.format('param', param)
    return ns_args


def add_params(d_yaml, sol_path, param=None):
    param = d_yaml.get('param', param)
    if param is None or 'addparams' not in d_yaml:
        return param

    p_name = param.rsplit('/')[-1]
    n_param = '%s/%s' % (sol_path, p_name)
    shutil.copyfile(param, n_param)
    b = '\n'
    for k, v in d_yaml['addparams'].items():
        b += '%s=%s\n' % (k, v)
    with open(n_param, 'a') as f:
        f.write(b)
    return n_param


def skip_instance(t, status, name):
    if 'status' in t and not re.search(status, t['status'], re.IGNORECASE):
        return True
    if not re.search(name, t['name'], re.IGNORECASE):
        return True
    return False


def run_benchmark(benchmark_file, exe, pipe=True,
                  status='.', name='.', break_on_error=False,
                  start=0, end=None, prefix=''):
    with open(benchmark_file, 'r') as stream:
        try:
            o_file_modified = False
            d = yaml.safe_load(stream)
            if not exe:
                exe = d['benchmark']['exe']
            n_replications = d['benchmark'].get('replications', 1)
            p_time = '/%s%d' % (prefix, int(time.time()))
            b_name = benchmark_file.rsplit('.', 1)[0]
            b_name = b_name.rsplit('/')[-1]
            o_file = 'results/%s%s.yml' % (b_name, p_time) if pipe else None
            if o_file:
                os.makedirs(o_file.rsplit('/', 1)[0], exist_ok=True)
            sol_path = 'outfiles/%s%s' % (b_name, p_time)
            os.makedirs(sol_path, exist_ok=True)
            i = 0
            d_args = ''
            param = None
            instances = None
            if 'steps' in d['benchmark']:
                d_args = get_args(d['benchmark'])
                param = add_params(d['benchmark'], sol_path, param)
                dsteps = d['benchmark']['steps']
                instances = d['benchmark'].get('instances', None)
            else:
                dsteps = [d['benchmark']]
            n_tests = sum(len(s.get('instances', instances)) for s in dsteps)
            for j, s in enumerate(dsteps):
                sol_step_path = ('%s/step%d' % (sol_path, j)) if 'steps' in d['benchmark'] else sol_path
                os.makedirs(sol_step_path, exist_ok=True)
                s_param = add_params(s, sol_step_path, param)
                ns_args = get_args(s, d_args, s_param)
                if 'instances' not in s:
                    s['instances'] = deepcopy(instances)
                for t in s['instances']:
                    if skip_instance(t, status, name):
                        continue
                    if i < start or (end is not None and i > end):
                        i += 1
                        continue
                    # check if inrc2 format
                    i_name = t['name']
                    if is_INRC2_instance(t['name']):
                        cmd = '--instance {a[0]} --his {a[1]} ' \
                              '--weeks {a[2]} '.format(a=t['name'].split('_'))
                    else:
                        cmd = '--instance {} '.format(t['name'])
                        i_name = t['name'].rsplit('.', 1)[0]
                    # add common args
                    cmd += ns_args
                    # add local dir if any
                    if 'dir' in t:
                        cmd += ' --dir ' + t['dir']
                    l_exe = exe
                    if 'exe' in t:
                        l_exe = t['exe']
                    for p in range(n_replications):
                        # add solution path
                        n_sol_path = '%s/%d' % (sol_step_path, p) if n_replications > 1 else sol_step_path
                        sol_inst_path = '%s/%s' % (n_sol_path, i_name)
                        os.makedirs(sol_inst_path, exist_ok=True)
                        cmd += ' --sol {0}/'.format(sol_inst_path)
                        try:
                            bds = run_cmd("{} {}".format(l_exe, cmd), i+1, n_tests,
                                          pipe, sol_inst_path, p if n_replications > 1 else None)
                            if 'ub' in bds and 'ub' in t and bds['ub'] != t['ub']:
                                print("The UB is different than the previous one.")
                            for k, v in enumerate(bds):
                                if p == 0:
                                    t[output_keys[k]] = v
                                elif isinstance(t[output_keys[k]], list):
                                    t[output_keys[k]] += v
                                elif isinstance(t[output_keys[k]], int) or \
                                        isinstance(t[output_keys[k]], float):
                                    t[output_keys[k]] = min(t[output_keys[k]], v)
                                else:
                                    t[output_keys[k]] = v
                        except KeyboardInterrupt:
                            return
                        except Exception as e:
                            if break_on_error:
                                raise
                            print(e)
                            t['status'] = 'ERROR'
                        if o_file:
                            with open(o_file, 'w+') as wstream:
                                yaml.dump(d, wstream, sort_keys=False)
                            o_file_modified = True
                    i += 1
            return o_file if o_file_modified else None
        except yaml.YAMLError as exc:
            print(exc)
            return None


def create_graph(lbs_line, ubs_line, x_max=None, file=None, gap=1.25, granularity=50):
    import matplotlib.pyplot as plt

    lbs = [[], []]  # time, value
    for p in lbs_line.split():
        lb, t = p.split(',')
        lbs[0].append(float(t))
        lbs[1].append(float(lb))

    ubs = [[], []]  # time, value
    markers = {k:[] for k in markers_plt}
    for i, p in enumerate(ubs_line.split()):
        ub, t, m = p.split(',')
        ubs[0].append(float(t))
        ubs[1].append(int(ub))
        markers[m].append(i)

    # # ends the lines
    # if ubs[0][-1] < lbs[0][-1]:
    #     # add a point on the last time
    #     ubs[0].append(lbs[0][-1])
    #     ubs[1].append(ubs[1][-1])
    # else:
    #     # add a point on the last time
    #     lbs[0].append(ubs[0][-1])
    #     lbs[1].append(lbs[1][-1])

    if x_max:
        lbs[0].append(x_max)
        ubs[0].append(x_max)
        lbs[1].append(lbs[1][-1])
        ubs[1].append(ubs[1][-1])


    plt.figure()
    plt.plot(lbs[0], lbs[1], ':g.', label="lb")
    # plt.plot(ubs[0], ubs[1], '-bx', label="ub")
    for n, m in markers_plt.items():
        plt.plot(ubs[0], ubs[1], '-b'+m, label="ub - "+markers_name[n], markevery=markers[n])

    ax = plt.gca()
    ax.set_xlim([0, x_max])
    y_min = (math.floor(lbs[1][0]) // granularity) * granularity
    y_max = (math.ceil(y_min * gap) // granularity) * granularity
    ax.set_ylim([y_min, y_max])

    plt.legend()
    if file:
        plt.savefig(file)
    else:
        plt.show()


def validate(name, b_name, dir, ub, env):
    sol_path = 'outfiles/%s' % b_name
    i_name = name.split('.')[0]
    if is_INRC2_instance(name):
        print("Not implemented.")
    else:
        inst = "%sinstances_xml/%s.xml" % (dir, i_name)
        sol = "%s/%s/%s.xml" % (sol_path, i_name, i_name)
        cmd = 'java -jar validatorINRC.jar %s %s' % (inst, sol)
        p = subprocess.Popen(cmd, shell=True, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if stderr:
            print(stderr.decode("utf-8"))
            return False
        try:
            buf = io.StringIO(stdout.decode("utf-8"))
            lines = buf.readlines()
            if len(lines) < 2:
                return False
            m = re.search("[0-9]+", lines[-2])
            if m.group():
                if ub == int(m.group()):
                    return True
                print("Found UB="+m.group())
            return False
        except Exception as e:
            print(e)
            return False


def run_validator(benchmark_file, env, status='.', name='.'):
    with open(benchmark_file, 'r') as stream:
        try:
            d = yaml.safe_load(stream)
            dsteps = d['benchmark']['steps'] if 'steps' in d['benchmark'] \
                else [d['benchmark']]
            instances = d['benchmark'].get('instances', None)
            b_path = benchmark_file.rsplit('.', 1)[0]
            b_name = b_path.rsplit('/')[-1]
            base_dir = d['benchmark'].get('dir', '')
            for i, s in enumerate(dsteps):
                b_step = ('%s/step%d' % (b_name, i)) if 'steps' in d['benchmark'] else b_name
                dir = s.get('dir', base_dir)
                if 'instances' not in s:
                    s['instances'] = deepcopy(instances)
                for t in s['instances']:
                    if skip_instance(t, status, name):
                        continue
                    if 'ub' in t:
                        if validate(t['name'], b_step, dir, t['ub'], env):
                            print('%s: Valid UB (%d)' % (t['name'], t['ub']))
                        else:
                            print('%s: Invalid UB (%d)' % (t['name'], t['ub']))
                    else:
                        print('%s: No UB available' % t['name'])
        except Exception as e:
            print(e)
            return False


def extract_results(results_file, plot=False, status='.', name='.'):
    with open(results_file, 'r') as stream:
        try:
            d = yaml.safe_load(stream)
            print('Results for '+results_file)
            print("Instance, Status, Time, rootLB, LB, UB, Time")
            if 'steps' in d['benchmark']:
                dsteps = d['benchmark']['steps']
            else:
                dsteps = [d['benchmark']]
            for i, s in enumerate(dsteps):
                step_name = "%d_" % i if len(dsteps) > 1 else ""
                if len(dsteps) > 1:
                    print("Step %d:" % i)
                if 'instances' not in s:
                    continue
                for r in s['instances']:
                    if skip_instance(r, status, name):
                        continue
                    try:
                        lb = r['lb']
                        res = "%s, %s, %s, %s, %s, %s" % \
                              (r['name'], r['status'], r['time'],
                               r['root-lb'], lb, r['ub'])
                        if 'ubs' in r and len(r['ubs']) > 1:
                            tub = 0
                            for l in r['ubs'][-2:]:
                                t = float((l.split()[-1]).split(',')[1])
                                if t > tub:
                                    tub = t
                            res += ", %.0f" % tub
                        print(res)
                        if plot and 'lbs' in r and 'ubs' in r:
                            f = results_file
                            s = results_file.split('/')
                            if s:
                                f = s[-1]
                            f = "images/%s_%s%s" % (r['name'], step_name, f.rsplit('.', 1)[0])
                            for j, l in enumerate(r['lbs']):
                                fn = "%s_%d.png" % (f, j)
                                create_graph(l, r['ubs'][j], r['time'], fn)
                    except KeyError:
                        print("%s,,,,,," % r['name'])
        except yaml.YAMLError as exc:
            print(exc)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--benchmark', dest='benchmark',
                        default='./benchmark/w4.yml',
                        help='Benchmark to run. Default: benchmark/w4.yml')
    parser.add_argument('-d', '--docker', dest='docker', action='store_true',
                        help='Set if using docker. Default: False')
    parser.add_argument('-s', '--status', dest='status', default='.',
                        help='Status of the instances to solve. '
                             'Default: match anything, thus "."')
    parser.add_argument('-np', '--no-pipe', dest='pipe', action='store_false',
                        help='Do not redirect.')
    parser.add_argument('-n', '--name', dest='name', default='.',
                        help='Name of the instances to solve. '
                             'Default: match anything, thus "."')
    parser.add_argument('-r', '--results', dest='results', action='store_true',
                        help='Set to true if only processing the results. Default: False')
    parser.add_argument('-p', '--plot', dest='plot', action='store_true',
                        help='Set to true if only plotting the results. Default: False')
    parser.add_argument('-v', '--validator', dest='validator', action='store_true',
                        help='Run the validator on the instance present if the file. Default: False')
    parser.add_argument('--break', dest='break_on_error', action='store_true',
                        help='Set to true if breaking the program when having an error. Default: False')
    parser.add_argument('-t', '--tests', dest='tests', default='',
                        help='Range of the tests to run. '
                             '3 means that all the tests after 3 (included) will be performed. '
                             '-6 means that all the tests before 6 (included) will be performed. '
                             '4-8 means that all the tests between 4 and 8 (included) will be performed'
                             'Default: all "".')
    parser.add_argument('--prefix', dest='prefix', default='',
                        help='Timestamp prefix used for the results '
                             'file/folder. Default: \"\".')
    args = parser.parse_args()

    i = 0
    j = None
    if args.tests and args.tests != '-':
        if '-' not in args.tests:
            i = int(args.tests)
        elif args.tests.startswith('-'):
            j = int(args.tests[1:])
        else:
            s = args.tests.split('-')
            i = int(s[0])
            if s[1]:
                j = int(s[1])
    if j is not None and j < i:
        raise ValueError("The end indice of the test needs to be "
                         "greater than its start.")

    if args.docker:
        run_cmd('docker build -t {} .'.format('ns'), pipe=args.pipe)
        run = 'docker run --rm {}'.format('ns')
    else:
        run = ''

    benchmark = []
    for p in args.benchmark.split(','):
        benchmark += [f for f in glob.glob(p) if f.endswith('.yml') or f.endswith('.yaml')]
    results = benchmark
    print("Run these benchmarks: %s" %benchmark)
    if args.validator:
        env = dict(os.environ)
        for b in benchmark:
            run_validator(b, env, status=args.status, name=args.name)
    else:
        if not args.results and not args.plot:
            results = []
            for b in benchmark:
                r = run_benchmark(b, run, status=args.status, name=args.name,
                                  break_on_error=args.break_on_error,
                                  pipe=args.pipe, start=i, end=j,
                                  prefix=args.prefix)
                if not r:
                    break
                results.append(r)

        for r in results:
            if r:
                extract_results(r, args.plot, status=args.status, name=args.name)
