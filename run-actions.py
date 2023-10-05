#!/usr/bin/env python
import os

import yaml
import subprocess
import sys
import argparse
import re
import cpplint


def run_cmd(cmd, i=1, n=1, j=0, pipe=subprocess.PIPE, name=''):
    retry = (" (retry %d)" % j) if j > 0 else ""
    print("Run test ({:2d}/{}) {}{}: {}".format(i, n, name, retry, cmd))
    p = subprocess.Popen(cmd.split(), stdout=pipe, stderr=pipe)
    stdout, stderr = p.communicate()
    if p.returncode:
        err = ""
        if stdout:
            err += stdout.decode("utf-8")
        if stderr:
            err += stderr.decode("utf-8")
        err += "Fail test {} on try {}: {}".format(i, j, cmd)
        return False, err
    return True, ""


def run_travis(build=True, tag=None, script=None, pipe=subprocess.PIPE):
    with open(".travis.yml", 'r') as stream:
        try:
            d = yaml.safe_load(stream)
            cmds = d['before_install'] + d['script'] if build else d['script']
            for i, cmd in enumerate(cmds):
                if tag:
                    s = cmd.split(tag)
                    cmd = script + s[1]
                passed, err = run_cmd(cmd, i, len(cmds), pipe)
                if not passed:
                    break
        except yaml.YAMLError as exc:
            print(exc)


def run_actions(name='.', tag='ns', docker=True, pipe=subprocess.PIPE, start=0, end=None):
    if docker:
        passed, err = run_cmd('docker build -t {} .'.format(tag))
        if not passed:
            return
        run = 'docker run --rm {}'.format(tag)
    else:
        run = './docker-entrypoint.sh'
    with open(".github/workflows/docker-scheduler-tests.yml", 'r') as stream:
        try:
            d = yaml.safe_load(stream)
            i = -1
            n = len([s for j in d['jobs'].values() for s in j['steps']
                    if 'uses' in s and s['uses'] == './'])
            for j in d['jobs'].values():
                for s in j['steps']:
                    if 'with' not in s or 'ns-args' not in s['with']:
                        continue
                    i += 1
                    if i < start:
                        continue
                    if end is not None and i > end:
                        break
                    if re.search(name, s['name'], re.IGNORECASE) or \
                       re.search(name, s['with']['ns-args'], re.IGNORECASE):
                        cmd = '-i {} {}'.format(s['with']['instance'], s['with']['ns-args'])
                        # cmd = "-i " + cmd.split('-i ')[-1]
                        k = 0
                        retries = s['with'].get('retries', 0)
                        passed = False
                        while not passed and k <= retries:
                            passed, err = run_cmd("{} {}".format(run, cmd), i, n, k, pipe, name=s['name'])
                            k += 1
                        if not passed:
                            print(err)
                            return
        except yaml.YAMLError as exc:
            print(exc)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--docker', dest='docker', action='store_true',
                        help='Set if using docker. Default: False')
    parser.add_argument('-n', '--name', dest='name', default='.',
                        help='Name of the instances to solve. '
                             'Default: match anything, thus "."')
    parser.add_argument('-t', '--tests', dest='tests', default='',
                        help='Range of the tests to perform. '
                             '3 means that all the tests after 3 (included) will be performed. '
                             '-6 means that all the tests before 6 (included) will be performed. '
                             '4-8 means that all the tests between 4 and 8 (included) will be performed'
                             'Default: all "".')
    parser.add_argument('-l', '--lint', dest='lint', action='store_true',
                        help='Use this flag to run cpp lint on the sources. Default: False')
    parser.add_argument('-jl', '--just_lint', dest='just_lint', action='store_true',
                        help='Use this flag to run cpp lint on the sources and deactivate the tests. Default: False')
    args = parser.parse_args()

    if args.just_lint:
        args.lint = True

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

    if args.lint:
        print("################### CPP LINT #########################")
        os.system("python3 -m cpplint --recursive ./src/")
        print("##################### DONE ###########################")

    if not args.just_lint:
        print("################### RUN TESTS ########################")
        if args.docker:
            run_actions(args.name, start=i, end=j)
        else:
            run_actions(args.name, docker=False, start=i, end=j)
        print("##################### DONE ###########################")
