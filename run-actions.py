#!/usr/bin/env python
import yaml
import subprocess
import sys
import argparse
import re


def run_cmd(cmd, i=1, n=1, j=0, pipe=subprocess.PIPE):
    retry = (", retry %d" % j) if j > 0 else ""
    print("Run test {:2d}/{}{}: {}".format(i, n, retry, cmd))
    p = subprocess.Popen(cmd.split(), stdout=pipe, stderr=pipe)
    stdout, stderr = p.communicate()
    if p.returncode:
        err = ""
        if stdout:
            err += stdout.decode("utf-8")
        if stderr:
            err += stderr.decode("utf-8")
        err = "Fail test {} on try {}: {}".format(i, j, cmd)
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


def run_actions(name='.', tag='ns', docker=True, pipe=subprocess.PIPE):
    if docker:
        passed, err = run_cmd('docker build -t {} .'.format(tag))
        if not passed:
            return
        run = 'docker run --rm {}'.format(tag)
    else:
        run = './docker-entrypoint.sh'
    with open(".github/workflows/docker-tests.yml", 'r') as stream:
        try:
            d = yaml.safe_load(stream)
            i = 0
            for j in d['jobs'].values():
                for s in j['steps']:
                    if 'with' in s and \
                            (re.search(name, s['name'], re.IGNORECASE) or
                             re.search(name, s['with']['ns-args'], re.IGNORECASE)):
                        cmd = '-i {} {}'.format(s['with']['instance'], s['with']['ns-args'])
                        # cmd = "-i " + cmd.split('-i ')[-1]
                        i += 1
                        k = 0
                        retries = s['with'].get('retries', 0)
                        passed = False
                        while not passed and k <= retries:
                            passed, err = run_cmd("{} {}".format(run, cmd), i, len(j['steps']), k, pipe)
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
    args = parser.parse_args()

    if args.docker:
        # run_travis()
        run_actions(args.name)
    else:
        # run_travis(False, 'ns', './docker-entrypoint.sh')
        run_actions(args.name, docker=False)
