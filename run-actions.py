#!/usr/bin/env python

import yaml
import subprocess
import sys
import argparse
import re


def run_cmd(cmd, i=0, pipe=subprocess.PIPE):
    print("Running test {}: {}".format(i, cmd))
    p = subprocess.Popen(cmd.split(), stdout=pipe, stderr=pipe)
    stdout, stderr = p.communicate()
    if p.returncode:
        if stdout:
            print(stdout.decode("utf-8"))
        if stderr:
            sys.stderr.write(stderr.decode("utf-8"))
        sys.stderr.write("Fail test {}: {}".format(i, cmd))
        return False
    return True


def run_travis(build=True, tag=None, script=None, pipe=subprocess.PIPE):
    with open(".travis.yml", 'r') as stream:
        try:
            d = yaml.safe_load(stream)
            cmds = d['before_install'] + d['script'] if build else d['script']
            for i, cmd in enumerate(cmds):
                if tag:
                    s = cmd.split(tag)
                    cmd = script + s[1]
                if not run_cmd(cmd, i, pipe):
                    break
        except yaml.YAMLError as exc:
            print(exc)


def run_actions(name='.', tag='ns', docker=True, pipe=subprocess.PIPE):
    if docker:
        run_cmd('docker build -t {} .'.format(tag))
        run = 'docker run --rm {}'.format(tag)
    else:
        run = './docker-entrypoint.sh'
    with open(".github/workflows/docker-tests.yml", 'r') as stream:
        try:
            d = yaml.safe_load(stream)
            i = 0
            for j in d['jobs'].values():
                for s in j['steps']:
                    if 'with' in s and re.search(name, s['name'], re.IGNORECASE):
                        cmd = '-i {} {}'.format(s['with']['instance'], s['with']['ns-args'])
                        # cmd = "-i " + cmd.split('-i ')[-1]
                        i += 1
                        if not run_cmd("{} {}".format(run, cmd), i, pipe):
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
