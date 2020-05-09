#!/usr/bin/env python

import yaml
import subprocess
import sys

with open(".travis.yml", 'r') as stream:
    try:
        d = yaml.safe_load(stream)
        for i, cmd in enumerate(d['before_install']+d['script']):
            print("Running test {}: {}".format(i, cmd))
            p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = p.communicate()
            if p.returncode:
                print(stdout.decode("utf-8"))
                sys.stderr.write(stderr.decode("utf-8"))
                sys.stderr.write("Fail test {}: {}".format(i, cmd))
                break
    except yaml.YAMLError as exc:
        print(exc)
