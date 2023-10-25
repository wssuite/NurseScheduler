#!/usr/bin/env python
import os
import sys
import argparse
import re
from datetime import timedelta, date

if __name__ == "__main__":
    real_path = os.path.realpath(__file__)
    dir_path = os.path.dirname(real_path)
    sys.path.insert(0, dir_path)

import ui

startDate = "2023-09-04"
defC = "def"
work = "Work"
rest = "Rest"

underCoverage = 30
consShifts = 15
consDaysWork = 30
consDaysOff = 30
pref = 10
complete = 30
totalShifts = 20
totalWeekends = 30
hard = "hard"

weekday_indexes = {
    "Mon": 0,
    "Tue": 1,
    "Wed": 2,
    "Thu": 3,
    "Fri": 4,
    "Sat": 5,
    "Sun": 6,
}

def get_bds(bds):
    bds = bds[1:-1].split(",")
    return int(bds[0]), int(bds[1])


def read_scenario(path, instance, iname):
    s_file = os.path.join(path, instance, "Sc-{}.txt".format(instance))
    nWeeks = int(instance.rsplit("w", 1)[-1])
    endDate = date.fromisoformat(startDate) + timedelta(weeks=nWeeks) - timedelta(days=1)
    iui = ui.ui(iname, startDate, endDate.isoformat())
    with open(s_file) as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if line.startswith("SKILLS"):
                n = int(line.rsplit("=", 1)[-1])
                for j in range(n):
                    i += 1
                    sk = lines[i].strip()
                    iui.add_skill(sk)
                    iui.forbid_skills("N_"+sk, [sk], hard)

            elif line.startswith("SHIFT_TYPES"):
                n = int(line.rsplit("=", 1)[-1])
                for j in range(n):
                    i += 1
                    s = lines[i].strip().split()
                    lb, ub = get_bds(s[-1])
                    iui.add_shift(s[0], "0:0", "0:0")
                    iui.add_shift_type("T_"+s[0], [s[0]])
                    iui.add_min_max_cons_shift(
                        defC, s[0], lb, consShifts,
                        ub, consShifts)

            elif line.startswith("FORBIDDEN_SHIFT_TYPES_SUCCESSIONS"):
                n = iui.n_shifts
                for j in range(n):
                    i += 1
                    succ = lines[i].strip().split()
                    if len(succ) > 2:
                        pat = ui.Pattern()
                        pat.add_shifts([succ[0]])
                        pat.add_shifts(succ[2:])
                        iui.add_unwanted_pattern(defC, pat, hard)

            elif line.startswith("CONTRACTS"):
                n = int(line.rsplit("=", 1)[-1])
                for j in range(n):
                    i += 1
                    s = lines[i].strip().split()
                    n4 = nWeeks / 4

                    lb, ub = get_bds(s[1])
                    iui.add_min_max_assignment_in_4_weeks(
                        s[0], work, round(lb / n4), totalShifts,
                        round(ub / n4), totalShifts)

                    lb, ub = get_bds(s[2])
                    iui.add_min_max_cons_shift(
                        s[0], work, lb, consDaysWork, ub, consDaysWork)

                    lb, ub = get_bds(s[3])
                    iui.add_min_max_cons_shift(
                        s[0], rest, lb, consDaysOff, ub, consDaysOff)

                    nw = int(s[4])
                    iui.add_min_max_weekend_in_4_weeks(
                        s[0], work, 0, 0, round(nw / n4), totalWeekends)

                    if int(s[5]) == 1:
                        iui.add_complete_weekend(s[0], complete)

            elif line.startswith("NURSES"):
                n = int(line.rsplit("=", 1)[-1])
                for j in range(n):
                    i += 1
                    s = lines[i].strip().split()
                    f_skills = set(iui.skills)
                    f_skills.difference_update(set(s[3:]))
                    contracts = [defC, s[1]] + ["N_"+sk for sk in f_skills]
                    iui.add_employee(s[0], contracts)

            i += 1

        return iui


def read_demand(path, instance, index, startDate, iui):
    s_file = os.path.join(path, instance, "WD-{}-{}.txt".format(instance, index))
    with open(s_file) as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if line.startswith("WEEK_DATA"):
                i += 1
                line = lines[i].strip()
                if not iui.name.startswith(line):
                    raise ValueError("Demand {} (for {}) is not for instance {}".format(index, line, iui.name))

            elif line.startswith("REQUIREMENTS"):
                i += 1
                line = lines[i].strip()
                while line:
                    line = line.split()
                    shift = line[0]
                    skill = line[1]
                    for j, l in enumerate(line[2:]):
                        d = startDate + timedelta(days=j)
                        lb, ub = get_bds(l)
                        iui.add_demand(d.isoformat(), shift, skill, lb, hard, lb, 0)
                        iui.add_demand(d.isoformat(), shift, skill, ub, underCoverage, ub, 0, 1)

                    i += 1
                    line = lines[i].strip() if i < len(lines) else ""

            elif line.startswith("SHIFT_OFF_REQUESTS"):
                n = int(line.rsplit("=", 1)[-1])
                for j in range(n):
                    i += 1
                    s = lines[i].strip().split()
                    d = startDate + timedelta(days=weekday_indexes[s[2]])
                    iui.add_preference_off(d.isoformat(), s[0], s[1], pref)

            i += 1


def read_history(path, instance, index, iui):
    s_file = os.path.join(path, instance, "H0-{}-{}.txt".format(instance, index))
    startDate = iui.start_date

    with open(s_file) as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if line.startswith("HISTORY"):
                i += 1
                ind, name = lines[i].strip().split()
                if not iui.name.startswith(name):
                    raise ValueError("History {} ({}) is not for instance {}".format(index, lines[i][:-1], iui.name))

            elif line.startswith("NURSE_HISTORY"):
                i += 1
                line = lines[i].strip()
                while line:
                    line = line.split()
                    name = line[0]
                    shift = line[3]
                    if shift == "None":
                        o_shift = iui.shifts[0]
                        n = int(line[6])
                        d = startDate + timedelta(days=-n-1)
                        iui.add_history(d.isoformat(), name, o_shift)
                    else:
                        o_shift = next(s for s in iui.shifts if s != shift)
                        n = int(line[5])
                        cons = int(line[4])
                        for j in range(0, n):
                            d = startDate + timedelta(days=j-n)
                            iui.add_history(d.isoformat(), name, shift if n - j <= cons else o_shift)

                    i += 1
                    line = lines[i].strip() if i < len(lines) else ""

            i += 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "Convert INRCII for ui format. You can either use the flag -n n005w4_0_1-2-3-3 or "
        "-i n005w4 -his 0 -d 1-2-3-3")
    parser.add_argument('-n', '--name', dest='name',
                        help='instance name, for example n005w4_0_1-2-3-3')
    parser.add_argument('-i', '--instance', dest='instance',
                        help='instance, for example n005w4')
    parser.add_argument('-p', '--path', dest='path',
                        help='directory path')
    parser.add_argument('-d', '--demand', dest='demand',
                        help='List of demand weeks separated by a dash (-), for example 1-2-3-3')
    parser.add_argument('-his', '--history', dest='history',
                        help='history index, for example 0')
    args = parser.parse_args()

    if args.name:
        n = args.name.split("_")
        args.instance = n[0]
        args.history = n[1]
        args.demand = n[2]
    else:
        args.name = "{}_{}_{}".format(
            args.instance, args.history, args.demand)

    iui = read_scenario(args.path, args.instance, args.name)

    startDate = iui.start_date
    for i in args.demand.split("-"):
        read_demand(args.path, args.instance, int(i), startDate, iui)
        startDate += timedelta(weeks=1)

    read_history(args.path, args.instance, int(args.history), iui)

    s_ui = iui.to_string()
    print(s_ui)
    f_ui = os.path.join(args.path, args.instance, "{}.txt".format(iui.name))
    with open(f_ui, "w") as f:
        f.write(s_ui)
        print("Instance {} has been saved to {}".format(iui.name, f_ui))
