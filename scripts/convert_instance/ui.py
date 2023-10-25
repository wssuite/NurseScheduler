from datetime import date


def get_list_string(list_, sep = ","):
    s = "{}".format(len(list_))
    if list_:
        s += ",{}".format(sep.join(list_))
    return s


class Pattern:
    def __init__(self):
        self.days = []
        self.shifts = []

    def add(self, days, shifts):
        self.days.append(days)
        self.shifts.append(shifts)

    def add_shifts(self, shifts):
        self.days.append(["*"])
        self.shifts.append(shifts)

    def to_string(self):
        s = "{}".format(len(self.days))
        for i, d in enumerate(self.days):
            s += ",{};{}".format("|".join(d), "|".join(self.shifts[i]))
        return s


class ui:
    def __init__(self, name, start, end):
        self.name = name
        self.__period = (start, end)
        self.__skills = []
        self.__shifts = {}
        self.__shift_types = {}
        self.__shift_groups = []
        self.__contracts = {}
        self.__contract_groups = []
        self.__employees = []
        self.__employees_index = {}
        self.__demands = []
        self.__preferences = []
        self.__history = []

    @property
    def n_shifts(self):
        return len(self.__shifts)

    @property
    def shifts(self):
        return list(self.__shifts.keys())

    @property
    def skills(self):
        return self.__skills

    @property
    def start_date(self):
        return date.fromisoformat(self.__period[0])

    def add_skill(self, skill):
        self.__skills.append(skill)

    def add_shift(self, shift, start, end):
        self.__shifts[shift] = (start, end)

    def add_shift_type(self, shift_type, shifts: list):
        self.__shift_types[shift_type] = \
            shift_type + "," + get_list_string(shifts)

    def add_shift_group(self, name, shifts: list, shift_types: list):
        s = name
        s += "," + get_list_string(shifts)
        s += "," + get_list_string(shift_types)
        self.__shift_groups.append(s)

    def add_rest_group(self):
        for g in self.__shift_groups:
            if g.statsWith("Rest"):
                return
        self.__shift_groups.append("Rest,0,0")

    def add_work_group(self):
        for g in self.__shift_groups:
            if g.startswith("Work"):
                return
        self.add_shift_group(
            "Work", list(self.__shifts), list(self.__shift_types))

    def add_constraint(self, name, ctr):
        if name in self.__contracts:
            self.__contracts[name].append(ctr)
        else:
            self.__contracts[name] = [ctr]

    def add_min_max_cons_shift(self, contract, shift_type, lb, clb, ub, cub):
        s = "MinMaxConsecutiveShiftType,{},{},{},{},{}".format(
            lb, clb, ub, cub, shift_type)
        self.add_constraint(contract, s)

    def add_min_max_assignment_in_4_weeks(self, contract, shift, lb, clb, ub, cub):
        s = "MinMaxNumAssignmentsInFourWeeks,{},{},{},{},{}".format(
            lb, clb, ub, cub, shift)
        self.add_constraint(contract, s)

    def add_min_max_weekend_in_4_weeks(self, contract, shift, lb, clb, ub, cub):
        s = "TotalWeekendsInFourWeeks,{},{},{},{},{}".format(
            lb, clb, ub, cub, shift)
        self.add_constraint(contract, s)

    def add_complete_weekend(self, contract, cost):
        self.add_constraint(contract, "CompleteWeekends,{}".format(cost))

    def forbid_shift(self, contract, shift, cost):
        self.add_constraint(contract, "UnwantedShift,{},{}".format(cost, shift))

    def forbid_skills(self, contract, skills, cost):
        self.add_constraint(contract, "UnwantedSkills,{},{}".format(cost, get_list_string(skills)))

    def add_unwanted_pattern(self, contract, pattern, cost):
        self.add_constraint(contract, "UnwantedPatterns,{},{}".format(cost, pattern.to_string()))

    def add_contract_group(self, name, contracts):
        self.__contract_groups[name] = name + "," + get_list_string(contracts)

    def add_employee(self, name, contracts, contract_groups=[]):
        self.__employees_index[name] = len(self.__employees)
        s = "{},{},{},{}".format(
            len(self.__employees), name,
            get_list_string(contracts),
            get_list_string(contract_groups))
        self.__employees.append(s)

    def add_demand(self, date, shift, skill, lb, clb, ub, cub, index=0):
        while len(self.__demands) <= index:
            self.__demands.append([])
        self.__demands[index].append("{},{},{},{},{},{},{}".format(
            date, shift, skill, lb, clb, ub, cub))

    def add_preference_on(self, date, employee, shift, cost):
        self.add_preference(date, employee, True, shift, cost)

    def add_preference_off(self, date, employee, shift, cost):
        self.add_preference(date, employee, False, shift, cost)

    def add_preference(self, date, employee, on, shift, cost):
        i = self.__employees_index[employee]
        self.__preferences.append("{},{},{},{},{}".format(
            date, i, "ON" if on else "OFF", shift, cost))

    def add_history(self, date, employee, shift):
        i = self.__employees_index[employee]
        self.__history.append("{},{},{}".format(date, i, shift))

    def to_string(self):
        e = "END\n"

        s = "SCHEDULING_PERIOD\n"
        s += "{},1.0,{},{}\n".format(self.name, self.__period[0], self.__period[1])
        s += e

        s += "SKILLS\n"
        for sk in self.__skills:
            s += sk+"\n"
        s += e

        s += "SHIFTS\n"
        for k, v in self.__shifts.items():
            s += "{},{},{}\n".format(k, v[0], v[1])
        s += e

        s += "SHIFT_TYPES\n"
        for v in self.__shift_types.values():
            s += "{}\n".format(v)
        s += e

        s += "SHIFT_GROUPS\n"
        self.add_rest_group()
        self.add_work_group()
        for v in self.__shift_groups:
            s += "{}\n".format(v)
        s += e

        s += "CONTRACTS\n"
        for k, v in self.__contracts.items():
            s += "{\ncontractName,"+k+"\nconstraints\n"
            for c in v:
                s += c + "\n"
            s += "}\n"
        s += e

        if self.__contract_groups:
            s += "CONTRACT_GROUPS\n"
            for v in self.__contract_groups.values():
                s += v + "\n"
            s += e

        s += "EMPLOYEES\n"
        for v in self.__employees:
            s += v + "\n"
        s += e

        for demand in self.__demands:
            s += "HOSPITAL_DEMAND\n"
            for v in demand:
                s += v + "\n"
            s += e

        s += "PREFERENCES\n"
        for v in self.__preferences:
            s += v + "\n"
        s += e

        if self.__history:
            s += "HISTORY\n"
            for v in self.__history:
                s += v + "\n"
            s += e

        return s
