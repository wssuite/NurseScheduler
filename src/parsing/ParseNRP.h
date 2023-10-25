/*
 * Copyright (C) 2020 Antoine Legrain, Jeremy Omer, and contributors.
 * All Rights Reserved.
 *
 * You may use, distribute and modify this code under the terms of the MIT
 * license.
 *
 * Please see the LICENSE file or visit https://opensource.org/licenses/MIT for
 * full license detail.
 */

#ifndef SRC_PARSING_PARSENRP_H_
#define SRC_PARSING_PARSENRP_H_

#include <string>

#include "tools/Tools.h"
#include "data/Scenario.h"

PScenario readNRPInstance(const string &fileName);

#endif  // SRC_PARSING_PARSENRP_H_
