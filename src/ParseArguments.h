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

#ifndef SRC_PARSEARGUMENTS_H_
#define SRC_PARSEARGUMENTS_H_

#include "tools/InputPaths.h"

// Read the arguments in noncompact/compact format
InputPaths *readArguments(int argc, char **argv);

#endif  // SRC_PARSEARGUMENTS_H_
