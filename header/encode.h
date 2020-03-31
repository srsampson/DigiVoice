/*
 * Copyright (C) 1993-2019 David Rowe
 *
 * All rights reserved
 *
 * Modified March 2020 by Steve Sampson
 * 
 * Licensed under GNU LGPL V2.1
 * See LICENSE file for information
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
    
#include "defines.h"

#define WO_MIN      (TAU / P_MAX)
#define WO_MAX      (TAU / P_MIN)
#define WO_DIFF     (log10f(WO_MAX) - log10f(WO_MIN))
#define WO_LEVELS   (1 << 6)
    
#define ENERGY_M    16

uint16_t encode_energy(float);
float decode_energy(int);
uint16_t encode_pitch(float);
float decode_pitch(int);

#ifdef __cplusplus
}
#endif
