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

#include <complex.h>
#include <stdint.h>
    
#include "digivoice.h"
#include "defines.h"
#include "fft.h"

#define AMP_K               20  /* rate K vector length                            */
#define AMP_M               512 /* number of elements in codebook                  */
#define MBEST_ENTRIES       5   /* how many candidates we keep for each VQ stage   */
    
void amp_model_to_index(uint16_t index[], MODEL *model);
void amp_index_to_models(MODEL [], uint16_t []);

#ifdef __cplusplus
}
#endif
