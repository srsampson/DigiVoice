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

int encode_create(void);
int decode_create(void);

void encode_destroy(void);
void decode_destroy(void);

void codec_encode(uint16_t [], int16_t []);

void codec_decode(int16_t [], uint16_t []);
float codec_decode_energy(uint16_t []);

int codec_indexes_per_frame(void);
int codec_samples_per_frame(void);

#ifdef __cplusplus
}
#endif
