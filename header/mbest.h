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
    
#define MBEST_STAGES 4

struct MBEST_LIST {
    uint16_t index[MBEST_STAGES];
    float error;
};

struct MBEST {
    struct MBEST_LIST *list;
};

struct MBEST *mbest_create(void);
void mbest_destroy(struct MBEST *);
void mbest_search(const float *, float [], struct MBEST *, uint16_t []);

#ifdef __cplusplus
}
#endif