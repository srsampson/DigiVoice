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

#include <math.h>

#include "encode.h"

static const float Energy_table[] = {   // Verified [srs]
    10.0f,
    12.5f,
    15.0f,
    17.5f,
    20.0f,
    22.5f,
    25.0f,
    27.5f,
    30.0f,
    32.5f,
    35.0f,
    37.5f,
    40.0f,
    42.5f,
    45.0f,
    47.5f
};

static const float Pitch_table[] = {    // Verified [srs]
    0.039270f,
    0.040567f,
    0.041907f,
    0.043290f,
    0.044720f,
    0.046197f,
    0.047723f,
    0.049299f,
    0.050927f,
    0.052609f,
    0.054346f,
    0.056141f,
    0.057995f,
    0.059910f,
    0.061889f,
    0.063932f,
    0.066044f,
    0.068225f,
    0.070478f,
    0.072806f,
    0.075210f,
    0.077694f,
    0.080260f,
    0.082910f,
    0.085648f,
    0.088477f,
    0.091399f,
    0.094417f,
    0.097535f,
    0.100756f,
    0.104084f,
    0.107521f,
    0.111072f,
    0.114740f,
    0.118529f,
    0.122444f,
    0.126488f,
    0.130665f,
    0.134980f,
    0.139438f,
    0.144043f,
    0.148800f,
    0.153714f,
    0.158790f,
    0.164034f,
    0.169451f,
    0.175047f,
    0.180828f,
    0.186800f,
    0.192969f,
    0.199342f,
    0.205925f,
    0.212726f,
    0.219751f,
    0.227008f,
    0.234505f,
    0.242250f,
    0.250250f,
    0.258515f,
    0.267052f,
    0.275871f,
    0.284982f,
    0.294394f,
    0.304116f
};

uint16_t encode_energy(float energy) {
    uint16_t bestindex = 0;
    float besterror = 1E32f;

    for (int i = 0; i < ENERGY_M; i++) {
        float diff = Energy_table[i] - energy;
        float error = (diff * diff);

        if (error < besterror) {
            besterror = error;
            bestindex = i;
        }
    }

    return bestindex & 0x0F;   // 4 bits
}

float decode_energy(int energy) {
    return Energy_table[energy];   // 4 bits
}

uint16_t encode_pitch(float wo) {
    uint16_t index = (uint16_t) floorf(WO_LEVELS *
                ((log10f(wo) - log10f(WO_MIN)) / WO_DIFF) + 0.5f);

    if (index < 0) {
        index = 0;
    } else if (index > (WO_LEVELS - 1)) {
        index = WO_LEVELS - 1;
    }

    return index & 0x3F;    // 6 bits
}

float decode_pitch(int pitch) {
    return Pitch_table[pitch];    // 6 bits
}
