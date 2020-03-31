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

#include <stdbool.h>
#include <complex.h>
#include <string.h>
#include <math.h>

#include "defines.h"
#include "sine.h"
#include "nlp.h"
#include "encode.h"
#include "digivoice.h"
#include "mbest.h"
#include "amp.h"

char copyright1[] = "Copyright (C) 1993-2019 David Rowe, All rights reserved";

static MODEL Encode_model;
static MODEL Decode_models[N_MODELS];

int codec_create() {
    if (sine_create() != 0) {
        return -1;
    }

    if (nlp_create() != 0) {
        return -2;
    }

    return 0;
}

void codec_destroy() {
    sine_destroy();
    nlp_destroy();
}

int codec_indexes_per_frame() {
    return 4;
}

int codec_samples_per_frame() {
    return 320;
}

/*
 * Encodes frames of 320 samples of 15-bit + sign PCM
 * speech into an array index of bits.
 * 
 * 40ms segments, or at a 25 Hz rate
 * 
 * The indexed values are encoded into 16 bits:
 * [xxxx | yyyy yyyy yyyy ] Where x = #bits and y = the bits
 *
 * index[0] = VQ magnitude1 (9 bits)
 * index[1] = VQ magnitude2 (9 bits)
 * index[2] = energy        (4 bits)
 * index[3] = pitch         (6 bits)
 */
void codec_encode(uint16_t index[], int16_t speech[]) {
    /*
     * Process each 10 ms segment and update model
     * Only last model gets used going forward
     */
    for (int i = 0; i < N_MODELS; i++) {
        analyze_one_frame(&Encode_model, &speech[N_SAMP * i]);
    }

    /*
     * Convert the model into the indexed bits
     */
    amp_model_to_index(index, &Encode_model);
}

/*
 * Decodes array of indexed bits into 320 samples of speech (40ms)
 */
void codec_decode(int16_t speech[], uint16_t index[]) {

    amp_index_to_models(Decode_models, index);

    for (int i = 0; i < N_MODELS; i++) {
        synthesize_one_frame(&speech[N_SAMP * i], &Decode_models[i]);
    }
}

/*
 * Decodes energy value from encoded bits
 * Jeroen Vreeken, 2017
 */
float codec_get_energy(uint16_t index[]) {   
    int bits = index[2] >> 12;
    int mask = (1 << bits) - 1;    
    int energy = index[2] & mask;  // Energy
    
    bits = index[3] >> 12;
    mask = (1 << bits) - 1;    
    int pitch = index[3] & mask;  // Pitch
    
    float mean = decode_energy(energy) - 10.0f;

    /* Decrease mean if unvoiced */

    if (pitch == 0) // pitch == 0 means unvoiced
        mean -= 10.0f;

    return powf(10.0f, mean / 10.0f);
}
