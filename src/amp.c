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

#include <complex.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>

#include "digivoice.h"
#include "defines.h"
#include "sine.h"
#include "encode.h"
#include "mbest.h"
#include "amp.h"

static void rate_K_mbest_encode(float [], uint16_t []);
static void interp_para(float [], float [], float [], int, float [], int);
static void resample_const_rate_f(float [], MODEL *);
static void post_filter_amp(float []);
static void interp_Wo_v(float [], int [], bool [], float, bool);
static void resample_rate_L(MODEL *, int);
static void determine_phase(MODEL *);
static void amp_index_to_rate_K_vec(float [], int, int, int);

static float Amp_freqs_kHz[] = {    // Verified [srs]
    0.199816f,
    0.278224f,
    0.363464f,
    0.456131f,
    0.556873f,
    0.666393f,
    0.785457f,
    0.914895f,
    1.055613f,
    1.208592f,
    1.374901f,
    1.555703f,
    1.752259f,
    1.965942f,
    2.198245f,
    2.450789f,
    2.725340f,
    3.023815f,
    3.348299f,
    3.701056f
};

/* postfilter 20 * log10(freq / 0.3) */

static float Amp_pre[] = {    // Verified [srs]
    -3.529820f,
    -0.654534f,
    1.666803f,
    3.639367f,
    5.372698f,
    6.932183f,
    8.360023f,
    9.685000f,
    10.927670f,
    12.103168f,
    13.223004f,
    14.296108f,
    15.329541f,
    16.328987f,
    17.299097f,
    18.243692f,
    19.165989f,
    20.068680f,
    20.954060f,
    21.824089f
};

static float Amp_interpolated_surface_[N_MODELS][AMP_K];
static float Amp_prev_rate_K_vec_[AMP_K];
static float Amp_Wo_left;
static bool Amp_voicing_left;

/*
 * Quantize the model parameters into indexed bits
 * This is the data that will be transmitted across the medium
 * 
 * The indexed values are encoded into 16 bits:
 * [xxxx | yyyy yyyy yyyy ] Where x = #bits and y = the bits
 * 
 * index[0] = VQ magnitude1 (9 bits)
 * index[1] = VQ magnitude2 (9 bits)
 * index[2] = energy        (4 bits)
 * index[3] = pitch         (6 bits)
 */
void amp_model_to_index(uint16_t index[], MODEL *model) {
    float vec[AMP_K];
    float vec_no_mean[AMP_K];

    /* convert variable rate L to fixed rate K */

    resample_const_rate_f(vec, model);

    /* remove mean and two stage VQ */

    float sum = 0.0f;
    
    for (int k = 0; k < AMP_K; k++) {
        sum += vec[k];
    }

    float mean = sum / AMP_K;
    
    /* scalar quantize mean (effectively the frame energy) */

    index[2] = (4 << 12) | encode_energy(mean);     // energy 4 bits
    
    for (int k = 0; k < AMP_K; k++) {
        vec_no_mean[k] = vec[k] - mean;
    }

    rate_K_mbest_encode(vec_no_mean, index);    // VQ pair, each 9 bits

    /*
     * We steal the smallest Wo index to signal an unvoiced frame
     */

    if (model->voiced) {
        uint16_t pitch = encode_pitch(model->Wo);

        if (pitch == 0) {
            pitch = 1;
        }

        index[3] = (6 << 12) | pitch;       // pitch 6 bits
    } else {
        index[3] = (6 << 12) | 0;
    }
}

/*
 * Convert the quantized and indexed data bits back into a model
 */
void amp_index_to_models(MODEL models[], uint16_t index[]) {
    float vec_[AMP_K];

    int bits = index[0] >> 12;
    int mask = (1 << bits) - 1;
    int n1 = index[0] & mask;  // VQ1 Magnitude
    
    bits = index[1] >> 12;
    mask = (1 << bits) - 1;    
    int n2 = index[1] & mask;  // VQ2 Magnitude
    
    bits = index[2] >> 12;
    mask = (1 << bits) - 1;    
    int energy = index[2] & mask;  // Energy

    bits = index[3] >> 12;
    mask = (1 << bits) - 1;    
    int pitch = index[3] & mask;  // Pitch
    
    /* extract latest rate K vector */
    
    amp_index_to_rate_K_vec(vec_, n1, n2, energy);

    float Wo_right;
    bool voicing_right;

    /* decode latest Wo and voicing */

    if (pitch == 0) {
        /*
         * If pitch is zero, it is a code
         * to signal an unvoiced frame
         */
        Wo_right = TAU / 100.0f;
        voicing_right = false;
    } else {
        Wo_right = decode_pitch(pitch);
        voicing_right = true;
    }

    /* (linearly) interpolate 25Hz amplitude vectors back to 100Hz */

    float c;
    int i, k;

    for (i = 0, c = 1.0f; i < N_MODELS; i++, c -= 1.0f / N_MODELS) {
        for (k = 0; k < AMP_K; k++) {
            Amp_interpolated_surface_[i][k] = 
                    Amp_prev_rate_K_vec_[k] * c + vec_[k] * (1.0f - c);
        }
    }

    /* interpolate 25Hz v and Wo back to 100Hz */

    float aWo_[N_MODELS];
    bool avoicing_[N_MODELS];
    int aL_[N_MODELS];

    interp_Wo_v(aWo_, aL_, avoicing_, Wo_right, voicing_right);

    /* back to rate L amplitudes, synthesis phase for each frame */

    for (int i = 0; i < N_MODELS; i++) {
        models[i].Wo = aWo_[i];
        models[i].L = aL_[i];
        models[i].voiced = avoicing_[i];

        resample_rate_L(&models[i], i);
        determine_phase(&models[i]);
    }

    /* update memories for next time */

    for (int i = 0; i < AMP_K; i++) {
        Amp_prev_rate_K_vec_[i] = vec_[i];
    }

    Amp_Wo_left = Wo_right;
    Amp_voicing_left = voicing_right;
}

/*
 * A post filter is the key to the (relatively) high quality at such low bit rates.
 * The way it works is a little mysterious - and also a good topic for research.
 */
static void post_filter_amp(float vec[]) {
    float e_before = 0.0f;
    float e_after = 0.0f;

    for (int k = 0; k < AMP_K; k++) {
        vec[k] += Amp_pre[k];
        e_before += (powf(10.0f, 2.0f * vec[k] / 20.0f));
        
        vec[k] *= 1.5f;
        e_after += (powf(10.0f, 2.0f * vec[k] / 20.0f));
    }

    float gaindB = 10.0f * log10f(e_after / e_before);

    for (int k = 0; k < AMP_K; k++) {
        vec[k] -= gaindB;
        vec[k] -= Amp_pre[k];
    }
}

static void interp_para(float result[], float xp[], float yp[], int np, float x[], int n) {
    int k = 0;
    
    for (int i = 0; i < n; i++) {
        float xi = x[i];

        /* k is index into xp of where we start 3 points used to form parabola */

        while ((xp[k + 1] < xi) && (k < (np - 3)))
            k++;

        float x1 = xp[k];
        float y1 = yp[k];
        float x2 = xp[k + 1];
        float y2 = yp[k + 1];
        float x3 = xp[k + 2];
        float y3 = yp[k + 2];

        float a = ((y3 - y2) / (x3 - x2) - (y2 - y1) / (x2 - x1)) / (x3 - x1);
        float b = ((y3 - y2) / (x3 - x2) * (x2 - x1) + (y2 - y1) / (x2 - x1) * (x3 - x2)) / (x3 - x1);

        result[i] = a * (xi - x2) * (xi - x2) + b * (xi - x2) + y2;
    }
}

static void amp_index_to_rate_K_vec(float vec_[], int n1, int n2, int energy) {
    float vec_no_mean_[AMP_K];
    
    for (int k = 0; k < AMP_K; k++) {
        vec_no_mean_[k] = codebook1[AMP_K * n1 + k] + codebook2[AMP_K * n2 + k];
    }

    post_filter_amp(vec_no_mean_);
    
    float mean = decode_energy(energy);     // 4 bits

    for (int k = 0; k < AMP_K; k++) {
        vec_[k] = vec_no_mean_[k] + mean;
    }
}

static void resample_const_rate_f(float vec[], MODEL *model) {
    float amdB[MAX_AMP + 1];
    float rate_L_sample_freqs_kHz[MAX_AMP + 1];

    /* convert rate L=pi/Wo amplitude samples to fixed rate K */

    float amdB_peak = -100.0f;
    float tval = model->Wo * 4.0f / M_PI;
    
    for (int m = 1; m <= model->L; m++) {
        amdB[m] = 20.0f * log10f(model->A[m] + 1E-16f);

        if (amdB[m] > amdB_peak) {
            amdB_peak = amdB[m];
        }

        rate_L_sample_freqs_kHz[m] = (float) m * tval;
    }

    /* clip between peak and peak -50dB, to reduce dynamic range */

    for (int m = 1; m <= model->L; m++) {
        if (amdB[m] < (amdB_peak - 50.0f)) {
            amdB[m] = amdB_peak - 50.0f;
        }
    }

    interp_para(vec, &rate_L_sample_freqs_kHz[1], &amdB[1], model->L, Amp_freqs_kHz, AMP_K);
}

static void rate_K_mbest_encode(float vec_no_mean[], uint16_t index[]) {
    uint16_t entry[MBEST_STAGES];

    for (int i = 0; i < MBEST_STAGES; i++) {
        entry[i] = 0;
    }
    
    /* codebook is compiled for a fixed K */

    struct MBEST *mbest_stage1 = mbest_create();
    struct MBEST *mbest_stage2 = mbest_create();

    /* Stage 1 */

    mbest_search(codebook1, vec_no_mean, mbest_stage1, entry);

    /* Stage 2 */

    float target[AMP_K];
    int n1;

    for (int j = 0; j < MBEST_ENTRIES; j++) {
        entry[1] = n1 = mbest_stage1->list[j].index[0];

        for (int i = 0; i < AMP_K; i++) {
            target[i] = vec_no_mean[i] - codebook1[AMP_K * n1 + i];
        }
        
        mbest_search(codebook2, target, mbest_stage2, entry);
    }

    index[0] = (9 << 12) | mbest_stage2->list[0].index[1];
    index[1] = (9 << 12) | mbest_stage2->list[0].index[0];

    mbest_destroy(mbest_stage1);
    mbest_destroy(mbest_stage2);
}

static void interp_Wo_v(float Wo_[], int L_[], bool voicing_[], float Wo2, bool voicing_right) {
    int i;

    float tval = TAU / 100.0f;
    
    for (i = 0; i < N_MODELS; i++) {
        voicing_[i] = false;
    }

    if (!Amp_voicing_left && !voicing_right) {
        for (i = 0; i < N_MODELS; i++) {
            Wo_[i] = tval;
        }
    }

    if (Amp_voicing_left && !voicing_right) {
        Wo_[0] = Wo_[1] = Amp_Wo_left;
        Wo_[2] = Wo_[3] = tval;
        
        voicing_[0] = voicing_[1] = true;
    }

    if (!Amp_voicing_left && voicing_right) {
        Wo_[0] = Wo_[1] = tval;
        Wo_[2] = Wo_[3] = Wo2;
        
        voicing_[2] = voicing_[3] = true;
    }

    if (Amp_voicing_left && voicing_right) {
        float c = 1.0f;
        
        for (i = 0; i < N_MODELS; i++) {
            Wo_[i] = Amp_Wo_left * c + Wo2 * (1.0f - c);
            
            voicing_[i] = true;
            c -= 0.025f;        // 1 / N_MODELS
        }
    }

    for (i = 0; i < N_MODELS; i++) {
        L_[i] = floorf(M_PI / Wo_[i]);
    }
}

static void resample_rate_L(MODEL *model, int index) {
    float rate_K_vec_term[AMP_K + 2];
    float rate_K_sample_freqs_kHz_term[AMP_K + 2];
    float amdB[MAX_AMP + 1];
    float rate_L_sample_freqs_kHz[MAX_AMP + 1];

    /* init to zero in case we dump later for debug 0..80 */
    
    for (int m = 0; m <= MAX_AMP; m++) {
        model->A[m] = 0.0f;
    }

    /* terminate either end of the rate K vecs with 0dB points */

    rate_K_vec_term[0] = rate_K_vec_term[AMP_K + 1] = 0.0f;
    
    rate_K_sample_freqs_kHz_term[0] = 0.0f;
    rate_K_sample_freqs_kHz_term[AMP_K + 1] = 4.0f;

    for (int k = 0; k < AMP_K; k++) {
        rate_K_vec_term[k + 1] = Amp_interpolated_surface_[index][k];
        rate_K_sample_freqs_kHz_term[k + 1] = Amp_freqs_kHz[k];
    }

    float tval = model->Wo * 4.0f / M_PI;
    
    for (int m = 1; m <= model->L; m++) {
        rate_L_sample_freqs_kHz[m] = m * tval;
    }

    interp_para(&amdB[1], rate_K_sample_freqs_kHz_term, rate_K_vec_term, AMP_K + 2,
            &rate_L_sample_freqs_kHz[1], model->L);

    for (int m = 1; m <= model->L; m++) {
        model->A[m] = powf(10.0f, amdB[m] / 20.0f);
    }
}

static void determine_phase(MODEL *model) {
    float rate_L_sample_freqs_kHz[MAX_AMP + 1];
    float amdB[MAX_AMP + 1];
    float Gdbfk[NS];
    float sample_freqs_kHz[NS];
    float phase[NS];

    float tval = model->Wo * 4.0f / M_PI;
    
    for (int m = 1; m <= model->L; m++) {
        amdB[m] = 20.0f * log10f(model->A[m]);
        rate_L_sample_freqs_kHz[m] = (float) m * tval;
    }

    for (int i = 0; i < NS; i++) {
        sample_freqs_kHz[i] = 8.0f * (float) i / PHASE_FFT_SIZE;
    }

    interp_para(Gdbfk, &rate_L_sample_freqs_kHz[1], &amdB[1], model->L, sample_freqs_kHz, NS);

    mag_to_phase(phase, Gdbfk);

    tval = model->Wo * (float) PHASE_FFT_SIZE / TAU;
    
    for (int m = 1; m <= model->L; m++) {
        int b = floorf(0.5f + (float) m * tval);
        model->H[m] = cmplx(phase[b]);
    }
}
