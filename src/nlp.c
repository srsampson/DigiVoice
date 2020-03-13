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
#include <complex.h>

#include "defines.h"
#include "nlp.h"
#include "fft.h"

static float cnormf(complex float);
static int post_process_sub_multiples(float [], float, int);

static const float Nlp_cosw[] = {   // Verified [srs]
    0.000000f,
    0.002485f,
    0.009914f,
    0.022214f,
    0.039262f,
    0.060889f,
    0.086881f,
    0.116978f,
    0.150882f,
    0.188255f,
    0.228727f,
    0.271895f,
    0.317330f,
    0.364580f,
    0.413176f,
    0.462635f,
    0.512465f,
    0.562172f,
    0.611260f,
    0.659243f,
    0.705644f,
    0.750000f,
    0.791872f,
    0.830843f,
    0.866526f,
    0.898566f,
    0.926645f,
    0.950484f,
    0.969846f,
    0.984539f,
    0.994415f,
    0.999378f,
    0.999378f,
    0.994415f,
    0.984539f,
    0.969846f,
    0.950484f,
    0.926645f,
    0.898566f,
    0.866526f,
    0.830843f,
    0.791872f,
    0.750000f,
    0.705644f,
    0.659243f,
    0.611261f,
    0.562172f,
    0.512465f,
    0.462635f,
    0.413176f,
    0.364580f,
    0.317329f,
    0.271895f,
    0.228727f,
    0.188255f,
    0.150882f,
    0.116978f,
    0.086881f,
    0.060889f,
    0.039262f,
    0.022214f,
    0.009914f,
    0.002485f,
    0.000000f
};

/* 48 tap 600Hz low pass FIR filter coefficients, size 48 */

static const float Nlp_fir[] = {    // Verified [srs]
    -0.001082f,
    -0.001101f,
    -0.000928f,
    -0.000423f,
    0.000550f,
    0.002003f,
    0.003706f,
    0.005145f,
    0.005592f,
    0.004304f,
    0.000803f,
    -0.004820f,
    -0.011706f,
    -0.018199f,
    -0.022065f,
    -0.020921f,
    -0.012809f,
    0.003220f,
    0.026684f,
    0.055521f,
    0.086306f,
    0.114802f,
    0.136742f,
    0.148676f,
    0.148676f,
    0.136742f,
    0.114802f,
    0.086306f,
    0.055521f,
    0.026684f,
    0.003220f,
    -0.012809f,
    -0.020921f,
    -0.022065f,
    -0.018199f,
    -0.011706f,
    -0.004820f,
    0.000803f,
    0.004304f,
    0.005592f,
    0.005145f,
    0.003706f,
    0.002003f,
    0.000550f,
    -0.000423f,
    -0.000928f,
    -0.001101f,
    -0.001082f
};

static float Nlp_sq[M_PITCH];       // 320
static float Nlp_mem_x;
static float Nlp_mem_y;
static float Nlp_mem_fir[NLP_NTAP]; // 48
static int Nlp_prev_f0;

static fft_cfg Nlp_fft_cfg;

static float cnormf(complex float val) {
    float realf = crealf(val);
    float imagf = cimagf(val);

    return realf * realf + imagf * imagf;
}

int nlp_create() {
    if ((Nlp_fft_cfg = fft_alloc(FFT_SIZE, 0, NULL, NULL)) == NULL) {
        return -1;
    }

    return 0;
}

void nlp_destroy() {
    free(Nlp_fft_cfg);
}

int nlp(float Sn[]) {
    complex float Fw[FFT_SIZE];
    float fw[FFT_SIZE];
    
    for (int i = 0; i < FFT_SIZE; i++) {
	Fw[i] = 0.0f;
    }

    /* Square, notch filter at DC, and LP filter vector */

    for (int i = (M_PITCH - N_SAMP); i < M_PITCH; i++) { /* square last 80 speech samples */
        Nlp_sq[i] = (Sn[i] * Sn[i]);
    }

    for (int i = (M_PITCH - N_SAMP); i < M_PITCH; i++) { /* notch filter at DC */
        float notch = (Nlp_sq[i] - Nlp_mem_x) + (COEFF * Nlp_mem_y);

        Nlp_mem_x = Nlp_sq[i];
        Nlp_mem_y = notch;

        Nlp_sq[i] = notch + 1.0f;
    }

    for (int i = (M_PITCH - N_SAMP); i < M_PITCH; i++) { /* FIR filter vector */
        for (int j = 0; j < NLP_NTAP - 1; j++)
            Nlp_mem_fir[j] = Nlp_mem_fir[j + 1];

        Nlp_mem_fir[NLP_NTAP - 1] = Nlp_sq[i];

        Nlp_sq[i] = 0.0f;

        for (int j = 0; j < NLP_NTAP; j++)
            Nlp_sq[i] += (Nlp_mem_fir[j] * Nlp_fir[j]);
    }

    /* Decimate and FFT */

    for (int i = 0; i < (M_PITCH / DEC); i++) {
        Fw[i] = Nlp_sq[DEC * i] * Nlp_cosw[i];
    }

    fft(Nlp_fft_cfg, Fw, Fw);

    for (int i = 0; i < FFT_SIZE; i++) {
        fw[i] = cnormf(Fw[i]);
    }

    /* find global peak over 16..128 FFT filters */

    float gmax = 0.0f;
    int gmax_bin = FFT_SIZE * DEC / P_MAX;

    for (int i = (FFT_SIZE * DEC) / P_MAX; i <= (FFT_SIZE * DEC) / P_MIN; i++) {
        if (fw[i] > gmax) {
            gmax = fw[i];
            gmax_bin = i;
        }
    }

    /* Save as previous on next pass */
    
    Nlp_prev_f0 = post_process_sub_multiples(fw, gmax, gmax_bin);
    
    /* Shift samples in buffer to make room for new samples */

    for (int i = 0; i < (M_PITCH - N_SAMP); i++) {
        Nlp_sq[i] = Nlp_sq[N_SAMP + i];
    }

    /* return pitch */

    return FS / Nlp_prev_f0;
}

static int post_process_sub_multiples(float fw[], float gmax, int gmax_bin) {
    float thresh;
    int mult = 2;
    int cmax_bin = gmax_bin;
    int prev_f0_bin = Nlp_prev_f0 * (FFT_SIZE * DEC) / FS;

    while ((gmax_bin / mult) >= MIN_BIN) {

        int b = gmax_bin / mult; /* determine search interval */
        int bmin = 0.8f * b;
        int bmax = 1.2f * b;

        if (bmin < MIN_BIN) {
            bmin = MIN_BIN;
        }

        /* lower threshold to favor previous frames pitch estimate,
            this is a form of pitch tracking */

        if ((prev_f0_bin > bmin) && (prev_f0_bin < bmax)) {
            thresh = CNLP * gmax * 0.5f;
        } else {
            thresh = CNLP * gmax;
        }

        float lmax = 0.0f;
        int lmax_bin = bmin;

        for (int i = bmin; i <= bmax; i++) { /* look for maximum in interval */
            if (fw[i] > lmax) {
                lmax = fw[i];
                lmax_bin = i;
            }
        }

        if (lmax > thresh) {
            if ((lmax > fw[lmax_bin - 1]) && (lmax > fw[lmax_bin + 1])) {
                cmax_bin = lmax_bin;
            }
        }

        mult++;
    }

    return cmax_bin * (FS / (FFT_SIZE * DEC));
}
