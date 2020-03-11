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
#include <math.h>
#include <complex.h>
#include <string.h>

#include "defines.h"
#include "sine.h"
#include "nlp.h"
#include "fft.h"

static void hs_pitch_refinement(MODEL *, float, float, float);
static void synthesize(MODEL *);
static void two_stage_pitch_refinement(MODEL *);
static void estimate_amplitudes(MODEL *);
static void est_voicing_mbe(MODEL *);
static void postfilter(MODEL *);
static void phase_synth_zero_order(MODEL *);
static int codec2_rand(void);
static float cnormf(complex float);

/* Window, size 160 */

static const float Parzen[] = { // Verified [srs]
    0.000000f,
    0.012500f,
    0.025000f,
    0.037500f,
    0.050000f,
    0.062500f,
    0.075000f,
    0.087500f,
    0.100000f,
    0.112500f,
    0.125000f,
    0.137500f,
    0.150000f,
    0.162500f,
    0.175000f,
    0.187500f,
    0.200000f,
    0.212500f,
    0.225000f,
    0.237500f,
    0.250000f,
    0.262500f,
    0.275000f,
    0.287500f,
    0.300000f,
    0.312500f,
    0.325000f,
    0.337500f,
    0.350000f,
    0.362500f,
    0.375000f,
    0.387500f,
    0.400000f,
    0.412500f,
    0.425000f,
    0.437500f,
    0.450000f,
    0.462500f,
    0.475000f,
    0.487500f,
    0.500000f,
    0.512500f,
    0.525000f,
    0.537500f,
    0.550000f,
    0.562500f,
    0.575000f,
    0.587500f,
    0.600000f,
    0.612500f,
    0.625000f,
    0.637500f,
    0.650000f,
    0.662500f,
    0.675000f,
    0.687500f,
    0.700000f,
    0.712500f,
    0.725000f,
    0.737500f,
    0.750000f,
    0.762500f,
    0.775000f,
    0.787500f,
    0.800000f,
    0.812500f,
    0.825000f,
    0.837499f,
    0.849999f,
    0.862499f,
    0.874999f,
    0.887499f,
    0.899999f,
    0.912499f,
    0.924999f,
    0.937499f,
    0.949999f,
    0.962499f,
    0.974999f,
    0.987499f,
    1.000000f,
    0.987500f,
    0.975000f,
    0.962500f,
    0.950000f,
    0.937500f,
    0.925000f,
    0.912500f,
    0.900000f,
    0.887500f,
    0.875000f,
    0.862500f,
    0.850000f,
    0.837500f,
    0.825000f,
    0.812500f,
    0.800000f,
    0.787500f,
    0.775000f,
    0.762500f,
    0.750000f,
    0.737500f,
    0.725000f,
    0.712500f,
    0.700000f,
    0.687500f,
    0.675000f,
    0.662500f,
    0.650000f,
    0.637500f,
    0.625000f,
    0.612500f,
    0.600000f,
    0.587500f,
    0.575000f,
    0.562500f,
    0.550000f,
    0.537500f,
    0.525000f,
    0.512500f,
    0.500000f,
    0.487500f,
    0.475001f,
    0.462501f,
    0.450001f,
    0.437501f,
    0.425001f,
    0.412501f,
    0.400001f,
    0.387501f,
    0.375001f,
    0.362501f,
    0.350001f,
    0.337501f,
    0.325001f,
    0.312501f,
    0.300001f,
    0.287501f,
    0.275001f,
    0.262501f,
    0.250001f,
    0.237501f,
    0.225001f,
    0.212501f,
    0.200001f,
    0.187501f,
    0.175001f,
    0.162501f,
    0.150001f,
    0.137501f,
    0.125001f,
    0.112501f,
    0.100001f,
    0.087501f,
    0.075001f,
    0.062501f,
    0.050001f,
    0.037501f,
    0.025001f,
    0.012501f
};

/* this is the old W[] array impulse waveform, size 512 */

static const float Hamming[] = {    // Verified [srs]
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000001f,
    0.000000f,
    -0.000001f,
    0.000000f,
    0.000001f,
    0.000000f,
    0.000000f,
    0.000001f,
    0.000000f,
    -0.000001f,
    0.000000f,
    0.000001f,
    0.000000f,
    -0.000001f,
    0.000000f,
    0.000001f,
    0.000000f,
    -0.000001f,
    0.000001f,
    0.000000f,
    -0.000001f,
    0.000000f,
    0.000001f,
    0.000000f,
    -0.000001f,
    0.000000f,
    0.000001f,
    -0.000001f,
    -0.000001f,
    0.000001f,
    0.000001f,
    -0.000001f,
    -0.000001f,
    0.000001f,
    0.000000f,
    -0.000001f,
    0.000000f,
    0.000002f,
    -0.000001f,
    -0.000001f,
    0.000001f,
    0.000001f,
    -0.000002f,
    -0.000001f,
    0.000002f,
    0.000000f,
    -0.000002f,
    0.000000f,
    0.000002f,
    -0.000001f,
    -0.000002f,
    0.000001f,
    0.000002f,
    -0.000002f,
    -0.000002f,
    0.000003f,
    0.000001f,
    -0.000003f,
    0.000000f,
    0.000004f,
    -0.000001f,
    -0.000004f,
    0.000002f,
    0.000003f,
    -0.000003f,
    -0.000003f,
    0.000004f,
    0.000002f,
    -0.000005f,
    -0.000001f,
    0.000006f,
    -0.000001f,
    -0.000007f,
    0.000003f,
    0.000006f,
    -0.000005f,
    -0.000006f,
    0.000007f,
    0.000004f,
    -0.000010f,
    -0.000002f,
    0.000012f,
    -0.000001f,
    -0.000013f,
    0.000005f,
    0.000013f,
    -0.000009f,
    -0.000013f,
    0.000014f,
    0.000011f,
    -0.000020f,
    -0.000007f,
    0.000025f,
    0.000000f,
    -0.000030f,
    0.000009f,
    0.000034f,
    -0.000020f,
    -0.000035f,
    0.000035f,
    0.000033f,
    -0.000053f,
    -0.000025f,
    0.000075f,
    0.000009f,
    -0.000099f,
    0.000019f,
    0.000124f,
    -0.000064f,
    -0.000148f,
    0.000135f,
    0.000163f,
    -0.000246f,
    -0.000158f,
    0.000421f,
    0.000102f,
    -0.000708f,
    0.000079f,
    0.001208f,
    -0.000597f,
    -0.002176f,
    0.002195f,
    0.004429f,
    -0.008645f,
    -0.012196f,
    0.065359f,
    0.262390f,
    0.495616f,
    0.601647f,
    0.495616f,
    0.262390f,
    0.065359f,
    -0.012196f,
    -0.008645f,
    0.004429f,
    0.002195f,
    -0.002176f,
    -0.000597f,
    0.001208f,
    0.000079f,
    -0.000708f,
    0.000102f,
    0.000421f,
    -0.000158f,
    -0.000246f,
    0.000163f,
    0.000135f,
    -0.000148f,
    -0.000064f,
    0.000124f,
    0.000019f,
    -0.000099f,
    0.000009f,
    0.000075f,
    -0.000025f,
    -0.000053f,
    0.000033f,
    0.000035f,
    -0.000035f,
    -0.000020f,
    0.000034f,
    0.000009f,
    -0.000030f,
    0.000000f,
    0.000025f,
    -0.000007f,
    -0.000020f,
    0.000011f,
    0.000014f,
    -0.000013f,
    -0.000009f,
    0.000013f,
    0.000005f,
    -0.000013f,
    -0.000001f,
    0.000012f,
    -0.000002f,
    -0.000010f,
    0.000004f,
    0.000007f,
    -0.000006f,
    -0.000005f,
    0.000006f,
    0.000003f,
    -0.000007f,
    -0.000001f,
    0.000006f,
    -0.000001f,
    -0.000005f,
    0.000002f,
    0.000004f,
    -0.000003f,
    -0.000003f,
    0.000003f,
    0.000002f,
    -0.000004f,
    -0.000001f,
    0.000004f,
    0.000000f,
    -0.000003f,
    0.000001f,
    0.000003f,
    -0.000002f,
    -0.000002f,
    0.000002f,
    0.000001f,
    -0.000002f,
    -0.000001f,
    0.000002f,
    0.000000f,
    -0.000002f,
    0.000000f,
    0.000002f,
    -0.000001f,
    -0.000002f,
    0.000001f,
    0.000001f,
    -0.000001f,
    -0.000001f,
    0.000002f,
    0.000000f,
    -0.000001f,
    0.000000f,
    0.000001f,
    -0.000001f,
    -0.000001f,
    0.000001f,
    0.000001f,
    -0.000001f,
    -0.000001f,
    0.000001f,
    0.000000f,
    -0.000001f,
    0.000000f,
    0.000001f,
    0.000000f,
    -0.000001f,
    0.000000f,
    0.000001f,
    -0.000001f,
    0.000000f,
    0.000001f,
    0.000000f,
    -0.000001f,
    0.000000f,
    0.000001f,
    0.000000f,
    -0.000001f,
    0.000000f,
    0.000001f,
    0.000000f,
    0.000000f,
    0.000001f,
    0.000000f,
    -0.000001f,
    0.000000f,
    0.000001f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f
};

/* This is the old w[] array, size = 320 */

static const float Hamming2[] = {   // Verified [srs]
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000001f,
    0.000002f,
    0.000005f,
    0.000009f,
    0.000014f,
    0.000020f,
    0.000027f,
    0.000035f,
    0.000045f,
    0.000055f,
    0.000067f,
    0.000079f,
    0.000093f,
    0.000107f,
    0.000123f,
    0.000140f,
    0.000158f,
    0.000177f,
    0.000196f,
    0.000217f,
    0.000239f,
    0.000262f,
    0.000286f,
    0.000311f,
    0.000336f,
    0.000363f,
    0.000391f,
    0.000419f,
    0.000448f,
    0.000479f,
    0.000510f,
    0.000542f,
    0.000575f,
    0.000608f,
    0.000643f,
    0.000678f,
    0.000714f,
    0.000750f,
    0.000788f,
    0.000826f,
    0.000865f,
    0.000904f,
    0.000944f,
    0.000985f,
    0.001026f,
    0.001068f,
    0.001110f,
    0.001153f,
    0.001197f,
    0.001241f,
    0.001285f,
    0.001330f,
    0.001376f,
    0.001421f,
    0.001468f,
    0.001514f,
    0.001561f,
    0.001608f,
    0.001655f,
    0.001703f,
    0.001751f,
    0.001799f,
    0.001847f,
    0.001896f,
    0.001944f,
    0.001993f,
    0.002042f,
    0.002091f,
    0.002140f,
    0.002189f,
    0.002238f,
    0.002286f,
    0.002335f,
    0.002384f,
    0.002433f,
    0.002481f,
    0.002529f,
    0.002577f,
    0.002625f,
    0.002673f,
    0.002720f,
    0.002768f,
    0.002814f,
    0.002861f,
    0.002907f,
    0.002953f,
    0.002998f,
    0.003043f,
    0.003087f,
    0.003131f,
    0.003175f,
    0.003218f,
    0.003260f,
    0.003302f,
    0.003344f,
    0.003384f,
    0.003424f,
    0.003464f,
    0.003503f,
    0.003541f,
    0.003578f,
    0.003615f,
    0.003651f,
    0.003686f,
    0.003720f,
    0.003754f,
    0.003787f,
    0.003819f,
    0.003850f,
    0.003880f,
    0.003909f,
    0.003938f,
    0.003965f,
    0.003992f,
    0.004018f,
    0.004043f,
    0.004066f,
    0.004089f,
    0.004111f,
    0.004132f,
    0.004152f,
    0.004171f,
    0.004188f,
    0.004205f,
    0.004221f,
    0.004236f,
    0.004249f,
    0.004262f,
    0.004273f,
    0.004284f,
    0.004293f,
    0.004301f,
    0.004309f,
    0.004315f,
    0.004320f,
    0.004323f,
    0.004326f,
    0.004328f,
    0.004328f,
    0.004328f,
    0.004326f,
    0.004323f,
    0.004320f,
    0.004315f,
    0.004309f,
    0.004301f,
    0.004293f,
    0.004284f,
    0.004273f,
    0.004262f,
    0.004249f,
    0.004236f,
    0.004221f,
    0.004205f,
    0.004188f,
    0.004171f,
    0.004152f,
    0.004132f,
    0.004111f,
    0.004089f,
    0.004066f,
    0.004043f,
    0.004018f,
    0.003992f,
    0.003965f,
    0.003938f,
    0.003909f,
    0.003880f,
    0.003850f,
    0.003819f,
    0.003787f,
    0.003754f,
    0.003720f,
    0.003686f,
    0.003651f,
    0.003615f,
    0.003578f,
    0.003541f,
    0.003503f,
    0.003464f,
    0.003424f,
    0.003384f,
    0.003344f,
    0.003302f,
    0.003260f,
    0.003218f,
    0.003175f,
    0.003131f,
    0.003087f,
    0.003043f,
    0.002998f,
    0.002953f,
    0.002907f,
    0.002861f,
    0.002814f,
    0.002768f,
    0.002720f,
    0.002673f,
    0.002625f,
    0.002577f,
    0.002529f,
    0.002481f,
    0.002433f,
    0.002384f,
    0.002335f,
    0.002286f,
    0.002238f,
    0.002189f,
    0.002140f,
    0.002091f,
    0.002042f,
    0.001993f,
    0.001944f,
    0.001896f,
    0.001847f,
    0.001799f,
    0.001751f,
    0.001703f,
    0.001655f,
    0.001608f,
    0.001561f,
    0.001514f,
    0.001468f,
    0.001421f,
    0.001376f,
    0.001330f,
    0.001285f,
    0.001241f,
    0.001197f,
    0.001153f,
    0.001110f,
    0.001068f,
    0.001026f,
    0.000985f,
    0.000944f,
    0.000904f,
    0.000865f,
    0.000826f,
    0.000788f,
    0.000750f,
    0.000714f,
    0.000678f,
    0.000643f,
    0.000608f,
    0.000575f,
    0.000542f,
    0.000510f,
    0.000479f,
    0.000448f,
    0.000419f,
    0.000391f,
    0.000363f,
    0.000336f,
    0.000311f,
    0.000286f,
    0.000262f,
    0.000239f,
    0.000217f,
    0.000196f,
    0.000177f,
    0.000158f,
    0.000140f,
    0.000123f,
    0.000107f,
    0.000093f,
    0.000079f,
    0.000067f,
    0.000055f,
    0.000045f,
    0.000035f,
    0.000027f,
    0.000020f,
    0.000014f,
    0.000009f,
    0.000005f,
    0.000002f,
    0.000001f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f,
    0.000000f
};

static fft_cfg Sine_fft_fwd_cfg;
static fft_cfg Sine_fft_inv_cfg;

static fftr_cfg Sine_fftr_fwd_cfg;
static fftr_cfg Sine_fftr_inv_cfg;

static complex float Sine_Sw[FFT_SIZE];
static float Sine_Sn[M_PITCH];
static float Sine_Sn_[N_SAMP * 2];
static float Sine_ex_phase;
static float Sine_bg_est;

static unsigned long Next = 1;

static float cnormf(complex float val) {
    float realf = crealf(val);
    float imagf = cimagf(val);

    return realf * realf + imagf * imagf;
}

int sine_create() {
    Sine_fftr_fwd_cfg = fftr_alloc(FFT_SIZE, 0, NULL, NULL);
    Sine_fftr_inv_cfg = fftr_alloc(FFT_SIZE, 1, NULL, NULL);

    Sine_fft_fwd_cfg  = fft_alloc(PHASE_FFT_SIZE, 0, NULL, NULL);
    Sine_fft_inv_cfg  = fft_alloc(PHASE_FFT_SIZE, 1, NULL, NULL);

    if (Sine_fftr_fwd_cfg == NULL || Sine_fftr_inv_cfg == NULL
        || Sine_fft_fwd_cfg == NULL || Sine_fft_inv_cfg == NULL) {
        return -1;
    }

    return 0;
}

void sine_destroy() {
    free(Sine_fft_inv_cfg);
    free(Sine_fft_fwd_cfg);
    free(Sine_fftr_inv_cfg);
    free(Sine_fftr_fwd_cfg);
}

static int codec2_rand() {
    Next = Next * 1103515245 + 12345;
    return((unsigned)(Next/65536) % 32768);
}

void analyze_one_frame(MODEL *encode_model, int16_t speech[]) {
    float sw[FFT_SIZE];

    for (int i = 0; i < FFT_SIZE; i++) {
	sw[i] = 0.0f;
    }

    /* process the new 80 samples of speech */

    for (int i = 0; i < (M_PITCH - N_SAMP); i++) {
        Sine_Sn[i] = Sine_Sn[N_SAMP + i];         // Left shift history 80 samples
    }

    for (int i = 0; i < N_SAMP; i++) {
        Sine_Sn[(M_PITCH - N_SAMP) + i] = (float) speech[i]; // Add new 80 samples to end
    }

    /* move 2nd half to start of FFT input vector */

    for (int i = 0; i < (NW / 2); i++) {
        int half = i + (M_PITCH / 2);

        sw[i] = Sine_Sn[half] * Hamming2[half];
    }

    /* move 1st half to end of FFT input vector */

    for (int i = 0; i < (NW / 2); i++) {
        int half = i + (M_PITCH / 2) - (NW / 2);

        sw[(FFT_SIZE - (NW / 2)) + i] = Sine_Sn[half] * Hamming2[half];
    }

    fftr(Sine_fftr_fwd_cfg, sw, Sine_Sw);

    encode_model->Wo = TAU / nlp(Sine_Sn);    // nlp returns pitch
    encode_model->L = M_PI / encode_model->Wo;

    /* fill-in the model values */

    two_stage_pitch_refinement(encode_model);   // operates on Sine_Sw
    estimate_amplitudes(encode_model);          // operates on Sine_Sw
    est_voicing_mbe(encode_model);              // operates on Sine_Sw
}

void synthesize_one_frame(int16_t speech[], MODEL *decode_model) {
    phase_synth_zero_order(decode_model);
    postfilter(decode_model);
    synthesize(decode_model);   // Populate Sine_Sn_

    /* Limit output audio */

    float max_sample = 0.0f;

    for (int i = 0; i < N_SAMP; i++) {
        if (Sine_Sn_[i] > max_sample) {
            max_sample = Sine_Sn_[i];
        }
    }

    float over = max_sample / 30000.0f;

    if (over > 1.0f) {
        float gain = 1.0f / (over * over);

        for (int i = 0; i < N_SAMP; i++) {
            Sine_Sn_[i] *= gain;
        }
    }
    
    // Mode 700c is a little weak over-all

    for (int i = 0; i < N_SAMP; i++) {
        Sine_Sn_[i] *= 1.5f;
    }
    
    for (int i = 0; i < N_SAMP; i++) {
        if (Sine_Sn_[i] > 32760.0f) {
            speech[i] = (int16_t) 32760;        // don't saturate it
        } else if (Sine_Sn_[i] < -32760.0f) {
            speech[i] = (int16_t) -32760;       // ditto
        } else {
            speech[i] = (int16_t) Sine_Sn_[i];
        }
    }
}

void mag_to_phase(float phase[], float mag[]) {
    complex float Sdb[PHASE_FFT_SIZE];
    complex float cf[PHASE_FFT_SIZE];

    Sdb[0] = mag[0];

    for (int i = 1; i < NS; i++) {                     // 1 - 64
        Sdb[i                 ] = mag[i] + 0.0f * I;   // 1 - 64
        Sdb[PHASE_FFT_SIZE - i] = mag[i] + 0.0f * I;   // 63 - 0
    }

    /* compute real cepstrum from log magnitude spectrum */

    complex float c[PHASE_FFT_SIZE];

    fft(Sine_fft_inv_cfg, Sdb, c);  // 128

    for (int i = 0; i < PHASE_FFT_SIZE; i++) {
        c[i] = c[i] / (float) PHASE_FFT_SIZE;
    }

    /* Fold cepstrum to reflect non-min-phase zeros inside unit circle */

    for (int i = 0; i < PHASE_FFT_SIZE; i++) {
	cf[i] = 0.0f;
    }

    cf[0] = c[0];

    for (int i = 1; i < (NS - 1); i++) {
        cf[i] = c[i] + c[PHASE_FFT_SIZE - i];
    }

    cf[NS - 1] = c[NS - 1];

    fft(Sine_fft_fwd_cfg, cf, cf);  // in-place 128

    for (int i = 0; i < NS; i++) {
        phase[i] = cimagf(cf[i]) / SCALE;
    }
}

static void phase_synth_zero_order(MODEL *model) {
    complex float ex[MAX_AMP + 1];

    Sine_ex_phase += ((model->Wo * N_SAMP) - (floorf(Sine_ex_phase / TAU + 0.5f) * TAU));

    for (int m = 1; m <= model->L; m++) {

        /* generate excitation */

        if (model->voiced == true) {
            ex[m] = cmplx((float) m * Sine_ex_phase);
        } else {
            ex[m] = cmplx(TAU * (float) codec2_rand() / CODEC2_RND_MAX);
        }

        /* filter using LPC filter */
        /* Note: H was populated during determine_phase() in Amp */

        ex[m] *= model->H[m];

        /* modify sinusoidal phase */

        model->phi[m] = atan2f(cimagf(ex[m]), crealf(ex[m]) + 1E-12f);
    }
}

static void two_stage_pitch_refinement(MODEL *model) {
    float tval = TAU / model->Wo; /* compute once for below */
    
    /* Coarse refinement */

    float pmax = tval + 5.0f;
    float pmin = tval - 5.0f;
    hs_pitch_refinement(model, pmin, pmax, 1.0f);

    tval = TAU / model->Wo; /* compute once for below */

    /* Fine refinement */

    pmax = tval + 1.0f;
    pmin = tval - 1.0f;
    hs_pitch_refinement(model, pmin, pmax, 0.25f);

    /* Limit range */

    if (model->Wo < (TAU / P_MAX)) {            // 0.039269875 (L = 80)
        model->Wo = (TAU / P_MAX);
    } else if (model->Wo > (TAU / P_MIN)) {     // 0.314159 (L = 10)
        model->Wo = (TAU / P_MIN);
    }

    model->L = floorf(M_PI / model->Wo);
    
    if (model->Wo * model->L >= FRACTPI) {
        model->L--;
    }
}

static void hs_pitch_refinement(MODEL *model, float pmin, float pmax, float pstep) {
    float pitch;

    model->L = M_PI / model->Wo; /* use initial pitch est. for L */

    float wom = model->Wo; /* Wo that maximizes E */
    float em = 0.0f; /* maximum energy */

    /* Determine harmonic sum for a range of Wo values */

    for (pitch = pmin; pitch <= pmax; pitch += pstep) {
        float e = 0.0f;
        float wo = TAU / pitch;

        /* Sum harmonic magnitudes */
        float tval = wo * ONE_ON_R;
        
        for (int m = 1; m <= model->L; m++) {
            int b = (int) ((float) m * tval + 0.5f);
            e += cnormf(Sine_Sw[b]);
        }

        /* Compare to see if this is a maximum */

        if (e > em) {
            em = e;
            wom = wo;
        }
    }

    model->Wo = wom;
}

static void estimate_amplitudes(MODEL *model) {
    float amp = model->Wo * ONE_ON_R;

    /* init to zero in case we dump later for debug 0..80 */
    
    for (int m = 0; m <= MAX_AMP; m++) {
        model->A[m] = 0.0f;
    }

    for (int m = 1; m <= model->L; m++) {
        int am = (int) ((m - 0.5f) * amp + 0.5f);
        int bm = (int) ((m + 0.5f) * amp + 0.5f);

        float den = 0.0f;

        for (int i = am; i < bm; i++) {
            den += cnormf(Sine_Sw[i]);
        }

        model->A[m] = sqrtf(den);
    }
}

static void est_voicing_mbe(MODEL *model) {
    float sig = 1E-4f;

    for (int l = 1; l <= (model->L / 4); l++) {
        sig += (model->A[l] * model->A[l]);
    }

    float wo = model->Wo * FFT_SIZE / TAU;
    float error = 1E-4f;

    /*
     * accumulated error between original and synthesized
     * Just test across the harmonics in the first 1000 Hz (L/4)
     */

    for (int l = 1; l <= (model->L / 4); l++) {
        complex float am = 0.0f;
        float den = 0.0f;

        int al = ceilf((l - 0.5f) * wo);
        int bl = ceilf((l + 0.5f) * wo);

        /* Estimate amplitude of harmonic assuming harmonic is totally voiced */

        int offset = (FFT_SIZE / 2) - l * wo + 0.5f;

        for (int m = al; m < bl; m++) {
            am += (Sine_Sw[m] * Hamming[offset + m]);
            den += (Hamming[offset + m] * Hamming[offset + m]);
        }

        am /= den;

        /* Determine error between estimated harmonic and original */

        for (int m = al; m < bl; m++) {
            error += cnormf(Sine_Sw[m] - (am * Hamming[offset + m]));
        }
    }

    float snr = 10.0f * log10f(sig / error);

    if (snr > V_THRESH)
        model->voiced = true;
    else
        model->voiced = false;

    float elow = 1E-4f;
    float ehigh = 1E-4f;

    for (int l = 1; l <= (model->L / 2); l++) {
        elow += (model->A[l] * model->A[l]);
    }

    for (int l = (model->L / 2); l <= model->L; l++) {
        ehigh += (model->A[l] * model->A[l]);
    }

    float eratio = 10.0f * log10f(elow / ehigh);

    if (model->voiced == false)
        if (eratio > 10.0f)
            model->voiced = true;

    if (model->voiced == true) {
        if (eratio < -10.0f)
            model->voiced = false;

        if ((eratio < -4.0f) && (model->Wo <= SIXTY))
            model->voiced = false;
    }
}

static void postfilter(MODEL *decode_model) {
    /* determine average energy across spectrum */

    float e = 1E-12f;

    for (int i = 1; i <= decode_model->L; i++)
        e += (decode_model->A[i] * decode_model->A[i]);

    e = 10.0f * log10f(e / decode_model->L);

    if ((e < BG_THRESH) && (decode_model->voiced == false))
        Sine_bg_est *= (1.0f - BG_BETA) + e * BG_BETA;

    float thresh = powf(10.0f, (Sine_bg_est + BG_MARGIN) / 20.0f);

    if (decode_model->voiced == true) {
        for (int i = 1; i <= decode_model->L; i++) {
            if (decode_model->A[i] < thresh) {
                decode_model->phi[i] = (TAU * (float) codec2_rand() / CODEC2_RND_MAX);
            }
        }
    }
}

/*
 * Synthesize a speech signal in the frequency domain from the sinusoidal
 * model parameters. Uses overlap-add with a trapezoidal window to smoothly
 * interpolate between frames.
 */
static void synthesize(MODEL *model) {
    complex float Sw_[FFT_SIZE / 2 + 1];
    float sw_[FFT_SIZE];

    /* Update memories */

    for (int i = 0; i < (N_SAMP - 1); i++) {
        Sine_Sn_[i] = Sine_Sn_[N_SAMP + i]; // Left shift history 80 samples
    }

    Sine_Sn_[N_SAMP - 1] = 0.0f;

    /* Now set up frequency domain synthesized speech */

    for (int i = 0; i < (FFT_SIZE / 2 + 1); i++) {
	Sw_[i] = 0.0f;
    }

    float wo = model->Wo * FFT_SIZE / TAU;

    for (int l = 1; l <= model->L; l++) {
        int b = (int) (l * wo + 0.5f);

        if (b > ((FFT_SIZE / 2) - 1)) {
            b = (FFT_SIZE / 2) - 1;
        }

        Sw_[b] = cmplx(model->phi[l]) * model->A[l];
    }

    /* Perform inverse real FFT */

    fftri(Sine_fftr_inv_cfg, Sw_, sw_);

    /* Overlap add to previous samples */

    for (int i = 0; i < (N_SAMP - 1); i++) {
        Sine_Sn_[i] += sw_[FFT_SIZE - N_SAMP + 1 + i] * Parzen[i];
    }

    /* put the new data on the end of the window */

    for (int i = (N_SAMP - 1), j = 0; i < (N_SAMP * 2); i++, j++) {
        Sine_Sn_[i] = sw_[j] * Parzen[i];
    }
}

