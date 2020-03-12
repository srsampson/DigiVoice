#### An implementation of the Codec2 Mode 700C Speech Vocoder
This is a highly modified implementation of the Codec2 Mode 700C Speech Vocoder. It is written in C99 using complex float, the ```-O2``` optimization, and ```-fsingle-precision-constant``` enabled. It was developed and compiled using Netbeans 8.2 and an export ```ZIP``` containing the build is included. You should be able to burst this file on a Raspberry Pi and type ```make```. It also should be 32-bit compatible but I haven't tested that yet.

This speech vocoder converts a block of 320 16-bit signed PCM integers sampled at 8 kHz, into an array of four (4) unsigned 16 bit integers which contain the quantized version of the compressed speech, and vice-versa, every 40 ms, or a 25 Hz rate.

The software is provided as a dynamic library ```libdigivoice.so``` and you should add it to your ```/usr/local/lib``` directory, and also include the ```digivoice.h``` and ```c3file.h``` files in your ```/usr/local/include``` and perform a ```sudo ldconfig``` to make sure your Linux library references are updated.

Also included, are a sample ```encode``` and ```decode``` programs, and a test voice ```raw``` file.

#### Simplified Theory
The program operates on blocks of PCM audio that is sampled at 8 kHz into 15-bits plus sign. A block of 320 samples (equivalent to 40ms) is then analyzed as four segments of 80 samples (10ms).

These four segments are filtered and converted to the frequency domain, where the pitch and number of harmonics is determined. This information is put into a Model structure.

The model structure is further populated with the harmonic phases and amplitudes, and decides whether the data represents a voiced (vowel sounding) or unvoiced (consonant, quiet, or maybe noise).

The number of harmonics can vary, so they are resampled into a fixed set of 20 harmonics, and these are compared with a codebook of the same number, and the index of the values that produce the smallest error is added to the Model. This codebook was developed from samples of English speaking male and female voices. It is very very large, and will consume a good bit of ROM.

The Model data is then quantized to an array of unsigned integers and output to the user.

Decoding the array of integers into PCM time samples, uses a similar procedure. In this case, however, four Models are created. One for each 10ms sound segment are synthesized. The quantized voice index is used to get the fixed size 20 harmonics from the codebook, and then using the phase, amplitude, and voicing information, convert back to a number of harmonics compatible with the original (up to 80). This is then converted back to the time domain and issued out as a block of PCM voice samples.

#### Disclaimer
This library is being used where compatability with the original is not the primary goal. The different FFT and rounding precision of this version, produces differences which sound and look about the same.

#### Binary Array Format
```
Index 1 VQ Magnitude index 1     (9 bits)
Index 2 VQ Magnitide index 2     (9 bits)
Index 3 Quantized energy         (4 bits)
Index 4 Quantized pitch          (6 bits)

A Pitch of 0 means it is an unvoiced frame
```
#### Library Calls
```
int codec_create(void);
```
Returns 0 on success, -1 if the Sinusoidal process init fails, and -2 if the Non Linear Pitch (NLP) init fails
```
void codec_destroy(void);
```
Always call this function before ending vocoder processing, which releases allocated memory.
```
void codec_encode(uint16_t [], int16_t []);
```
Called with a buffer of 320 16-bit signed PCM samples with an 8 kHz sample rate.   
Returns a buffer of four indexed integers containing quantized speech-vocoder data.
```
void codec_decode(int16_t [], uint16_t []);
```
Called with a buffer of four indexed integers of quantized speech-vocoder data.   
Returns a buffer of 320 16-bit signed PCM samples with an 8 kHz sample rate.
```
int codec_indexes_per_frame(void);
```
Always returns the number 4.
```
int codec_samples_per_frame(void);
```
Always returns the number 320.
```
float codec_get_energy(uint16_t []);
```
Returns the mean energy of the given buffer of indexed quantized integers.
