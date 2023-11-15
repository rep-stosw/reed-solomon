# reed-solomon
Optimized Reed-Solomon GF(2^12) codec for Allwinner T113-s3 SoC

Based on the codec from Phil Karn KA9Q

--------------------------------------------------------------------------------

Codec parameters when measuring:

K=3189 - number of data symbols

E=779  - number of parity symbols

1 symbol - 12 bit

--------------------------------------------------------------------------------

Unoptimized codec:

Encoding: 20.8 FPS

Decoding (maximum error correction): 11.8 FPS

Search for syndromes (decoding without errors): 19.1 FPS

Chien Search: 59.4 FPS

--------------------------------------------------------------------------------

Optimized codec:

Encoded: 37.1 FPS

Syndrome search (without Decomposition): 25.0 FPS => Decoding: 14.7 FPS

Syndrome Search(4-Step Decomposition):  258.5 FPS => Decoding: 31.6 FPS

Chien Search: 67.8 FPS

--------------------------------------------------------------------------------

GCC Toolchain: arm-none-eabi v. 10.3

Compilation flags: -Ofast -marm -mcpu=cortex-a7 -mfloat-abi=hard -mfpu=vfpv4 -mfpu=neon -ftree-vectorize -fno-math-errno -ffunction-sections -fdata-sections

--------------------------------------------------------------------------------

Bibliography:

1) https://ipnpr.jpl.nasa.gov/progress_report/42-43/43Q.PDF
2) https://ipnpr.jpl.nasa.gov/progress_report2/42-52/52J.PDF
3) https://ipnpr.jpl.nasa.gov/progress_report2/42-54/54K.PDF

--------------------------------------------------------------------------------

