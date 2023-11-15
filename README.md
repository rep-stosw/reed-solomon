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

Encoding: 23.1 FPS

Decoding (maximum error correction): 15.2 FPS

Search for syndromes (decoding without errors): 25.2 FPS

--------------------------------------------------------------------------------

Optimized codec:

Encoded: 45.9 FPS

Syndrome search (without Decomposition): 26.7 FPS => Decoding: 15.8 FPS

Syndrome Search(4-Step Decomposition):  285.5 FPS => Decoding: 34.4 FPS

--------------------------------------------------------------------------------

GCC Toolchain: arm-none-eabi v. 10.3

Compilation flags: -Ofast -mthumb -mcpu=cortex-a7 -mfloat-abi=hard -mfpu=vfpv4 -mfpu=neon -ftree-vectorize -fno-math-errno -ffunction-sections -fdata-sections

--------------------------------------------------------------------------------

Bibliography:

1) https://ipnpr.jpl.nasa.gov/progress_report/42-43/43Q.PDF
2) https://ipnpr.jpl.nasa.gov/progress_report2/42-52/52J.PDF
3) https://ipnpr.jpl.nasa.gov/progress_report2/42-54/54K.PDF

--------------------------------------------------------------------------------

