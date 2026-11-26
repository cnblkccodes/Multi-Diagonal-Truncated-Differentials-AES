#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include <math.h>
#include <omp.h>

#ifndef NUM_THREADS
#define NUM_THREADS 12
#endif
#ifndef RNG_SEED
#define RNG_SEED ((uint64_t)time(NULL))
#endif

#ifndef NUM_KEYS
#define NUM_KEYS 100u
#endif

#ifndef NUM_PLAINTEXTS_LOG2
#define NUM_PLAINTEXTS_LOG2 20
#endif

#define NUM_PLAINTEXTS  (UINT64_C(1) << NUM_PLAINTEXTS_LOG2)
#define AES_BLOCK_BYTES 16

/* xoshiro256** PRNG, initialized with SplitMix64. */
static uint64_t rng_state[4];

static inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t next_u64(void) {
    const uint64_t result = rotl64(rng_state[1] * 5, 7) * 9;
    const uint64_t t = rng_state[1] << 17;

    rng_state[2] ^= rng_state[0];
    rng_state[3] ^= rng_state[1];
    rng_state[1] ^= rng_state[2];
    rng_state[0] ^= rng_state[3];
    rng_state[2] ^= t;
    rng_state[3] = rotl64(rng_state[3], 45);

    return result;
}

static void seed_rng(uint64_t seed) {
    for (int i = 0; i < 4; ++i) {
        seed += UINT64_C(0x9E3779B97F4A7C15);
        uint64_t z = seed;
        z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
        rng_state[i] = z ^ (z >> 31);
    }
}

static void random_bytes(uint8_t *out, size_t n) {
    size_t i = 0;
    while (i < n) {
        uint64_t x = next_u64();
        for (int j = 0; j < 8 && i < n; ++j, ++i) {
            out[i] = (uint8_t)x;
            x >>= 8;
        }
    }
}

static const uint8_t SBOX[256] = {
    0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
    0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0,
    0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15,
    0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75,
    0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84,
    0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf,
    0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8,
    0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2,
    0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73,
    0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb,
    0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79,
    0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08,
    0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a,
    0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e,
    0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf,
    0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16
};

static const uint8_t RCON[10] =
    {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,0x1B,0x36};

/* AES-128 key expansion: 176 bytes, corresponding to 11 round keys. */
static void aes128_key_expand(const uint8_t key[16], uint8_t rk[176]) {
    memcpy(rk, key, 16);

    int generated = 16;
    int rcon_index = 0;
    uint8_t t[4];

    while (generated < 176) {
        memcpy(t, rk + generated - 4, 4);

        if ((generated % 16) == 0) {
            uint8_t first = t[0];
            t[0] = SBOX[t[1]] ^ RCON[rcon_index++];
            t[1] = SBOX[t[2]];
            t[2] = SBOX[t[3]];
            t[3] = SBOX[first];
        }

        for (int i = 0; i < 4; ++i) {
            rk[generated] = rk[generated - 16] ^ t[i];
            ++generated;
        }
    }
}

static inline uint8_t xtime(uint8_t x) {
    return (uint8_t)((x << 1) ^ ((x & 0x80) ? 0x1B : 0x00));
}

static void add_round_key(uint8_t s[16], const uint8_t *rk) {
    for (int i = 0; i < 16; ++i)
        s[i] ^= rk[i];
}

static void sub_bytes(uint8_t s[16]) {
    for (int i = 0; i < 16; ++i)
        s[i] = SBOX[s[i]];
}

/* Standard AES state storage: index = row + 4*column. */
static void shift_rows(uint8_t s[16]) {
    uint8_t t[16];

    t[0]=s[0];   t[4]=s[4];   t[8]=s[8];   t[12]=s[12];
    t[1]=s[5];   t[5]=s[9];   t[9]=s[13];  t[13]=s[1];
    t[2]=s[10];  t[6]=s[14];  t[10]=s[2];  t[14]=s[6];
    t[3]=s[15];  t[7]=s[3];   t[11]=s[7];  t[15]=s[11];

    memcpy(s, t, 16);
}

static void mix_columns(uint8_t s[16]) {
    for (int c = 0; c < 4; ++c) {
        uint8_t *a = s + 4*c;
        uint8_t x0=a[0], x1=a[1], x2=a[2], x3=a[3];
        uint8_t all = x0 ^ x1 ^ x2 ^ x3;

        a[0] ^= all ^ xtime((uint8_t)(x0 ^ x1));
        a[1] ^= all ^ xtime((uint8_t)(x1 ^ x2));
        a[2] ^= all ^ xtime((uint8_t)(x2 ^ x3));
        a[3] ^= all ^ xtime((uint8_t)(x3 ^ x0));
    }
}

static void aes128_encrypt(const uint8_t in[16], uint8_t out[16],
                           const uint8_t rk[176]) {
    uint8_t s[16];
    memcpy(s, in, 16);

    add_round_key(s, rk);

    for (int round = 1; round <= 4; ++round) {
        sub_bytes(s);
        shift_rows(s);
        mix_columns(s);
        add_round_key(s, rk + 16*round);
    }

    /* The final AES round omits MixColumns. */
    sub_bytes(s);
    shift_rows(s);
    add_round_key(s, rk + 16*5);

    memcpy(out, s, 16);
}


/* Counts only EXACTLY one or EXACTLY two zero-difference bytes.
   Bit positions are AES internal column-major byte indices 0..15. */
static void classify_pair(const uint8_t *a, const uint8_t *b,
                          uint64_t one[16], uint64_t two[16][16]) {
    int passive0 = -1, passive1 = -1;
    int zeros = 0;
    for (int k = 0; k < 16; ++k) {
        if (a[k] == b[k]) {
            if (zeros == 0) passive0 = k;
            else if (zeros == 1) passive1 = k;
            if (++zeros > 2) return;
        }
    }
    if (zeros == 1) ++one[passive0];
    else if (zeros == 2) ++two[passive0][passive1];
}

static void print_pattern(uint64_t count, uint64_t n, long double p) {
    long double obs = (long double)count / (long double)n;
    long double ratio = obs / p;
    printf("2^%.15f,2^%.15f,%+.10f,%.15f,%" PRIu64 "\n",
           (double)(count ? log2l(obs) : -INFINITY), (double)log2l(p),
           (double)(100.0L*(ratio-1.0L)), (double)ratio, count);
}

static void print_sum(const char *label, uint64_t count, uint64_t n, long double p) {
    long double obs = (long double)count / (long double)n;
    long double ratio = obs / p;
    printf("%s: 2^%.15f (theoretical 2^%.15f)\n", label,
           (double)(count ? log2l(obs) : -INFINITY), (double)log2l(p));
    printf("Ratio: %.15f | Relative deviation: %+.10f%% | Count: %" PRIu64 "\n",
           (double)ratio, (double)(100.0L*(ratio-1.0L)), count);
}

static void print_results(const char *heading,
                          const uint64_t one[16],
                          const uint64_t two[16][16],
                          uint64_t num_pairs) {
    const long double q = 1.0L / 256.0L;
    const long double correction = 1.0L - ldexpl(1.0L, -128);
    const long double p1 = q * powl(1.0L-q, 15) / correction;
    const long double p2 = q*q * powl(1.0L-q, 14) / correction;
    uint64_t sum1 = 0, sum2 = 0;
    printf("\n%s\n", heading);
    printf("Total tested pairs: %" PRIu64 "\n", num_pairs);
    printf("\nEXACTLY 1 passive, 15 active bytes (16 patterns):\n");
    printf("passive_byte,observed_probability_2pow,theoretical_probability_2pow,relative_deviation_percent,observed_over_theoretical,count\n");
    for (int i = 0; i < 16; ++i) {
        printf("%d,", i);
        print_pattern(one[i], num_pairs, p1);
        sum1 += one[i];
    }
    printf("\nEXACTLY 2 passive, 14 active bytes (120 patterns):\n");
    printf("passive_byte_1,passive_byte_2,observed_probability_2pow,theoretical_probability_2pow,relative_deviation_percent,observed_over_theoretical,count\n");
    for (int i = 0; i < 16; ++i) {
        for (int j = i+1; j < 16; ++j) {
            printf("%d,%d,", i, j);
            print_pattern(two[i][j], num_pairs, p2);
            sum2 += two[i][j];
        }
    }
    printf("\n");
    print_sum("SUM exactly 1 passive probability", sum1, num_pairs, 16.0L*p1);
    print_sum("SUM exactly 2 passive probability", sum2, num_pairs, 120.0L*p2);
}

int main(void) {
    seed_rng(RNG_SEED);
    omp_set_dynamic(0);
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        #pragma omp single
        printf("OpenMP threads: %d\n", omp_get_num_threads());
    }
    if (NUM_PLAINTEXTS < 3 || NUM_PLAINTEXTS > SIZE_MAX / AES_BLOCK_BYTES) {
        fprintf(stderr, "Invalid plaintext count / buffer size.\n");
        return EXIT_FAILURE;
    }

    /* Exclude first OFB input X_1 = IV from pair comparisons. */
    const uint64_t kept = NUM_PLAINTEXTS - 1;
    const uint64_t pairs_per_key = kept * (kept - 1) / 2;
    const size_t buffer_bytes = (size_t)kept * AES_BLOCK_BYTES;
    uint8_t *inputs = malloc(buffer_bytes);
    if (!inputs) { perror("malloc"); return EXIT_FAILURE; }

    uint64_t total_one[16] = {0};
    uint64_t total_two[16][16] = {{0}};
    printf("AES-128-OFB128: keys=%u, blocks/key=2^%d, retained inputs/key=%" PRIu64 "\n",
           NUM_KEYS, NUM_PLAINTEXTS_LOG2, kept);
    printf("Pairs/key=%" PRIu64 ", total pairs=%" PRIu64 "\n",
           pairs_per_key, pairs_per_key * (uint64_t)NUM_KEYS);
    fflush(stdout);

    for (unsigned key_index = 0; key_index < NUM_KEYS; ++key_index) {
        uint8_t key[16], rk[176], prev[16];
        random_bytes(key, 16);
        random_bytes(prev, 16); /* IV */
        aes128_key_expand(key, rk);

        printf("KEY %u/%u | Generating OFB inputs...\n", key_index+1, NUM_KEYS);
        fflush(stdout);
        for (uint64_t i = 0; i < NUM_PLAINTEXTS; ++i) {
            uint8_t next[16];
            /* OFB-128: X_i = O_{i-1}; O_i = AES_K(O_{i-1}), O_0 = IV.
               Plaintext/ciphertext do not affect the feedback state. */
            if (i > 0) memcpy(inputs + (size_t)(i-1)*16, prev, 16);
            aes128_encrypt(prev, next, rk);
            memcpy(prev, next, 16);
        }

        uint64_t one[16] = {0}, two[16][16] = {{0}};
        printf("KEY %u/%u | Comparing pairs...\n", key_index+1, NUM_KEYS);
        fflush(stdout);
        const double start = omp_get_wtime();
        uint64_t completed = 0;
        unsigned milestone = 1;
        const uint64_t chunks = (kept - 1 + 63) / 64;
        #pragma omp parallel
        {
            uint64_t local_one[16] = {0}, local_two[16][16] = {{0}};
            #pragma omp for schedule(dynamic, 1)
            for (uint64_t chunk = 0; chunk < chunks; ++chunk) {
                uint64_t begin = chunk * 64;
                uint64_t end = begin + 64;
                if (end > kept - 1) end = kept - 1;
                uint64_t done = 0;
                for (uint64_t i = begin; i < end; ++i) {
                    const uint8_t *a = inputs + (size_t)i*16;
                    for (uint64_t j = i+1; j < kept; ++j) {
                        const uint8_t *b = inputs + (size_t)j*16;
                        classify_pair(a, b, local_one, local_two);
                    }
                    done += kept - i - 1;
                }
                #pragma omp critical(progress)
                {
                    completed += done;
                    if (milestone <= 4 &&
                        (long double)completed / pairs_per_key >= milestone / 4.0L) {
                        double elapsed = omp_get_wtime() - start;
                        double fraction = (double)completed / pairs_per_key;
                        printf("KEY %u/%u | Progress: %.2f%% | Pairs: %" PRIu64 "/%" PRIu64
                               " | Elapsed: %.1f s | ETA: %.1f s\n",
                               key_index+1, NUM_KEYS, 100.0*fraction,
                               completed, pairs_per_key, elapsed,
                               elapsed*(1.0/fraction-1.0));
                        fflush(stdout);
                        while (milestone <= 4 &&
                               (long double)completed / pairs_per_key >= milestone / 4.0L)
                            ++milestone;
                    }
                }
            }
            #pragma omp critical(merge_counts)
            {
                for (int i = 0; i < 16; ++i) {
                    one[i] += local_one[i];
                    for (int j = i+1; j < 16; ++j)
                        two[i][j] += local_two[i][j];
                }
            }
        }
        if (completed != pairs_per_key) {
            fprintf(stderr, "Pair count mismatch!\n");
            free(inputs);
            return EXIT_FAILURE;
        }
        for (int i = 0; i < 16; ++i) {
            total_one[i] += one[i];
            for (int j = i+1; j < 16; ++j) total_two[i][j] += two[i][j];
        }
        char label[80];
        snprintf(label, sizeof label, "KEY %u", key_index+1);
        print_results(label, one, two, pairs_per_key);
        fflush(stdout);
    }
    print_results("TOTAL OVER ALL KEYS", total_one, total_two,
                  pairs_per_key * (uint64_t)NUM_KEYS);
    free(inputs);
    return EXIT_SUCCESS;
}
