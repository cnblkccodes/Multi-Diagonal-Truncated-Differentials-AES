#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <inttypes.h>

#ifndef NUM_KEYS
#define NUM_KEYS 100
#endif
#ifndef NUM_THREADS
#define NUM_THREADS 12
#endif
#ifndef RNG_SEED
#define RNG_SEED ((uint32_t)time(NULL))
#endif

#define STATE_NIBBLES        16
#define BLOCK_NIBBLES        16
#define KEY_BYTES            8
#define ROUNDS               5

#ifndef STRUCTURES_LOG2
#define STRUCTURES_LOG2 12
#endif
#define NUM_STRUCTURES (UINT64_C(1) << STRUCTURES_LOG2)
#define PAIRS_PER_STRUCTURE UINT64_C(1658880000)
#define TEXTS_PER_STRUCTURE  65536
#define NUM_COUNTERS         14
#ifndef BRUTE_FORCE
#define BRUTE_FORCE 0
#endif
#if NUM_KEYS < 1 || NUM_THREADS < 1 || STRUCTURES_LOG2 < 2 || STRUCTURES_LOG2 > 48
#error Invalid experiment parameters
#endif

static uint32_t rng_state[4];

static inline uint32_t rotl(const uint32_t x, int k)
{
    return (x << k) | (x >> (32 - k));
}

static uint32_t next_prng(void)
{
    const uint32_t result =
        rotl(rng_state[0] + rng_state[3], 7) + rng_state[0];

    const uint32_t t = rng_state[1] << 9;

    rng_state[2] ^= rng_state[0];
    rng_state[3] ^= rng_state[1];
    rng_state[1] ^= rng_state[2];
    rng_state[0] ^= rng_state[3];

    rng_state[2] ^= t;
    rng_state[3] = rotl(rng_state[3], 11);

    return result;
}

static void init_prng(uint32_t seed)
{
    rng_state[0] = seed ^ 0x12345678;
    rng_state[1] = seed ^ 0x87654321;
    rng_state[2] = seed ^ 0xDEADBEEF;
    rng_state[3] = seed ^ 0xCAFEBABE;

    for (int i = 0; i < 20; ++i)
        next_prng();
}

static const uint8_t SMALL_AES_SBOX[16] = {
    0x6, 0xB, 0x5, 0x4,
    0x2, 0xE, 0x7, 0xA,
    0x9, 0xD, 0xF, 0xC,
    0x3, 0x1, 0x0, 0x8
};

static uint8_t gf16_mul(uint8_t a, uint8_t b)
{
    uint8_t result = 0;
    a &= 0xF;
    b &= 0xF;

    for (int i = 0; i < 4; ++i)
    {
        if (b & 1)
            result ^= a;

        b >>= 1;
        uint8_t carry = a & 0x8;
        a <<= 1;

        if (carry)
            a ^= 0x13;

        a &= 0xF;
    }

    return result & 0xF;
}

static void sub_bytes(uint8_t state[STATE_NIBBLES])
{
    for (int i = 0; i < STATE_NIBBLES; ++i)
        state[i] = SMALL_AES_SBOX[state[i] & 0xF];
}

static void shift_rows(uint8_t state[STATE_NIBBLES])
{
    uint8_t tmp[STATE_NIBBLES];

    tmp[0]  = state[0];
    tmp[1]  = state[5];
    tmp[2]  = state[10];
    tmp[3]  = state[15];

    tmp[4]  = state[4];
    tmp[5]  = state[9];
    tmp[6]  = state[14];
    tmp[7]  = state[3];

    tmp[8]  = state[8];
    tmp[9]  = state[13];
    tmp[10] = state[2];
    tmp[11] = state[7];

    tmp[12] = state[12];
    tmp[13] = state[1];
    tmp[14] = state[6];
    tmp[15] = state[11];

    memcpy(state, tmp, STATE_NIBBLES);
}

static void mix_columns(uint8_t state[STATE_NIBBLES])
{
    uint8_t tmp[STATE_NIBBLES];

    for (int col = 0; col < 4; ++col)
    {
        int base = 4 * col;

        uint8_t c0 = state[base + 0] & 0xF;
        uint8_t c1 = state[base + 1] & 0xF;
        uint8_t c2 = state[base + 2] & 0xF;
        uint8_t c3 = state[base + 3] & 0xF;

        tmp[base + 0] =
            gf16_mul(0x2, c0) ^ gf16_mul(0x3, c1) ^ c2 ^ c3;

        tmp[base + 1] =
            c0 ^ gf16_mul(0x2, c1) ^ gf16_mul(0x3, c2) ^ c3;

        tmp[base + 2] =
            c0 ^ c1 ^ gf16_mul(0x2, c2) ^ gf16_mul(0x3, c3);

        tmp[base + 3] =
            gf16_mul(0x3, c0) ^ c1 ^ c2 ^ gf16_mul(0x2, c3);

        tmp[base + 0] &= 0xF;
        tmp[base + 1] &= 0xF;
        tmp[base + 2] &= 0xF;
        tmp[base + 3] &= 0xF;
    }

    memcpy(state, tmp, STATE_NIBBLES);
}

static void add_round_key(
    uint8_t state[STATE_NIBBLES],
    const uint8_t RK[STATE_NIBBLES])
{
    for (int i = 0; i < STATE_NIBBLES; ++i)
        state[i] ^= RK[i] & 0xF;
}

static void small_aes_key_schedule(
    const uint8_t key[KEY_BYTES],
    uint8_t round_keys[ROUNDS + 1][STATE_NIBBLES])
{
    for (int i = 0; i < KEY_BYTES; ++i)
    {
        round_keys[0][2 * i] = (key[i] >> 4) & 0xF;
        round_keys[0][2 * i + 1] = key[i] & 0xF;
    }

    uint8_t kappa = 0x1;

    for (int r = 1; r <= ROUNDS; ++r)
    {
        const uint8_t *prev = round_keys[r - 1];
        uint8_t *curr = round_keys[r];
        uint8_t temp[4];

        temp[0] = SMALL_AES_SBOX[prev[13]];
        temp[1] = SMALL_AES_SBOX[prev[14]];
        temp[2] = SMALL_AES_SBOX[prev[15]];
        temp[3] = SMALL_AES_SBOX[prev[12]];

        temp[0] ^= kappa;

        for (int j = 0; j < 4; ++j)
            curr[j] = (prev[j] ^ temp[j]) & 0xF;

        for (int q = 1; q < 4; ++q)
        {
            for (int j = 0; j < 4; ++j)
            {
                curr[4 * q + j] =
                    (prev[4 * q + j] ^ curr[4 * (q - 1) + j]) & 0xF;
            }
        }

        kappa = gf16_mul(0x2, kappa);
    }
}

static void small_aes_round(
    uint8_t state[STATE_NIBBLES],
    const uint8_t RK[STATE_NIBBLES])
{
    sub_bytes(state);
    shift_rows(state);
    mix_columns(state);
    add_round_key(state, RK);
}

static void small_aes_final_round(
    uint8_t state[STATE_NIBBLES],
    const uint8_t RK[STATE_NIBBLES])
{
    sub_bytes(state);
    shift_rows(state);
    add_round_key(state, RK);
}

static void small_aes_encrypt(
    uint8_t state[STATE_NIBBLES],
    const uint8_t round_keys[ROUNDS + 1][STATE_NIBBLES])
{
    add_round_key(state, round_keys[0]);

    for (int r = 1; r < ROUNDS; ++r)
        small_aes_round(state, round_keys[r]);

    small_aes_final_round(state, round_keys[ROUNDS]);
}

static void generate_fixed_part(uint8_t fixed[16], uint64_t index, uint64_t mul, uint64_t offset)
{
    const uint64_t mask = (UINT64_C(1)<<48)-1;
    uint64_t x=(index+offset)&mask;
    x ^= x>>27; x=(x*mul)&mask;
    x ^= x>>23; x=(x*UINT64_C(0x9e3779b97f4a7c15))&mask;
    x ^= x>>29;
    for(int i=0;i<16;++i) {
        if(i==0 || i==5 || i==10 || i==15) fixed[i]=0;
        else { fixed[i]=(uint8_t)(x&15); x>>=4; }
    }
}

static void generate_plaintext(
    uint8_t plaintext[BLOCK_NIBBLES],
    const uint8_t fixed[BLOCK_NIBBLES],
    uint8_t value0,
    uint8_t value5,
    uint8_t value10, uint8_t value15)
{
    memcpy(plaintext, fixed, BLOCK_NIBBLES);
    plaintext[0] = value0 & 0xF;
    plaintext[5] = value5 & 0xF;
    plaintext[10] = value10 & 0xF;
    plaintext[15] = value15 & 0xF;
}

static inline uint8_t get_passive_inverse_diagonal_mask(
    const uint8_t c1[BLOCK_NIBBLES],
    const uint8_t c2[BLOCK_NIBBLES])
{
    uint8_t mask = 0;

    if (((c1[0]  ^ c2[0])  == 0) &&
        ((c1[7]  ^ c2[7])  == 0) &&
        ((c1[10] ^ c2[10]) == 0) &&
        ((c1[13] ^ c2[13]) == 0))
        mask |= 0x1;

    if (((c1[1]  ^ c2[1])  == 0) &&
        ((c1[4]  ^ c2[4])  == 0) &&
        ((c1[11] ^ c2[11]) == 0) &&
        ((c1[14] ^ c2[14]) == 0))
        mask |= 0x2;

    if (((c1[2]  ^ c2[2])  == 0) &&
        ((c1[5]  ^ c2[5])  == 0) &&
        ((c1[8]  ^ c2[8])  == 0) &&
        ((c1[15] ^ c2[15]) == 0))
        mask |= 0x4;

    if (((c1[3]  ^ c2[3])  == 0) &&
        ((c1[6]  ^ c2[6])  == 0) &&
        ((c1[9]  ^ c2[9])  == 0) &&
        ((c1[12] ^ c2[12]) == 0))
        mask |= 0x8;

    return mask;
}

static inline int popcount4(uint8_t x)
{
    return ((x >> 0) & 1) + ((x >> 1) & 1) +
           ((x >> 2) & 1) + ((x >> 3) & 1);
}

enum
{
    C_D0 = 0, C_D1, C_D2, C_D3,
    C_D0_D1, C_D0_D2, C_D0_D3, C_D1_D2, C_D1_D3, C_D2_D3,
    C_EXACTLY_1_PASSIVE, C_EXACTLY_2_PASSIVE,
    C_AT_LEAST_1_PASSIVE, C_AT_LEAST_2_PASSIVE
};

static const char *COUNTER_NAMES[NUM_COUNTERS] =
{
    "D0 passive; D1,D2,D3 free",
    "D1 passive; D0,D2,D3 free",
    "D2 passive; D0,D1,D3 free",
    "D3 passive; D0,D1,D2 free",
    "D0,D1 passive; D2,D3 free",
    "D0,D2 passive; D1,D3 free",
    "D0,D3 passive; D1,D2 free",
    "D1,D2 passive; D0,D3 free",
    "D1,D3 passive; D0,D2 free",
    "D2,D3 passive; D0,D1 free",
    "Exactly 1 passive + 3 active",
    "Exactly 2 passive + 2 active",
    "At least 1 passive",
    "At least 2 passive"
};

static inline void update_counters(
    uint64_t counters[NUM_COUNTERS],
    uint8_t mask)
{
    int number_of_passive = popcount4(mask);

    if (mask & 0x1) ++counters[C_D0];
    if (mask & 0x2) ++counters[C_D1];
    if (mask & 0x4) ++counters[C_D2];
    if (mask & 0x8) ++counters[C_D3];

    if ((mask & 0x3) == 0x3) ++counters[C_D0_D1];
    if ((mask & 0x5) == 0x5) ++counters[C_D0_D2];
    if ((mask & 0x9) == 0x9) ++counters[C_D0_D3];
    if ((mask & 0x6) == 0x6) ++counters[C_D1_D2];
    if ((mask & 0xA) == 0xA) ++counters[C_D1_D3];
    if ((mask & 0xC) == 0xC) ++counters[C_D2_D3];

    if (number_of_passive == 1)
        ++counters[C_EXACTLY_1_PASSIVE];

    if (number_of_passive == 2)
        ++counters[C_EXACTLY_2_PASSIVE];

    if (number_of_passive >= 1)
        ++counters[C_AT_LEAST_1_PASSIVE];

    if (number_of_passive >= 2)
        ++counters[C_AT_LEAST_2_PASSIVE];
}

static int valid_pair(int a,int b)
{
    unsigned d=(unsigned)(a^b);
    return (d&0xf000) && (d&0x0f00) && (d&0x00f0) && (d&0x000f);
}
static int self_test(void)
{
    uint8_t fixed[16], p[16], q[16];
    uint64_t hist[5]={0};
    generate_fixed_part(fixed,0,3,7);
    generate_plaintext(p,fixed,0,0,0,0);
    for(int a=0;a<65536;++a) {
        generate_plaintext(q,fixed,a>>12,(a>>8)&15,(a>>4)&15,a&15);
        unsigned mask=0;
        for(int i=0;i<16;++i) if(p[i]!=q[i]) mask|=1u<<i;
        unsigned expected=0;
        if(a&0xf000) expected|=1u;
        if(a&0x0f00) expected|=1u<<5;
        if(a&0x00f0) expected|=1u<<10;
        if(a&0x000f) expected|=1u<<15;
        if(mask!=expected) return 0;
        int n=!!(a&0xf000)+!!(a&0xf00)+!!(a&0xf0)+!!(a&15);
        if(valid_pair(0,a)!=(n==4)) return 0;
        ++hist[n];
    }
    return hist[0]==1 && hist[1]==60 && hist[2]==1350 && hist[3]==13500 && hist[4]==50625;
}

static void count_structure(const uint8_t ct[65536][16],uint64_t c[14])
{
    static const unsigned pos[4][4]={{0,7,10,13},{1,4,11,14},{2,5,8,15},{3,6,9,12}};
    int *head=malloc(65536*sizeof *head), *next=malloc(65536*sizeof *next);
    if(!head || !next) abort();
    for(unsigned d=0;d<4;++d) {
        for(int i=0;i<65536;++i) head[i]=-1;
        for(int a=0;a<65536;++a) {
            unsigned v=0;
            for(unsigned j=0;j<4;++j) v=(v<<4)|ct[a][pos[d][j]];
            for(int b=head[v];b>=0;b=next[b]) {
                if(!valid_pair(a,b)) continue;
                unsigned mask=get_passive_inverse_diagonal_mask(ct[a],ct[b]);
                if(mask==15 || !(mask&(1u<<d))) abort();
                if(mask&((1u<<d)-1)) continue;
                update_counters(c,(uint8_t)mask);
            }
            next[a]=head[v]; head[v]=a;
        }
    }
    free(head); free(next);
}
static void report(const char *label,const uint64_t c[NUM_COUNTERS],uint64_t n) {
    printf("\n%s | Pairs=%" PRIu64 "\n",label,n);
    for(int i=0;i<17;++i) {
        const char *name;
        uint64_t count;
        if(i==0) {name="0 specified passive; all free";count=n;}
        else if(i==11) {name="Exactly 0 passive + 4 active";count=n-c[C_AT_LEAST_1_PASSIVE];}
        else if(i==14) {name="At least 0 passive";count=n;}
        else {int j=i-(i>0)-(i>11)-(i>14);name=COUNTER_NAMES[j];count=c[j];}
        double p=(double)count/(double)n;
        printf("%2d. %-38s count=%" PRIu64 " p=%.17g log2(p)=%.15f\n",i+1,name,count,p,count?log2(p):-INFINITY);
    }
    fflush(stdout);
}
static void brute(const uint8_t ct[65536][16], uint64_t c[14]) {
    for(int a=0;a<65536;++a) for(int b=a+1;b<65536;++b)
        if(valid_pair(a,b)) update_counters(c,get_passive_inverse_diagonal_mask(ct[a],ct[b]));
}
static int validate_counting(void) {
    uint8_t key[8]={0}, rk[ROUNDS+1][16], fixed[16];
    uint8_t (*ct)[16]=malloc(65536*sizeof *ct);
    if(!ct) return 0;
    uint64_t a[14]={0}, b[14]={0};
    small_aes_key_schedule(key,rk);
    generate_fixed_part(fixed,0,3,7);
    for(int v=0;v<65536;++v) {
        generate_plaintext(ct[v],fixed,v>>12,(v>>8)&15,(v>>4)&15,v&15);
        small_aes_encrypt(ct[v],rk);
    }
    count_structure(ct,a); brute(ct,b); free(ct);
    return memcmp(a,b,sizeof a)==0;
}
int main(void)
{
    if(!self_test()) {fprintf(stderr,"Structure self-test failed\n");return 1;}
    printf("Structure self-test passed: 50625 exact m=4 differences per plaintext; 1658880000 pairs/structure\n");
    if(NUM_STRUCTURES > UINT64_MAX/PAIRS_PER_STRUCTURE/(uint64_t)NUM_KEYS) {fprintf(stderr,"Pair count overflow\n");return 1;}
    printf("Checking grouping against brute-force (one full structure)...\n"); fflush(stdout);
    if(!validate_counting()) {fprintf(stderr,"Grouping self-test failed\n");return 1;}
    printf("Grouping self-test passed: one full structure vs brute-force\n");
    printf("Counting method: %s\n",BRUTE_FORCE?"brute-force":"grouped"); fflush(stdout);
    uint32_t seed=RNG_SEED;
    init_prng(seed);
    omp_set_dynamic(0); omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel
    {
        #pragma omp single
        printf("OpenMP threads: %d | seed=%" PRIu32 "\n",omp_get_num_threads(),seed);
    }
    const uint64_t pairs_per_key=NUM_STRUCTURES*PAIRS_PER_STRUCTURE;
    uint64_t total[14]={0}, total_pairs=0;
    printf("Small-AES m=4 | rounds=%d keys=%d structures/key=%" PRIu64 " pairs/key=%" PRIu64 "\n",ROUNDS,NUM_KEYS,NUM_STRUCTURES,pairs_per_key);
    for(int k=0;k<NUM_KEYS;++k) {
        uint8_t key[8],rk[ROUNDS+1][16];
        for(int i=0;i<8;++i) key[i]=(uint8_t)next_prng();
        small_aes_key_schedule(key,rk);
        uint64_t mul=(uint64_t)next_prng()<<32; mul|=next_prng(); mul|=1;
        uint64_t offset=(uint64_t)next_prng()<<32; offset|=next_prng();
        uint64_t counts[14]={0}, tested=0;
        double start=omp_get_wtime();
        printf("KEY %d/%d | master key=",k+1,NUM_KEYS);
        for(int i=0;i<8;++i) printf("%02x",(unsigned)key[i]);
        printf("\n"); fflush(stdout);
        for(unsigned stage=0;stage<4;++stage) {
            uint64_t begin=NUM_STRUCTURES*stage/4,end=NUM_STRUCTURES*(stage+1)/4;
            #pragma omp parallel
            {
                uint64_t local[14]={0}, local_pairs=0;
                #pragma omp for schedule(static)
                for(uint64_t j=begin;j<end;++j) {
                    uint8_t fixed[16];
                    uint8_t (*ct)[16]=malloc(65536*sizeof *ct);
                    if(!ct) abort();
                    generate_fixed_part(fixed,j,mul,offset);
                    for(int a=0;a<65536;++a) {
                        generate_plaintext(ct[a],fixed,a>>12,(a>>8)&15,(a>>4)&15,a&15);
                        small_aes_encrypt(ct[a],rk);
                    }
                    uint64_t one[14]={0};
                    if(BRUTE_FORCE) brute(ct,one); else count_structure(ct,one);
#ifdef VERIFY_BRUTE
                    uint64_t reference[14]={0};
                    for(int a=0;a<65536;++a) for(int b=a+1;b<65536;++b)
                        if(valid_pair(a,b)) update_counters(reference,get_passive_inverse_diagonal_mask(ct[a],ct[b]));
                    for(int i=0;i<14;++i) if(reference[i]!=one[i]) abort();
#endif
                    for(int i=0;i<14;++i) local[i]+=one[i];
                    free(ct);
                    local_pairs+=PAIRS_PER_STRUCTURE;
                }
                #pragma omp critical
                {
                    for(int i=0;i<14;++i) counts[i]+=local[i];
                    tested+=local_pairs;
                }
            }
            if(tested!=end*PAIRS_PER_STRUCTURE) {fprintf(stderr,"Pair count mismatch\n");return 1;}
            double elapsed=omp_get_wtime()-start;
            printf("KEY %d/%d | Progress: %u%% | Structures=%" PRIu64 "/%" PRIu64 " | Pairs=%" PRIu64 " | Elapsed=%.1f s | ETA=%.1f s\n",k+1,NUM_KEYS,25*(stage+1),end,NUM_STRUCTURES,tested,elapsed,elapsed*(3-stage)/(stage+1));
            fflush(stdout);
        }
        if(counts[C_AT_LEAST_1_PASSIVE]>tested || tested!=pairs_per_key || counts[C_AT_LEAST_1_PASSIVE]!=counts[C_EXACTLY_1_PASSIVE]+counts[C_AT_LEAST_2_PASSIVE] || counts[C_AT_LEAST_2_PASSIVE]<counts[C_EXACTLY_2_PASSIVE]) {fprintf(stderr,"Counter consistency failure\n");return 1;}
        for(int i=0;i<14;++i) total[i]+=counts[i];
        total_pairs+=tested;
        char label[64]; snprintf(label,sizeof label,"KEY %d FINISHED",k+1);
        report(label,counts,tested);
    }
    report("ALL KEYS FINISHED",total,total_pairs);
    return 0;
}
