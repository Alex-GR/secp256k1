/**********************************************************************
 * Copyright (c) 2013, 2014 Pieter Wuille                             *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _SECP256K1_SCALAR_REPR_IMPL_H_
#define _SECP256K1_SCALAR_REPR_IMPL_H_

/* Limbs of the secp256k1 order. */
#define SECP256K1_N_0 ((uint64_t)0xBFD25E8CD0364141ULL)
#define SECP256K1_N_1 ((uint64_t)0xBAAEDCE6AF48A03BULL)
#define SECP256K1_N_2 ((uint64_t)0xFFFFFFFFFFFFFFFEULL)
#define SECP256K1_N_3 ((uint64_t)0xFFFFFFFFFFFFFFFFULL)

/* Limbs of 2^256 minus the secp256k1 order. */
#define SECP256K1_N_C_0 (~SECP256K1_N_0 + 1)
#define SECP256K1_N_C_1 (~SECP256K1_N_1)
#define SECP256K1_N_C_2 (1)

/* Limbs of half the secp256k1 order. */
#define SECP256K1_N_H_0 ((uint64_t)0xDFE92F46681B20A0ULL)
#define SECP256K1_N_H_1 ((uint64_t)0x5D576E7357A4501DULL)
#define SECP256K1_N_H_2 ((uint64_t)0xFFFFFFFFFFFFFFFFULL)
#define SECP256K1_N_H_3 ((uint64_t)0x7FFFFFFFFFFFFFFFULL)

SECP256K1_INLINE static void secp256k1_scalar_clear(secp256k1_scalar *r) {
    r->d[0] = 0;
    r->d[1] = 0;
    r->d[2] = 0;
    r->d[3] = 0;
}

SECP256K1_INLINE static void secp256k1_scalar_set_int(secp256k1_scalar *r, unsigned int v) {
    r->d[0] = v;
    r->d[1] = 0;
    r->d[2] = 0;
    r->d[3] = 0;
}

SECP256K1_INLINE static unsigned int secp256k1_scalar_get_bits(const secp256k1_scalar *a, unsigned int offset, unsigned int count) {
    VERIFY_CHECK((offset + count - 1) >> 6 == offset >> 6);
    return (a->d[offset >> 6] >> (offset & 0x3F)) & ((((uint64_t)1) << count) - 1);
}

SECP256K1_INLINE static unsigned int secp256k1_scalar_get_bits_var(const secp256k1_scalar *a, unsigned int offset, unsigned int count) {
    VERIFY_CHECK(count < 32);
    VERIFY_CHECK(offset + count <= 256);
    if ((offset + count - 1) >> 6 == offset >> 6) {
        return secp256k1_scalar_get_bits(a, offset, count);
    } else {
        VERIFY_CHECK((offset >> 6) + 1 < 4);
        return ((a->d[offset >> 6] >> (offset & 0x3F)) | (a->d[(offset >> 6) + 1] << (64 - (offset & 0x3F)))) & ((((uint64_t)1) << count) - 1);
    }
}

SECP256K1_INLINE static int secp256k1_scalar_check_overflow(const secp256k1_scalar *a) {
    int yes = 0;
    int no = 0;
    no |= (a->d[3] < SECP256K1_N_3); /* No need for a > check. */
    no |= (a->d[2] < SECP256K1_N_2);
    yes |= (a->d[2] > SECP256K1_N_2) & ~no;
    no |= (a->d[1] < SECP256K1_N_1);
    yes |= (a->d[1] > SECP256K1_N_1) & ~no;
    yes |= (a->d[0] >= SECP256K1_N_0) & ~no;
    return yes;
}

SECP256K1_INLINE static int secp256k1_scalar_reduce(secp256k1_scalar *r, unsigned int overflow) {
    uint128_t t;
    VERIFY_CHECK(overflow <= 1);
    t = (uint128_t)r->d[0] + overflow * SECP256K1_N_C_0;
    r->d[0] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
    t += (uint128_t)r->d[1] + overflow * SECP256K1_N_C_1;
    r->d[1] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
    t += (uint128_t)r->d[2] + overflow * SECP256K1_N_C_2;
    r->d[2] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
    t += (uint64_t)r->d[3];
    r->d[3] = t & 0xFFFFFFFFFFFFFFFFULL;
    return overflow;
}

static int secp256k1_scalar_add(secp256k1_scalar *r, const secp256k1_scalar *a, const secp256k1_scalar *b) {
    int overflow;
    uint128_t t = (uint128_t)a->d[0] + b->d[0];
    r->d[0] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
    t += (uint128_t)a->d[1] + b->d[1];
    r->d[1] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
    t += (uint128_t)a->d[2] + b->d[2];
    r->d[2] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
    t += (uint128_t)a->d[3] + b->d[3];
    r->d[3] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
    overflow = t + secp256k1_scalar_check_overflow(r);
    VERIFY_CHECK(overflow == 0 || overflow == 1);
    secp256k1_scalar_reduce(r, overflow);
    return overflow;
}

static void secp256k1_scalar_cadd_bit(secp256k1_scalar *r, unsigned int bit, int flag) {
    uint128_t t;
    VERIFY_CHECK(bit < 256);
    bit += ((uint32_t) flag - 1) & 0x100;  /* forcing (bit >> 6) > 3 makes this a noop */
    t = (uint128_t)r->d[0] + (((uint64_t)((bit >> 6) == 0)) << (bit & 0x3F));
    r->d[0] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
    t += (uint128_t)r->d[1] + (((uint64_t)((bit >> 6) == 1)) << (bit & 0x3F));
    r->d[1] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
    t += (uint128_t)r->d[2] + (((uint64_t)((bit >> 6) == 2)) << (bit & 0x3F));
    r->d[2] = t & 0xFFFFFFFFFFFFFFFFULL; t >>= 64;
    t += (uint128_t)r->d[3] + (((uint64_t)((bit >> 6) == 3)) << (bit & 0x3F));
    r->d[3] = t & 0xFFFFFFFFFFFFFFFFULL;
#ifdef VERIFY
    VERIFY_CHECK((t >> 64) == 0);
    VERIFY_CHECK(secp256k1_scalar_check_overflow(r) == 0);
#endif
}

static void secp256k1_scalar_set_b32(secp256k1_scalar *r, const unsigned char *b32, int *overflow) {
    int over;
    r->d[0] = (uint64_t)b32[31] | (uint64_t)b32[30] << 8 | (uint64_t)b32[29] << 16 | (uint64_t)b32[28] << 24 | (uint64_t)b32[27] << 32 | (uint64_t)b32[26] << 40 | (uint64_t)b32[25] << 48 | (uint64_t)b32[24] << 56;
    r->d[1] = (uint64_t)b32[23] | (uint64_t)b32[22] << 8 | (uint64_t)b32[21] << 16 | (uint64_t)b32[20] << 24 | (uint64_t)b32[19] << 32 | (uint64_t)b32[18] << 40 | (uint64_t)b32[17] << 48 | (uint64_t)b32[16] << 56;
    r->d[2] = (uint64_t)b32[15] | (uint64_t)b32[14] << 8 | (uint64_t)b32[13] << 16 | (uint64_t)b32[12] << 24 | (uint64_t)b32[11] << 32 | (uint64_t)b32[10] << 40 | (uint64_t)b32[9] << 48 | (uint64_t)b32[8] << 56;
    r->d[3] = (uint64_t)b32[7] | (uint64_t)b32[6] << 8 | (uint64_t)b32[5] << 16 | (uint64_t)b32[4] << 24 | (uint64_t)b32[3] << 32 | (uint64_t)b32[2] << 40 | (uint64_t)b32[1] << 48 | (uint64_t)b32[0] << 56;
    over = secp256k1_scalar_reduce(r, secp256k1_scalar_check_overflow(r));
    if (overflow) {
        *overflow = over;
    }
}

static void secp256k1_scalar_get_b32(unsigned char *bin, const secp256k1_scalar* a) {
    bin[0] = a->d[3] >> 56; bin[1] = a->d[3] >> 48; bin[2] = a->d[3] >> 40; bin[3] = a->d[3] >> 32; bin[4] = a->d[3] >> 24; bin[5] = a->d[3] >> 16; bin[6] = a->d[3] >> 8; bin[7] = a->d[3];
    bin[8] = a->d[2] >> 56; bin[9] = a->d[2] >> 48; bin[10] = a->d[2] >> 40; bin[11] = a->d[2] >> 32; bin[12] = a->d[2] >> 24; bin[13] = a->d[2] >> 16; bin[14] = a->d[2] >> 8; bin[15] = a->d[2];
    bin[16] = a->d[1] >> 56; bin[17] = a->d[1] >> 48; bin[18] = a->d[1] >> 40; bin[19] = a->d[1] >> 32; bin[20] = a->d[1] >> 24; bin[21] = a->d[1] >> 16; bin[22] = a->d[1] >> 8; bin[23] = a->d[1];
    bin[24] = a->d[0] >> 56; bin[25] = a->d[0] >> 48; bin[26] = a->d[0] >> 40; bin[27] = a->d[0] >> 32; bin[28] = a->d[0] >> 24; bin[29] = a->d[0] >> 16; bin[30] = a->d[0] >> 8; bin[31] = a->d[0];
}

SECP256K1_INLINE static int secp256k1_scalar_is_zero(const secp256k1_scalar *a) {
    return (a->d[0] | a->d[1] | a->d[2] | a->d[3]) == 0;
}

static void secp256k1_scalar_negate(secp256k1_scalar *r, const secp256k1_scalar *a) {
    uint64_t nonzero = 0xFFFFFFFFFFFFFFFFULL * (secp256k1_scalar_is_zero(a) == 0);
    uint128_t t = (uint128_t)(~a->d[0]) + SECP256K1_N_0 + 1;
    r->d[0] = t & nonzero; t >>= 64;
    t += (uint128_t)(~a->d[1]) + SECP256K1_N_1;
    r->d[1] = t & nonzero; t >>= 64;
    t += (uint128_t)(~a->d[2]) + SECP256K1_N_2;
    r->d[2] = t & nonzero; t >>= 64;
    t += (uint128_t)(~a->d[3]) + SECP256K1_N_3;
    r->d[3] = t & nonzero;
}

SECP256K1_INLINE static int secp256k1_scalar_is_one(const secp256k1_scalar *a) {
    return ((a->d[0] ^ 1) | a->d[1] | a->d[2] | a->d[3]) == 0;
}

static int secp256k1_scalar_is_high(const secp256k1_scalar *a) {
    int yes = 0;
    int no = 0;
    no |= (a->d[3] < SECP256K1_N_H_3);
    yes |= (a->d[3] > SECP256K1_N_H_3) & ~no;
    no |= (a->d[2] < SECP256K1_N_H_2) & ~yes; /* No need for a > check. */
    no |= (a->d[1] < SECP256K1_N_H_1) & ~yes;
    yes |= (a->d[1] > SECP256K1_N_H_1) & ~no;
    yes |= (a->d[0] > SECP256K1_N_H_0) & ~no;
    return yes;
}

static int secp256k1_scalar_cond_negate(secp256k1_scalar *r, int flag) {
    /* If we are flag = 0, mask = 00...00 and this is a no-op;
     * if we are flag = 1, mask = 11...11 and this is identical to secp256k1_scalar_negate */
    uint64_t mask = !flag - 1;
    uint64_t nonzero = (secp256k1_scalar_is_zero(r) != 0) - 1;
    uint128_t t = (uint128_t)(r->d[0] ^ mask) + ((SECP256K1_N_0 + 1) & mask);
    r->d[0] = t & nonzero; t >>= 64;
    t += (uint128_t)(r->d[1] ^ mask) + (SECP256K1_N_1 & mask);
    r->d[1] = t & nonzero; t >>= 64;
    t += (uint128_t)(r->d[2] ^ mask) + (SECP256K1_N_2 & mask);
    r->d[2] = t & nonzero; t >>= 64;
    t += (uint128_t)(r->d[3] ^ mask) + (SECP256K1_N_3 & mask);
    r->d[3] = t & nonzero;
    return 2 * (mask == 0) - 1;
}

/* Inspired by the macros in OpenSSL's crypto/bn/asm/x86_64-gcc.c. */

/** Add a*b to the number defined by (c0,c1,c2). c2 must never overflow. */
#define muladd(a,b) { \
    uint64_t tl, th; \
    { \
        uint128_t t = (uint128_t)a * b; \
        th = t >> 64;         /* at most 0xFFFFFFFFFFFFFFFE */ \
        tl = t; \
    } \
    c0 += tl;                 /* overflow is handled on the next line */ \
    th += (c0 < tl) ? 1 : 0;  /* at most 0xFFFFFFFFFFFFFFFF */ \
    c1 += th;                 /* overflow is handled on the next line */ \
    c2 += (c1 < th) ? 1 : 0;  /* never overflows by contract (verified in the next line) */ \
    VERIFY_CHECK((c1 >= th) || (c2 != 0)); \
}

/** Add a*b to the number defined by (c0,c1). c1 must never overflow. */
#define muladd_fast(a,b) { \
    uint64_t tl, th; \
    { \
        uint128_t t = (uint128_t)a * b; \
        th = t >> 64;         /* at most 0xFFFFFFFFFFFFFFFE */ \
        tl = t; \
    } \
    c0 += tl;                 /* overflow is handled on the next line */ \
    th += (c0 < tl) ? 1 : 0;  /* at most 0xFFFFFFFFFFFFFFFF */ \
    c1 += th;                 /* never overflows by contract (verified in the next line) */ \
    VERIFY_CHECK(c1 >= th); \
}

/** Add 2*a*b to the number defined by (c0,c1,c2). c2 must never overflow. */
#define muladd2(a,b) { \
    uint64_t tl, th, th2, tl2; \
    { \
        uint128_t t = (uint128_t)a * b; \
        th = t >> 64;               /* at most 0xFFFFFFFFFFFFFFFE */ \
        tl = t; \
    } \
    th2 = th + th;                  /* at most 0xFFFFFFFFFFFFFFFE (in case th was 0x7FFFFFFFFFFFFFFF) */ \
    c2 += (th2 < th) ? 1 : 0;       /* never overflows by contract (verified the next line) */ \
    VERIFY_CHECK((th2 >= th) || (c2 != 0)); \
    tl2 = tl + tl;                  /* at most 0xFFFFFFFFFFFFFFFE (in case the lowest 63 bits of tl were 0x7FFFFFFFFFFFFFFF) */ \
    th2 += (tl2 < tl) ? 1 : 0;      /* at most 0xFFFFFFFFFFFFFFFF */ \
    c0 += tl2;                      /* overflow is handled on the next line */ \
    th2 += (c0 < tl2) ? 1 : 0;      /* second overflow is handled on the next line */ \
    c2 += (c0 < tl2) & (th2 == 0);  /* never overflows by contract (verified the next line) */ \
    VERIFY_CHECK((c0 >= tl2) || (th2 != 0) || (c2 != 0)); \
    c1 += th2;                      /* overflow is handled on the next line */ \
    c2 += (c1 < th2) ? 1 : 0;       /* never overflows by contract (verified the next line) */ \
    VERIFY_CHECK((c1 >= th2) || (c2 != 0)); \
}

/** Add a to the number defined by (c0,c1,c2). c2 must never overflow. */
#define sumadd(a) { \
    unsigned int over; \
    c0 += (a);                  /* overflow is handled on the next line */ \
    over = (c0 < (a)) ? 1 : 0; \
    c1 += over;                 /* overflow is handled on the next line */ \
    c2 += (c1 < over) ? 1 : 0;  /* never overflows by contract */ \
}

/** Add a to the number defined by (c0,c1). c1 must never overflow, c2 must be zero. */
#define sumadd_fast(a) { \
    c0 += (a);                 /* overflow is handled on the next line */ \
    c1 += (c0 < (a)) ? 1 : 0;  /* never overflows by contract (verified the next line) */ \
    VERIFY_CHECK((c1 != 0) | (c0 >= (a))); \
    VERIFY_CHECK(c2 == 0); \
}

/** Extract the lowest 64 bits of (c0,c1,c2) into n, and left shift the number 64 bits. */
#define extract(n) { \
    (n) = c0; \
    c0 = c1; \
    c1 = c2; \
    c2 = 0; \
}

/** Extract the lowest 64 bits of (c0,c1,c2) into n, and left shift the number 64 bits. c2 is required to be zero. */
#define extract_fast(n) { \
    (n) = c0; \
    c0 = c1; \
    c1 = 0; \
    VERIFY_CHECK(c2 == 0); \
}


static void secp256k1_scalar_reduce_512(secp256k1_scalar *r, const uint64_t *l) {
#ifdef USE_ASM_X86_64
  __asm__ __volatile__(
    /* Preload. */
    "movq 32(%%rsi), %%r11\n" 
    "movq 40(%%rsi), %%r12\n" 
    "movq 0(%%rsi), %%rbx\n"  
    "movq %3, %%rax\n"
    "movq %%rax, %%r10\n"
    "xor %%ecx, %%ecx\n"  
    "xorq %%r15, %%r15\n"
    "xorq %%r9, %%r9\n"
    "xorq %%r8, %%r8\n"
    "mulq %%r11\n"
    "movq %4, %%r14\n" 
    "addq %%rax, %%rbx\n" /*q0 into rbx*/
    "adcq %%rdx, %%rcx\n"
    "addq 8(%%rsi), %%rcx\n" 
    "adcq %%r9, %%r15\n"
    "movq %%r10, %%rax\n"
    "movq 48(%%rsi), %%r13\n" 
    "mulq %%r12\n"
    "addq %%rax, %%rcx\n" /*q1 stored to rcx*/
    "adcq %%rdx, %%r15\n"
    "adcq %%r9, %%r8\n"
    "movq %%r14, %%rax\n" 
    "mulq %%r11\n"
    "addq %%rax, %%rcx\n"
    "adcq %%rdx, %%r15\n"
    "adcq %%r9, %%r8\n"
    "addq 16(%%rsi), %%r15\n"
    "adcq %%r9, %%r8\n"
    "adcq %%r9, %%r9\n"
    "movq %%r10, %%rax\n"
    "mulq %%r13\n"
    "addq %%rax, %%r15\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    "movq %%r14, %%rax\n"
    "movq 56(%%rsi), %%r14\n" 
    "mulq %%r12\n"
    "addq %%rax, %%r15\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    "movq %%r10, %%rax\n"
    "movq $0, %%r10\n"
    "addq %%r11, %%r15\n" /*q2 into r15*/
    "adcq $0, %%r8\n"
    "adcq $0, %%r9\n"
    "addq 24(%%rsi), %%r8\n"
    "adcq $0, %%r9\n"
    "adcq %%r10, %%r10\n"
    "mulq %%r14\n"
    "movq %4, %%rsi\n"  
    "addq %%rax, %%r8\n"
    "movq %%rsi, %%rax\n"  
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    "mulq %%r13\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    "addq %%r8, %%r12\n" /* q3 into r12*/
    "adcq $0, %%r9\n"
    "adcq $0, %%r10\n"
    "movq %%rsi, %%rax\n" 
    "xorq %%r8, %%r8\n"
    "mulq %%r14\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq %%r8, %%r8\n"
    "addq %%r9, %%r13\n" /*q4 into r13*/
    "adcq $0, %%r10\n" 
    "adcq $0, %%r8\n" 
    "addq %%r14, %%r10\n" /* q5 into r10 */ 
    "movq %3, %%rax\n"
    "movq %%rax, %%r9\n"
    "adcq $0, %%r8\n" /*q6 into r8*/
  
/* %q5 input for second operation is %q0 output from first / RBX as the connecting link
    %q6 input for second operation is %q1 output from first / RCX as the connecting link
    %q7 input for second operation is %q2 output from first / R15 as the connecting link
    %q8 input for second operation is %q3 output from first / R12 as the connecting link
    %q9  input for second operation is %q4 output from first / R13 as the connecting link*
    %q10 input for second operation is %q5 output from first / R10 as the connecting link*
    %q11 input for second operation is %q6 output from first  / R8 as the connecting link */    
    
    /* Reduce 385 bits into 258. */

    "mulq %%r13\n"
    "xorq %%r14, %%r14\n"
    "xorq %%r11, %%r11\n"
    "addq %%rax, %%rbx\n" /* q0 output*/
    "adcq %%rdx, %%r14\n"
    "addq %%rcx, %%r14\n" 
    "adcq %%r11, %%r11\n"
    "xor %%ecx, %%ecx\n"  
    "movq %%r9, %%rax\n"
    "mulq %%r10\n"
    "addq %%rax, %%r14\n"
    "adcq %%rdx, %%r11\n"
    "adcq %%rcx, %%rcx\n"
    "movq %%rsi, %%rax\n"
    "mulq %%r13\n"
    "addq %%rax, %%r14\n" /* q1 output */
    "movq %%r9, %%rax\n"
    "adcq %%rdx, %%r11\n"
    "adcq $0, %%rcx\n"
    "xorq %%r9, %%r9\n"
    "addq %%r15, %%r11\n" 
    "adcq %%r9, %%rcx\n"
    "movq %%rax, %%r15\n"
    "adcq %%r9, %%r9\n"
    "mulq %%r8\n"
    "addq %%rax, %%r11\n"
    "adcq %%rdx, %%rcx\n"
    "adcq $0, %%r9\n"
    "movq %%rsi, %%rax\n"
    "mulq %%r10\n"
    "addq %%rax, %%r11\n"
    "adcq %%rdx, %%rcx\n"
    "adcq $0, %%r9\n"
    "addq %%r13, %%r11\n" /* q2 output */
    "adcq $0, %%rcx\n"
    "adcq $0, %%r9\n"
    "addq %%r12, %%rcx\n"
    "adcq $0, %%r9\n"
    "movq %%rsi, %%rax\n"
    "mulq %%r8\n"
    "addq %%rax, %%rcx\n"
    "adcq %%rdx, %%r9\n"
    "addq %%r10, %%rcx\n"    /* q3 output */
    "adcq $0, %%r9\n"
    "addq %%r8, %%r9\n" /* q4 output */
    
/* %q1 input for next operation is %q0 output from prior / RBX as the connecting link
    %q2 input for next operation is %q1 output from prior / R14 as the connecting link 
    %q3 input for next operation is %q2 output from prior / R11 as the connecting link  
    %q4 input for next operation is %q3 output from prior / RCX as the connecting link
    %q5 input for next operation is %q4 output from prior / R9 as the connecting link   */
        
    /* Reduce 258 bits into 256. */

    "movq %%r15, %%rax\n"
    "mulq %%r9\n"   
    "addq %%rbx, %%rax\n"
    "adcq $0, %%rdx\n"
    "movq %%rax, %%r8\n"  /* 0(q2) output */
    "movq %%rdx, %%r12\n" 
    "xorq %%r13, %%r13\n"
    "addq %%r14, %%r12\n"
    "adcq %%r13, %%r13\n"
    "movq %%rsi, %%rax\n"
    "mulq %%r9\n"
    "addq %%rax, %%r12\n" /* 8(q2) output */
    "adcq %%rdx, %%r13\n" 
    "xor %%ebx, %%ebx\n"
    "addq %%r9, %%r13\n"
    "adcq %%rbx, %%rbx\n"
    "movq $0xffffffffffffffff, %%r14\n"
    "addq %%r11, %%r13\n" /* 16(q2) output */
    "movq $0, %%r11\n"
    "adcq $0, %%rbx\n"
    "addq %%rcx, %%rbx\n"  /* 24(q2) output */
    "adcq $0, %%r11\n" /* c  output */

    
    /*FINAL REDUCTION*/
    
/*    r8 carries ex 0(%%rdi), 
       r12 carries ex 8(%%rdi),
       r13 carries ex 16(%%rdi), 
       rbx carries ex 24(%%rdi)
       r11 carries c */
    "movq $0xbaaedce6af48a03b,%%r9\n"
    "movq $0xbaaedce6af48a03a,%%rcx\n"
    "movq $0xbfd25e8cd0364140,%%r10\n"
    "cmp   %%r14 ,%%rbx\n"
    "setne %%dl\n"
    "cmp   $0xfffffffffffffffd,%%r13\n"
    "setbe %%al\n"
    "or     %%eax,%%edx\n"
    "cmp  %%rcx,%%r12\n"
    "setbe %%cl\n"
    "or     %%edx,%%ecx\n"
    "cmp  %%r9,%%r12\n"
    "movzbl %%dl,%%edx\n"
    "seta  %%r9b\n"
    "cmp  %%r10,%%r8\n"
    "movzbl %%cl,%%ecx\n"
    "seta  %%r10b\n"
    "not   %%ecx\n"
    "not   %%edx\n"
    "or     %%r10d,%%r9d\n"
    "movzbl %%r9b,%%r9d\n"
    "and   %%r9d,%%ecx\n"
    "xor    %%r9d,%%r9d\n"
    "cmp   %%r14,%%r13\n"
    "sete  %%r9b\n"
    "xor   %%r10d,%%r10d\n"
    "and   %%r9d,%%edx\n"
    "or     %%edx,%%ecx\n"
    "xor   %%edx,%%edx\n"
    "add  %%ecx,%%r11d\n"
    "imulq %%r11,%%r15\n"
    "addq  %%r15,%%r8\n"
    "adcq  %%rdx,%%r10\n"  
    "imulq %%r11,%%rsi\n"
    "xorq %%r15,%%r15\n"
    "xor   %%eax,%%eax\n"
    "movq  %%r8,0(%%rdi)\n"
    "xor   %%edx,%%edx\n"
    "addq %%r12,%%rsi\n"
    "adcq %%rdx,%%rdx\n" 
    "addq %%rsi,%%r10\n"
    "movq %%r10,8(%%rdi)\n"
    "adcq %%rdx,%%r15\n"
    "addq %%r11,%%r13\n"
    "adcq %%rax,%%rax\n" 
    "addq %%r15,%%r13\n"
    "movq %%r13,16(%%rdi)\n"
    "adcq $0,%%rax\n"
    "addq %%rbx,%%rax\n"
    "movq %%rax,24(%%rdi)\n"

    : "=D"(r)
    : "S"(l), "D"(r), "n"(SECP256K1_N_C_0), "n"(SECP256K1_N_C_1)
    : "rax", "rbx", "rcx", "rdx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "cc", "memory"); 
#else
    uint128_t c;
    uint64_t c0, c1, c2;
    uint64_t n0 = l[4], n1 = l[5], n2 = l[6], n3 = l[7];
    uint64_t m0, m1, m2, m3, m4, m5;
    uint32_t m6;
    uint64_t p0, p1, p2, p3;
    uint32_t p4;

    /* Reduce 512 bits into 385. */
    /* m[0..6] = l[0..3] + n[0..3] * SECP256K1_N_C. */
    c0 = l[0]; c1 = 0; c2 = 0;
    muladd_fast(n0, SECP256K1_N_C_0);
    extract_fast(m0);
    sumadd_fast(l[1]);
    muladd(n1, SECP256K1_N_C_0);
    muladd(n0, SECP256K1_N_C_1);
    extract(m1);
    sumadd(l[2]);
    muladd(n2, SECP256K1_N_C_0);
    muladd(n1, SECP256K1_N_C_1);
    sumadd(n0);
    extract(m2);
    sumadd(l[3]);
    muladd(n3, SECP256K1_N_C_0);
    muladd(n2, SECP256K1_N_C_1);
    sumadd(n1);
    extract(m3);
    muladd(n3, SECP256K1_N_C_1);
    sumadd(n2);
    extract(m4);
    sumadd_fast(n3);
    extract_fast(m5);
    VERIFY_CHECK(c0 <= 1);
    m6 = c0;

    /* Reduce 385 bits into 258. */
    /* p[0..4] = m[0..3] + m[4..6] * SECP256K1_N_C. */
    c0 = m0; c1 = 0; c2 = 0;
    muladd_fast(m4, SECP256K1_N_C_0);
    extract_fast(p0);
    sumadd_fast(m1);
    muladd(m5, SECP256K1_N_C_0);
    muladd(m4, SECP256K1_N_C_1);
    extract(p1);
    sumadd(m2);
    muladd(m6, SECP256K1_N_C_0);
    muladd(m5, SECP256K1_N_C_1);
    sumadd(m4);
    extract(p2);
    sumadd_fast(m3);
    muladd_fast(m6, SECP256K1_N_C_1);
    sumadd_fast(m5);
    extract_fast(p3);
    p4 = c0 + m6;
    VERIFY_CHECK(p4 <= 2);

    /* Reduce 258 bits into 256. */
    /* r[0..3] = p[0..3] + p[4] * SECP256K1_N_C. */
    c = p0 + (uint128_t)SECP256K1_N_C_0 * p4;
    r->d[0] = c & 0xFFFFFFFFFFFFFFFFULL; c >>= 64;
    c += p1 + (uint128_t)SECP256K1_N_C_1 * p4;
    r->d[1] = c & 0xFFFFFFFFFFFFFFFFULL; c >>= 64;
    c += p2 + (uint128_t)p4;
    r->d[2] = c & 0xFFFFFFFFFFFFFFFFULL; c >>= 64;
    c += p3;
    r->d[3] = c & 0xFFFFFFFFFFFFFFFFULL; c >>= 64;
    
        /* Final reduction of r. */
    secp256k1_scalar_reduce(r, c + secp256k1_scalar_check_overflow(r));
    
#endif

    /* Final reduction of r. */
    /*secp256k1_scalar_reduce(r, c + secp256k1_scalar_check_overflow(r));  Only active in the non-asm version. The asm version doesn't need it anymore*/
}

static void secp256k1_scalar_mul_512(uint64_t l[8], const secp256k1_scalar *a, const secp256k1_scalar *b) {
#ifdef USE_ASM_X86_64
    const uint64_t *pb = b->d;
    __asm__ __volatile__(
    /* Preload */
    "movq 0(%%rdi), %%r15\n"
    "movq 8(%%rdi), %%rbx\n"
    "movq 16(%%rdi), %%rcx\n"
    "movq 0(%%rdx), %%r11\n"
    "movq 8(%%rdx), %%r12\n"
    "movq 16(%%rdx), %%r13\n"
    "movq 24(%%rdx), %%r14\n"
    /* (rax,rdx) = a0 * b0 */
    "movq %%r15, %%rax\n"
    "mulq %%r11\n"
    /* Extract l0 */
    "movq %%rax, 0(%%rsi)\n"
    /* (r8,r9,r10) = (rdx) */
    "movq %%rdx, %%r8\n"
    "xorq %%r9, %%r9\n"
    "xorq %%r10, %%r10\n"
    /* (r8,r9,r10) += a0 * b1 */
    "movq %%r15, %%rax\n"
    "mulq %%r12\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    /* (r8,r9,r10) += a1 * b0 */
    "movq %%rbx, %%rax\n"
    "mulq %%r11\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    /* Extract l1 */
    "movq %%r8, 8(%%rsi)\n"
    "xorq %%r8, %%r8\n"
    /* (r9,r10,r8) += a0 * b2 */
    "movq %%r15, %%rax\n"
    "mulq %%r13\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    /* (r9,r10,r8) += a1 * b1 */
    "movq %%rbx, %%rax\n"
    "mulq %%r12\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    /* (r9,r10,r8) += a2 * b0 */
    "movq %%rcx, %%rax\n"
    "mulq %%r11\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    /* Extract l2 */
    "movq %%r9, 16(%%rsi)\n"
    "xorq %%r9, %%r9\n"
    /* (r10,r8,r9) += a0 * b3 */
    "movq %%r15, %%rax\n"
    "mulq %%r14\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    /* Preload a3 */
    "movq 24(%%rdi), %%r15\n"
    /* (r10,r8,r9) += a1 * b2 */
    "movq %%rbx, %%rax\n"
    "mulq %%r13\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    /* (r10,r8,r9) += a2 * b1 */
    "movq %%rcx, %%rax\n"
    "mulq %%r12\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    /* (r10,r8,r9) += a3 * b0 */
    "movq %%r15, %%rax\n"
    "mulq %%r11\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    /* Extract l3 */
    "movq %%r10, 24(%%rsi)\n"
    "xorq %%r10, %%r10\n"
    /* (r8,r9,r10) += a1 * b3 */
    "movq %%rbx, %%rax\n"
    "mulq %%r14\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    /* (r8,r9,r10) += a2 * b2 */
    "movq %%rcx, %%rax\n"
    "mulq %%r13\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    /* (r8,r9,r10) += a3 * b1 */
    "movq %%r15, %%rax\n"
    "mulq %%r12\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    /* Extract l4 */
    "movq %%r8, 32(%%rsi)\n"
    "xorq %%r8, %%r8\n"
    /* (r9,r10,r8) += a2 * b3 */
    "movq %%rcx, %%rax\n"
    "mulq %%r14\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    /* (r9,r10,r8) += a3 * b2 */
    "movq %%r15, %%rax\n"
    "mulq %%r13\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    /* Extract l5 */
    "movq %%r9, 40(%%rsi)\n"
    /* (r10,r8) += a3 * b3 */
    "movq %%r15, %%rax\n"
    "mulq %%r14\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    /* Extract l6 */
    "movq %%r10, 48(%%rsi)\n"
    /* Extract l7 */
    "movq %%r8, 56(%%rsi)\n"
    : "+d"(pb)
    : "S"(l), "D"(a->d)
    : "rax", "rbx", "rcx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "cc", "memory");
#else
    /* 160 bit accumulator. */
    uint64_t c0 = 0, c1 = 0;
    uint32_t c2 = 0;

    /* l[0..7] = a[0..3] * b[0..3]. */
    muladd_fast(a->d[0], b->d[0]);
    extract_fast(l[0]);
    muladd(a->d[0], b->d[1]);
    muladd(a->d[1], b->d[0]);
    extract(l[1]);
    muladd(a->d[0], b->d[2]);
    muladd(a->d[1], b->d[1]);
    muladd(a->d[2], b->d[0]);
    extract(l[2]);
    muladd(a->d[0], b->d[3]);
    muladd(a->d[1], b->d[2]);
    muladd(a->d[2], b->d[1]);
    muladd(a->d[3], b->d[0]);
    extract(l[3]);
    muladd(a->d[1], b->d[3]);
    muladd(a->d[2], b->d[2]);
    muladd(a->d[3], b->d[1]);
    extract(l[4]);
    muladd(a->d[2], b->d[3]);
    muladd(a->d[3], b->d[2]);
    extract(l[5]);
    muladd_fast(a->d[3], b->d[3]);
    extract_fast(l[6]);
    VERIFY_CHECK(c1 == 0);
    l[7] = c0;
#endif
}

static void secp256k1_scalar_sqr_512(uint64_t l[8], const secp256k1_scalar *a) {
#ifdef USE_ASM_X86_64
    __asm__ __volatile__(
    /* Preload */
    "movq 0(%%rdi), %%r11\n"
    "movq 8(%%rdi), %%r12\n"
    "movq 16(%%rdi), %%r13\n"
    "movq 24(%%rdi), %%r14\n"
    /* (rax,rdx) = a0 * a0 */
    "movq %%r11, %%rax\n"
    "mulq %%r11\n"
    /* Extract l0 */
    "movq %%rax, 0(%%rsi)\n"
    /* (r8,r9,r10) = (rdx,0) */
    "movq %%rdx, %%r8\n"
    "xorq %%r9, %%r9\n"
    "xorq %%r10, %%r10\n"
    /* (r8,r9,r10) += 2 * a0 * a1 */
    "movq %%r11, %%rax\n"
    "mulq %%r12\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    /* Extract l1 */
    "movq %%r8, 8(%%rsi)\n"
    "xorq %%r8, %%r8\n"
    /* (r9,r10,r8) += 2 * a0 * a2 */
    "movq %%r11, %%rax\n"
    "mulq %%r13\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    /* (r9,r10,r8) += a1 * a1 */
    "movq %%r12, %%rax\n"
    "mulq %%r12\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    /* Extract l2 */
    "movq %%r9, 16(%%rsi)\n"
    "xorq %%r9, %%r9\n"
    /* (r10,r8,r9) += 2 * a0 * a3 */
    "movq %%r11, %%rax\n"
    "mulq %%r14\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    /* (r10,r8,r9) += 2 * a1 * a2 */
    "movq %%r12, %%rax\n"
    "mulq %%r13\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    /* Extract l3 */
    "movq %%r10, 24(%%rsi)\n"
    "xorq %%r10, %%r10\n"
    /* (r8,r9,r10) += 2 * a1 * a3 */
    "movq %%r12, %%rax\n"
    "mulq %%r14\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    /* (r8,r9,r10) += a2 * a2 */
    "movq %%r13, %%rax\n"
    "mulq %%r13\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    /* Extract l4 */
    "movq %%r8, 32(%%rsi)\n"
    "xorq %%r8, %%r8\n"
    /* (r9,r10,r8) += 2 * a2 * a3 */
    "movq %%r13, %%rax\n"
    "mulq %%r14\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    /* Extract l5 */
    "movq %%r9, 40(%%rsi)\n"
    /* (r10,r8) += a3 * a3 */
    "movq %%r14, %%rax\n"
    "mulq %%r14\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    /* Extract l6 */
    "movq %%r10, 48(%%rsi)\n"
    /* Extract l7 */
    "movq %%r8, 56(%%rsi)\n"
    :
    : "S"(l), "D"(a->d)
    : "rax", "rdx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "cc", "memory");
#else
    /* 160 bit accumulator. */
    uint64_t c0 = 0, c1 = 0;
    uint32_t c2 = 0;

    /* l[0..7] = a[0..3] * b[0..3]. */
    muladd_fast(a->d[0], a->d[0]);
    extract_fast(l[0]);
    muladd2(a->d[0], a->d[1]);
    extract(l[1]);
    muladd2(a->d[0], a->d[2]);
    muladd(a->d[1], a->d[1]);
    extract(l[2]);
    muladd2(a->d[0], a->d[3]);
    muladd2(a->d[1], a->d[2]);
    extract(l[3]);
    muladd2(a->d[1], a->d[3]);
    muladd(a->d[2], a->d[2]);
    extract(l[4]);
    muladd2(a->d[2], a->d[3]);
    extract(l[5]);
    muladd_fast(a->d[3], a->d[3]);
    extract_fast(l[6]);
    VERIFY_CHECK(c1 == 0);
    l[7] = c0;
#endif
}

#undef sumadd
#undef sumadd_fast
#undef muladd
#undef muladd_fast
#undef muladd2
#undef extract
#undef extract_fast

static void secp256k1_scalar_mul(secp256k1_scalar *r, const secp256k1_scalar *a, const secp256k1_scalar *b) {
 #ifdef USE_ASM_X86_64
    uint64_t l[8];
    const uint64_t *pb = b->d;
    
    __asm__ __volatile__(
    /* Preload */
    "movq 0(%%rdi), %%r15\n"
    "movq 8(%%rdi), %%rbx\n"
    "movq 16(%%rdi), %%rcx\n"
    "movq 0(%%rdx), %%r11\n"
    "movq 8(%%rdx), %%r9\n"
    "movq 16(%%rdx), %%r10\n"
    "movq 24(%%rdx), %%r8\n"
    /* (rax,rdx) = a0 * b0 */
    "movq %%r15, %%rax\n"
    "mulq %%r11\n"
    /* Extract l0 */
    "movq %%rax, 0(%%rsi)\n"
    /* (r14,r12,r13) = (rdx) */
    "movq %%rdx, %%r14\n"
    "xorq %%r12, %%r12\n"
    "xorq %%r13, %%r13\n"
    /* (r14,r12,r13) += a0 * b1 */
    "movq %%r15, %%rax\n"
    "mulq %%r9\n"
    "addq %%rax, %%r14\n"
    "adcq %%rdx, %%r12\n"
    "movq %%rbx, %%rax\n"
    "adcq $0, %%r13\n"
    /* (r14,r12,r13) += a1 * b0 */
    "mulq %%r11\n"
    "addq %%rax, %%r14\n"
    "adcq %%rdx, %%r12\n"
    "adcq $0, %%r13\n"
    /* Extract l1 */
    "movq %%r14, 8(%%rsi)\n"
    "xorq %%r14, %%r14\n"
    /* (r12,r13,r14) += a0 * b2 */
    "movq %%r15, %%rax\n"
    "mulq %%r10\n"
    "addq %%rax, %%r12\n"
    "adcq %%rdx, %%r13\n"
    "movq %%rbx, %%rax\n"
    "adcq $0, %%r14\n"
    /* (r12,r13,r14) += a1 * b1 */
    "mulq %%r9\n"
    "addq %%rax, %%r12\n"
    "adcq %%rdx, %%r13\n"
    "movq %%rcx, %%rax\n"
    "adcq $0, %%r14\n"
    /* (r12,r13,r14) += a2 * b0 */
    "mulq %%r11\n"
    "addq %%rax, %%r12\n"
    "adcq %%rdx, %%r13\n"
    "adcq $0, %%r14\n"
    /* Extract l2 */
    "movq %%r12, 16(%%rsi)\n"
    "xorq %%r12, %%r12\n"
    /* (r13,r14,r12) += a0 * b3 */
    "movq %%r15, %%rax\n"
    "mulq %%r8\n"
    "addq %%rax, %%r13\n"
    "adcq %%rdx, %%r14\n"
    "adcq $0, %%r12\n"
    /* Preload a3 */
    "movq 24(%%rdi), %%r15\n"
    /* (r13,r14,r12) += a1 * b2 */
    "movq %%rbx, %%rax\n"
    "mulq %%r10\n"
    "addq %%rax, %%r13\n"
    "adcq %%rdx, %%r14\n"
    "movq %%rcx, %%rax\n"
    "adcq $0, %%r12\n"
    /* (r13,r14,r12) += a2 * b1 */
    "mulq %%r9\n"
    "addq %%rax, %%r13\n"
    "adcq %%rdx, %%r14\n"
    "movq %%r15, %%rax\n"
    "adcq $0, %%r12\n"
    /* (r13,r14,r12) += a3 * b0 */
    "mulq %%r11\n"
    "addq %%rax, %%r13\n"
    "adcq %%rdx, %%r14\n"
    "adcq $0, %%r12\n"
    /* Extract l3 */
    "movq %%r13, 24(%%rsi)\n"
    "xorq %%r13, %%r13\n"
    /* (r14,r12,r13) += a1 * b3 */
    "movq %%rbx, %%rax\n"
    "mulq %%r8\n"
    "addq %%rax, %%r14\n"
    "adcq %%rdx, %%r12\n"
    "movq %%rcx, %%rax\n"
    "adcq $0, %%r13\n"
    /* (r14,r12,r13) += a2 * b2 */
    "mulq %%r10\n"
    "addq %%rax, %%r14\n"
    "adcq %%rdx, %%r12\n"
    "movq %%r15, %%rax\n"
    "adcq $0, %%r13\n"
    /* (r14,r12,r13) += a3 * b1 */
    "mulq %%r9\n"
    "addq %%rax, %%r14\n"
    "adcq %%rdx, %%r12\n"
    "movq %%rcx, %%rax\n"
    "adcq $0, %%r13\n"
    /* Extract l4 */
   /* "movq %%r14, 32(%%rsi)\n"*/
    /* (r12,r13,r14) += a2 * b3 */
    "mulq %%r8\n"
    "movq %%r14, %%r11\n"
    "xorq %%r14, %%r14\n"
    "addq %%rax, %%r12\n"
    "movq %%r15, %%rax\n"
    "adcq %%rdx, %%r13\n"
    "adcq $0, %%r14\n"
    /* (r12,r13,r14) += a3 * b2 */
    "mulq %%r10\n"
    "addq %%rax, %%r12\n"
    "adcq %%rdx, %%r13\n"
    "movq %%r15, %%rax\n"
    "adcq $0, %%r14\n"
    /* Extract l5 */
    /*"movq %%r12, 40(%%rsi)\n"*/
    /* (r13,r14) += a3 * b3 */
    "mulq %%r8\n"
    "addq %%rax, %%r13\n"
    "adcq %%rdx, %%r14\n"
    /* Extract l6 */
    /*"movq %%r13, 48(%%rsi)\n"*/
    /* Extract l7 */
    /*"movq %%r14, 56(%%rsi)\n"*/
    : "+d"(pb)
    : "S"(l), "D"(a->d)
    : "rax", "rbx", "rcx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "cc", "memory");
    

      __asm__ __volatile__(
    /* Preload. */
  /*  "movq 32(%%rsi), %%r11\n" */
  /*  "movq 40(%%rsi), %%r12\n" */
   /*"movq 48(%%rsi), %%r13\n" */
  /*   "movq 56(%%rsi), %%r14\n" */
    "movq 0(%%rsi), %%rbx\n"  
    "movq %3, %%rax\n"
    "movq %%rax, %%r10\n"
    "xor %%ecx, %%ecx\n"  
    "xorq %%r15, %%r15\n"
    "xorq %%r9, %%r9\n"
    "xorq %%r8, %%r8\n"
    "mulq %%r11\n"
    "addq %%rax, %%rbx\n" /*q0 into rbx*/
    "adcq %%rdx, %%rcx\n"
    "addq 8(%%rsi), %%rcx\n" 
    "movq %%r10, %%rax\n"
    "adcq %%r9, %%r15\n"
    "mulq %%r12\n"
    "addq %%rax, %%rcx\n" /*q1 stored to rcx*/
    "adcq %%rdx, %%r15\n"
    "movq %4, %%rax\n" 
    "adcq %%r9, %%r8\n"
    "mulq %%r11\n"
    "addq %%rax, %%rcx\n"
    "adcq %%rdx, %%r15\n"
    "adcq %%r9, %%r8\n"
    "addq 16(%%rsi), %%r15\n"
    "adcq %%r9, %%r8\n"
    "movq %%r10, %%rax\n"
    "adcq %%r9, %%r9\n"
    "mulq %%r13\n"
    "addq %%rax, %%r15\n"
    "adcq %%rdx, %%r8\n"
    "movq %4, %%rax\n"
    "adcq $0, %%r9\n"
    "mulq %%r12\n"
    "addq %%rax, %%r15\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    "movq %%r10, %%rax\n"
    "movq $0, %%r10\n"
    "addq %%r11, %%r15\n" /*q2 into r15*/
    "adcq $0, %%r8\n"
    "adcq $0, %%r9\n"
    "addq 24(%%rsi), %%r8\n"
    "adcq $0, %%r9\n"
    "adcq %%r10, %%r10\n"
    "mulq %%r14\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "movq %4, %%rax\n"  
    "movq %%rax, %%rsi\n"  
    "adcq $0, %%r10\n"
    "mulq %%r13\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    "addq %%r8, %%r12\n" /* q3 into r12*/
    "adcq $0, %%r9\n"
    "movq $0, %%r8\n"
    "movq %%rsi, %%rax\n" 
    "adcq $0, %%r10\n"
    "mulq %%r14\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq %%r8, %%r8\n"
    "addq %%r9, %%r13\n" /*q4 into r13*/
    "adcq $0, %%r10\n" 
    "adcq $0, %%r8\n" 
    "addq %%r14, %%r10\n" /* q5 into r10 */ 
    "movq %3, %%rax\n"
    "movq %%rax, %%r9\n"
    "adcq $0, %%r8\n" /*q6 into r8*/
  
/* %q5 input for second operation is %q0 output from first / RBX as the connecting link
    %q6 input for second operation is %q1 output from first / RCX as the connecting link
    %q7 input for second operation is %q2 output from first / R15 as the connecting link
    %q8 input for second operation is %q3 output from first / R12 as the connecting link
    %q9  input for second operation is %q4 output from first / R13 as the connecting link*
    %q10 input for second operation is %q5 output from first / R10 as the connecting link*
    %q11 input for second operation is %q6 output from first  / R8 as the connecting link */    
    
    /* Reduce 385 bits into 258. */

    "mulq %%r13\n"
    "xorq %%r14, %%r14\n"
    "xorq %%r11, %%r11\n"
    "addq %%rax, %%rbx\n" /* q0 output*/
    "adcq %%rdx, %%r14\n"
    "addq %%rcx, %%r14\n" 
    "mov $0, %%ecx\n"  
    "movq %%r9, %%rax\n"
    "adcq %%r11, %%r11\n"
    "mulq %%r10\n"
    "addq %%rax, %%r14\n"
    "adcq %%rdx, %%r11\n"
    "movq %%rsi, %%rax\n"
    "adcq %%rcx, %%rcx\n"
    "mulq %%r13\n"
    "addq %%rax, %%r14\n" /* q1 output */
    "movq %%r9, %%rax\n"
    "adcq %%rdx, %%r11\n"
    "adcq $0, %%rcx\n"
    "xorq %%r9, %%r9\n"
    "addq %%r15, %%r11\n" 
    "adcq %%r9, %%rcx\n"
    "movq %%rax, %%r15\n"
    "adcq %%r9, %%r9\n"
    "mulq %%r8\n"
    "addq %%rax, %%r11\n"
    "adcq %%rdx, %%rcx\n"
    "movq %%rsi, %%rax\n"
    "adcq $0, %%r9\n"
    "mulq %%r10\n"
    "addq %%rax, %%r11\n"
    "adcq %%rdx, %%rcx\n"
    "adcq $0, %%r9\n"
    "addq %%r13, %%r11\n" /* q2 output */
    "adcq $0, %%rcx\n"
    "adcq $0, %%r9\n"
    "addq %%r12, %%rcx\n"
    "movq %%rsi, %%rax\n"
    "adcq $0, %%r9\n"
    "mulq %%r8\n"
    "addq %%rax, %%rcx\n"
    "adcq %%rdx, %%r9\n"
    "addq %%r10, %%rcx\n"    /* q3 output */
    "adcq $0, %%r9\n"
    "movq %%r15, %%rax\n"
    "addq %%r8, %%r9\n" /* q4 output */
    
/* %q1 input for next operation is %q0 output from prior / RBX as the connecting link
    %q2 input for next operation is %q1 output from prior / R14 as the connecting link 
    %q3 input for next operation is %q2 output from prior / R11 as the connecting link  
    %q4 input for next operation is %q3 output from prior / RCX as the connecting link
    %q5 input for next operation is %q4 output from prior / R9 as the connecting link   */
        
    /* Reduce 258 bits into 256. */

    "mulq %%r9\n"   
    "addq %%rbx, %%rax\n"
    "adcq $0, %%rdx\n"
    "movq %%rax, %%r8\n"  /* 0(q2) output */
    "movq %%rdx, %%r12\n" 
    "xorq %%r13, %%r13\n"
    "addq %%r14, %%r12\n"
    "movq %%rsi, %%rax\n"
    "adcq %%r13, %%r13\n"
    "mulq %%r9\n"
    "addq %%rax, %%r12\n" /* 8(q2) output */
    "adcq %%rdx, %%r13\n" 
    "xor %%ebx, %%ebx\n"
    "addq %%r9, %%r13\n"
    "adcq %%rbx, %%rbx\n"
    "movq $0xffffffffffffffff, %%r14\n"
    "addq %%r11, %%r13\n" /* 16(q2) output */
    "movq $0, %%r11\n"
    "adcq $0, %%rbx\n"
    "addq %%rcx, %%rbx\n"  /* 24(q2) output */
    "adcq $0, %%r11\n" /* c  output */

    
    /*FINAL REDUCTION*/
    
/*    r8 carries ex 0(%%rdi), 
       r12 carries ex 8(%%rdi),
       r13 carries ex 16(%%rdi), 
       rbx carries ex 24(%%rdi)
       r11 carries c */
    "movq $0xbaaedce6af48a03b,%%r9\n"
    "movq $0xbaaedce6af48a03a,%%rcx\n"
    "movq $0xbfd25e8cd0364140,%%r10\n"
    "cmp   %%r14 ,%%rbx\n"
    "setne %%dl\n"
    "cmp   $0xfffffffffffffffd,%%r13\n"
    "setbe %%al\n"
    "or     %%eax,%%edx\n"
    "cmp  %%rcx,%%r12\n"
    "setbe %%cl\n"
    "or     %%edx,%%ecx\n"
    "cmp  %%r9,%%r12\n"
    "movzbl %%dl,%%edx\n"
    "seta  %%r9b\n"
    "cmp  %%r10,%%r8\n"
    "movzbl %%cl,%%ecx\n"
    "seta  %%r10b\n"
    "not   %%ecx\n"
    "not   %%edx\n"
    "or     %%r10d,%%r9d\n"
    "movzbl %%r9b,%%r9d\n"
    "and   %%r9d,%%ecx\n"
    "xor    %%r9d,%%r9d\n"
    "cmp   %%r14,%%r13\n"
    "sete  %%r9b\n"
    "xor   %%r10d,%%r10d\n"
    "and   %%r9d,%%edx\n"
    "or     %%edx,%%ecx\n"
    "xor   %%edx,%%edx\n"
    "add  %%ecx,%%r11d\n"
    "imulq %%r11,%%r15\n"
    "addq  %%r15,%%r8\n"
    "adcq  %%rdx,%%r10\n"  
    "imulq %%r11,%%rsi\n"
    "xorq %%r15,%%r15\n"
    "xor   %%eax,%%eax\n"
    "movq  %%r8,0(%q2)\n"
    "xor   %%edx,%%edx\n"
    "addq %%r12,%%rsi\n"
    "adcq %%rdx,%%rdx\n" 
    "addq %%rsi,%%r10\n"
    "movq %%r10,8(%q2)\n"
    "adcq %%rdx,%%r15\n"
    "addq %%r11,%%r13\n"
    "adcq %%rax,%%rax\n" 
    "addq %%r15,%%r13\n"
    "movq %%r13,16(%q2)\n"
    "adcq $0,%%rax\n"
    "addq %%rbx,%%rax\n"
    "movq %%rax,24(%q2)\n"
    : "=D"(r)
    : "S"(l), "D"(r), "n"(SECP256K1_N_C_0), "n"(SECP256K1_N_C_1)
    : "rax", "rbx", "rcx", "rdx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "cc", "memory"); 
     

#else
    uint64_t l[8];
    secp256k1_scalar_mul_512(l, a, b);
    secp256k1_scalar_reduce_512(r, l);
#endif   
}

static int secp256k1_scalar_shr_int(secp256k1_scalar *r, int n) {
    int ret;
    VERIFY_CHECK(n > 0);
    VERIFY_CHECK(n < 16);
    ret = r->d[0] & ((1 << n) - 1);
    r->d[0] = (r->d[0] >> n) + (r->d[1] << (64 - n));
    r->d[1] = (r->d[1] >> n) + (r->d[2] << (64 - n));
    r->d[2] = (r->d[2] >> n) + (r->d[3] << (64 - n));
    r->d[3] = (r->d[3] >> n);
    return ret;
}

static void secp256k1_scalar_sqr(secp256k1_scalar *r, const secp256k1_scalar *a) {
 #ifdef USE_ASM_X86_64
    uint64_t l[8];
    
    __asm__ __volatile__(
    /* Preload */
    "movq 0(%%rdi), %%r11\n"
    "movq 8(%%rdi), %%r12\n"
    "movq 16(%%rdi), %%rcx\n"
    "movq 24(%%rdi), %%r14\n"
    /* (rax,rdx) = a0 * a0 */
    "movq %%r11, %%rax\n"
    "mulq %%r11\n"
    /* Extract l0 */
    "movq %%rax, %%rbx\n" /*0(%%rsi)\n"*/
    /* (r8,r9,r10) = (rdx,0) */
    "movq %%rdx, %%r15\n"
    "xorq %%r9, %%r9\n"
    "xorq %%r10, %%r10\n"
    "xorq %%r8, %%r8\n"
    /* (r8,r9,r10) += 2 * a0 * a1 */
    "movq %%r11, %%rax\n"
    "mulq %%r12\n"
    "addq %%rax, %%r15\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    "addq %%rax, %%r15\n" /*8 rsi in r15*/
    "adcq %%rdx, %%r9\n"
    "movq %%r11, %%rax\n"
    "adcq $0, %%r10\n"
    /* Extract l1 */
   /* 8(rsi) in r15*/
    /* (r9,r10,r8) += 2 * a0 * a2 */
    "mulq %%rcx\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "movq %%r12, %%rax\n"
    "adcq $0, %%r8\n"
    /* (r9,r10,r8) += a1 * a1 */
    "mulq %%r12\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    /* Extract l2 */
    "movq %%r9, 16(%%rsi)\n"
    "movq %%r11, %%rax\n"
    "xorq %%r9, %%r9\n"
    /* (r10,r8,r9) += 2 * a0 * a3 */
    "mulq %%r14\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    "movq %%r12, %%rax\n"
    "adcq $0, %%r9\n"
    /* (r10,r8,r9) += 2 * a1 * a2 */
    "mulq %%rcx\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    "adcq $0, %%r9\n"
    "addq %%rax, %%r10\n"
    "adcq %%rdx, %%r8\n"
    "movq %%r10, %%r13\n"
    "movq %%r12, %%rax\n"
    "adcq $0, %%r9\n"
    /* Extract l3 */
    /*"movq %%r10, 24(%%rsi)\n"*/

    /* (r8,r9,r10) += 2 * a1 * a3 */
    "mulq %%r14\n"
    "xorq %%r10, %%r10\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "movq %%rcx, %%rax\n"
    "adcq $0, %%r10\n"
    /* (r8,r9,r10) += a2 * a2 */
    "mulq %%rcx\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r9\n"
    "adcq $0, %%r10\n"
    /* Extract l4 */
    /*"movq %%r8, 32(%%rsi)\n"*/
    "movq %%r8, %%r11\n"
    "movq %%rcx, %%rax\n"
    "xorq %%r8, %%r8\n"
    /* (r9,r10,r8) += 2 * a2 * a3 */
    "mulq %%r14\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "adcq $0, %%r8\n"
    "addq %%rax, %%r9\n"
    "adcq %%rdx, %%r10\n"
    "movq %%r14, %%rax\n"
    "adcq $0, %%r8\n"
    /* Extract l5 */
    /*"movq %%r9, 40(%%rsi)\n"*/
 /*   "movq %%r9, %%r12\n"*/
    /* (r10,r8) += a3 * a3 */
    "mulq %%r14\n"
    "addq %%rax, %%r10\n"
    /* Extract l6 */
    /*"movq %%r10, 48(%%rsi)\n"*/
    /*"movq %%r10, %%rcx\n"*/
    /* Extract l7 */
    /*"movq %%r8, 56(%%rsi)\n"*/
    /*"movq %%r8, %%r14\n"*/
    :
    : "S"(l), "D"(a->d)
    : "rax", "rbx", "rcx", "rdx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "cc", "memory");
        
      __asm__ __volatile__(
    /* Preload. */
  /*  "movq 32(%%rsi), %%r11\n" */
  /*  "movq 40(%%rsi), %%r9\n" */
  /*   "movq 48(%%rsi), %%r10\n" */
  /*   "movq 56(%%rsi), %%r8\n" */
  /*  "movq 0(%%rsi), %%rbx\n"  */
 /*   "movq %%rcx, %%r13\n"*/
    "movq %3, %%rax\n"
    "adcq %%rdx, %%r8\n"
    "mulq %%r11\n"
    "xor %%ecx, %%ecx\n" 
    "xorq %%r12, %%r12\n"
    "xorq %%r14, %%r14\n"
    "addq %%rax, %%rbx\n" /*q0 into rbx*/
    "adcq %%rdx, %%rcx\n"
 /*   "addq 8(%%rsi), %%rcx\n" */
    "addq %%r15, %%rcx\n" 
    "mov $0, %%r15d\n"
    "movq %3, %%rax\n"
    "adcq %%r12, %%r15\n"
    "mulq %%r9\n"
    "addq %%rax, %%rcx\n" /*q1 stored to rcx*/
    "adcq %%rdx, %%r15\n"
    "movq %4, %%rax\n" 
    "adcq %%r12, %%r14\n"
    "mulq %%r11\n"
    "addq %%rax, %%rcx\n"
    "adcq %%rdx, %%r15\n"
    "adcq %%r12, %%r14\n"
    "addq 16(%%rsi), %%r15\n"
    "adcq %%r12, %%r14\n"
    "movq %3, %%rax\n"
    "adcq %%r12, %%r12\n"
    "mulq %%r10\n"
    "movq %4, %%rsi\n"  
    "addq %%rax, %%r15\n"
    "adcq %%rdx, %%r14\n"
    "movq %%rsi, %%rax\n"
    "adcq $0, %%r12\n"
    "mulq %%r9\n"
    "addq %%rax, %%r15\n"
    "adcq %%rdx, %%r14\n"
    "adcq $0, %%r12\n"
    "movq %3, %%rax\n"
    "addq %%r11, %%r15\n" /*q2 into r15*/
    "adcq $0, %%r14\n"
    "adcq $0, %%r12\n"
    "addq %%r13, %%r14\n"
    "movq $0, %%r13\n"
    "adcq $0, %%r12\n"
    "adcq $0, %%r13\n"
    "mulq %%r8\n"
    "addq %%rax, %%r14\n"
    "movq %%rsi, %%rax\n"  
    "adcq %%rdx, %%r12\n"
    "adcq $0, %%r13\n"
    "mulq %%r10\n"
    "addq %%rax, %%r14\n"
    "adcq %%rdx, %%r12\n"
    "adcq $0, %%r13\n"
    "addq %%r14, %%r9\n" /* q3 into r9*/
    "adcq $0, %%r12\n"
    "movq %%rsi, %%rax\n" 
    "movq $0, %%r14\n"
    "adcq $0, %%r13\n"
    "mulq %%r8\n"
    "addq %%rax, %%r12\n"
    "adcq %%rdx, %%r13\n"
    "adcq %%r14, %%r14\n"
    "addq %%r12, %%r10\n" /*q4 into r10*/
    "adcq $0, %%r13\n" 
    "adcq $0, %%r14\n" 
    "addq %%r8, %%r13\n" /* q5 into r13 */ 
    "movq %3, %%rax\n"
    "movq %%rax, %%r12\n"
    "adcq $0, %%r14\n" /*q6 into r14*/
  
/* %q5 input for second operation is %q0 output from first / RBX as the connecting link
    %q6 input for second operation is %q1 output from first / RCX as the connecting link
    %q7 input for second operation is %q2 output from first / R15 as the connecting link
    %q8 input for second operation is %q3 output from first / r9 as the connecting link
    %q9  input for second operation is %q4 output from first / r10 as the connecting link*
    %q10 input for second operation is %q5 output from first / r13 as the connecting link*
    %q11 input for second operation is %q6 output from first  / r14 as the connecting link */    
    
    /* Reduce 385 bits into 258. */

    "mulq %%r10\n"
    "xorq %%r8, %%r8\n"
    "xorq %%r11, %%r11\n"
    "addq %%rax, %%rbx\n" /* q0 output*/
    "adcq %%rdx, %%r8\n"
    "addq %%rcx, %%r8\n" 
    "movq %%r12, %%rax\n"
    "adcq %%r11, %%r11\n"
    "xor %%ecx, %%ecx\n"  
    "mulq %%r13\n"
    "addq %%rax, %%r8\n"
    "adcq %%rdx, %%r11\n"
    "movq %%rsi, %%rax\n"
    "adcq %%rcx, %%rcx\n"
    "mulq %%r10\n"
    "addq %%rax, %%r8\n" /* q1 output */
    "movq %%r12, %%rax\n"
    "adcq %%rdx, %%r11\n"
    "adcq $0, %%rcx\n"
    "xorq %%r12, %%r12\n"
    "addq %%r15, %%r11\n" 
    "adcq %%r12, %%rcx\n"
    "movq %%rax, %%r15\n"
    "adcq %%r12, %%r12\n"
    "mulq %%r14\n"
    "addq %%rax, %%r11\n"
    "adcq %%rdx, %%rcx\n"
    "movq %%rsi, %%rax\n"
    "adcq $0, %%r12\n"
    "mulq %%r13\n"
    "addq %%rax, %%r11\n"
    "adcq %%rdx, %%rcx\n"
    "adcq $0, %%r12\n"
    "addq %%r10, %%r11\n" /* q2 output */
    "adcq $0, %%rcx\n"
    "adcq $0, %%r12\n"
    "addq %%r9, %%rcx\n"
    "movq %%rsi, %%rax\n"
    "adcq $0, %%r12\n"
    "mulq %%r14\n"
    "addq %%rax, %%rcx\n"
    "movq %%r15, %%rax\n"
    "adcq %%rdx, %%r12\n"
    "addq %%r13, %%rcx\n"    /* q3 output */
    "adcq $0, %%r12\n"
    "addq %%r14, %%r12\n" /* q4 output */
    
/* %q1 input for next operation is %q0 output from prior / RBX as the connecting link
    %q2 input for next operation is %q1 output from prior / r8 as the connecting link 
    %q3 input for next operation is %q2 output from prior / R11 as the connecting link  
    %q4 input for next operation is %q3 output from prior / RCX as the connecting link
    %q5 input for next operation is %q4 output from prior / r12 as the connecting link   */
        
    /* Reduce 258 bits into 256. */

    "mulq %%r12\n"   
    "addq %%rbx, %%rax\n"
    "adcq $0, %%rdx\n"
    "movq %%rax, %%r14\n"  /* 0(q2) output */
    "movq %%rdx, %%r9\n" 
    "xorq %%r10, %%r10\n"
    "addq %%r8, %%r9\n"
    "movq %%rsi, %%rax\n"
    "adcq %%r10, %%r10\n"
    "mulq %%r12\n"
    "addq %%rax, %%r9\n" /* 8(q2) output */
    "adcq %%rdx, %%r10\n" 
    "xor %%ebx, %%ebx\n"
    "addq %%r12, %%r10\n"
    "adcq %%rbx, %%rbx\n"
    "movq $0xffffffffffffffff, %%r8\n"
    "addq %%r11, %%r10\n" /* 16(q2) output */
    "movq $0, %%r11\n"
    "adcq $0, %%rbx\n"
    "addq %%rcx, %%rbx\n"  /* 24(q2) output */
    "adcq $0, %%r11\n" /* c  output */

    
    /*FINAL REDUCTION*/
    
/*    r14 carries ex 0(%%rdi), 
       r9 carries ex 8(%%rdi),
       r10 carries ex 16(%%rdi), 
       rbx carries ex 24(%%rdi)
       r11 carries c */
    "movq $0xbaaedce6af48a03b,%%r12\n"
    "movq $0xbaaedce6af48a03a,%%rcx\n"
    "movq $0xbfd25e8cd0364140,%%r13\n"
    "cmp   %%r8 ,%%rbx\n"
    "setne %%dl\n"
    "cmp   $0xfffffffffffffffd,%%r10\n"
    "setbe %%al\n"
    "or     %%eax,%%edx\n"
    "cmp  %%rcx,%%r9\n"
    "setbe %%cl\n"
    "or     %%edx,%%ecx\n"
    "cmp  %%r12,%%r9\n"
    "movzbl %%dl,%%edx\n"
    "seta  %%r12b\n"
    "cmp  %%r13,%%r14\n"
    "movzbl %%cl,%%ecx\n"
    "seta  %%r13b\n"
    "not   %%ecx\n"
    "not   %%edx\n"
    "or     %%r13d,%%r12d\n"
    "movzbl %%r12b,%%r12d\n"
    "and   %%r12d,%%ecx\n"
    "xor    %%r12d,%%r12d\n"
    "cmp   %%r8,%%r10\n"
    "sete  %%r12b\n"
    "xor   %%r13d,%%r13d\n"
    "and   %%r12d,%%edx\n"
    "or     %%edx,%%ecx\n"
    "xor   %%edx,%%edx\n"
    "add  %%ecx,%%r11d\n"
    "imulq %%r11,%%r15\n"
    "addq  %%r15,%%r14\n"
    "adcq  %%rdx,%%r13\n"  
    "imulq %%r11,%%rsi\n"
    "xorq %%r15,%%r15\n"
    "xor   %%eax,%%eax\n"
    "movq  %%r14,0(%q2)\n"
    "xor   %%edx,%%edx\n"
    "addq %%r9,%%rsi\n"
    "adcq %%rdx,%%rdx\n" 
    "addq %%rsi,%%r13\n"
    "movq %%r13,8(%q2)\n"
    "adcq %%rdx,%%r15\n"
    "addq %%r11,%%r10\n"
    "adcq %%rax,%%rax\n" 
    "addq %%r15,%%r10\n"
    "movq %%r10,16(%q2)\n"
    "adcq $0,%%rax\n"
    "addq %%rbx,%%rax\n"
    "movq %%rax,24(%q2)\n"
    : "=D"(r)
    : "S"(l), "D"(r), "n"(SECP256K1_N_C_0), "n"(SECP256K1_N_C_1)
    : "rax", "rbx", "rcx", "rdx", "r8", "r9", "r10", "r11", "r12", "r13", "r14", "r15", "cc", "memory");      
     
#else
    uint64_t l[8];
    secp256k1_scalar_sqr_512(l, a);
    secp256k1_scalar_reduce_512(r, l);
#endif    
}

#ifdef USE_ENDOMORPHISM
static void secp256k1_scalar_split_128(secp256k1_scalar *r1, secp256k1_scalar *r2, const secp256k1_scalar *a) {
    r1->d[0] = a->d[0];
    r1->d[1] = a->d[1];
    r1->d[2] = 0;
    r1->d[3] = 0;
    r2->d[0] = a->d[2];
    r2->d[1] = a->d[3];
    r2->d[2] = 0;
    r2->d[3] = 0;
}
#endif

SECP256K1_INLINE static int secp256k1_scalar_eq(const secp256k1_scalar *a, const secp256k1_scalar *b) {
    return ((a->d[0] ^ b->d[0]) | (a->d[1] ^ b->d[1]) | (a->d[2] ^ b->d[2]) | (a->d[3] ^ b->d[3])) == 0;
}

SECP256K1_INLINE static void secp256k1_scalar_mul_shift_var(secp256k1_scalar *r, const secp256k1_scalar *a, const secp256k1_scalar *b, unsigned int shift) {
    uint64_t l[8];
    unsigned int shiftlimbs;
    unsigned int shiftlow;
    unsigned int shifthigh;
    VERIFY_CHECK(shift >= 256);
    secp256k1_scalar_mul_512(l, a, b);
    shiftlimbs = shift >> 6;
    shiftlow = shift & 0x3F;
    shifthigh = 64 - shiftlow;
    r->d[0] = shift < 512 ? (l[0 + shiftlimbs] >> shiftlow | (shift < 448 && shiftlow ? (l[1 + shiftlimbs] << shifthigh) : 0)) : 0;
    r->d[1] = shift < 448 ? (l[1 + shiftlimbs] >> shiftlow | (shift < 384 && shiftlow ? (l[2 + shiftlimbs] << shifthigh) : 0)) : 0;
    r->d[2] = shift < 384 ? (l[2 + shiftlimbs] >> shiftlow | (shift < 320 && shiftlow ? (l[3 + shiftlimbs] << shifthigh) : 0)) : 0;
    r->d[3] = shift < 320 ? (l[3 + shiftlimbs] >> shiftlow) : 0;
    secp256k1_scalar_cadd_bit(r, 0, (l[(shift - 1) >> 6] >> ((shift - 1) & 0x3f)) & 1);
}

#endif
