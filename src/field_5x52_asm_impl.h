/**********************************************************************
 * Copyright (c) 2013-2014 Diederik Huys, Pieter Wuille               *
 * Distributed under the MIT software license, see the accompanying   *
 * file COPYING or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/**
 * Changelog:
 * - March 2013, Diederik Huys:    original version
 * - November 2014, Pieter Wuille: updated to use Peter Dettman's parallel multiplication algorithm
 * - December 2014, Pieter Wuille: converted from YASM to GCC inline assembly
 */

#ifndef _SECP256K1_FIELD_INNER5X52_IMPL_H_
#define _SECP256K1_FIELD_INNER5X52_IMPL_H_

SECP256K1_INLINE static void secp256k1_fe_mul_inner(uint64_t *r, const uint64_t *a, const uint64_t * SECP256K1_RESTRICT b) {
/**
 * Registers: rdx:rax = multiplication accumulator
 *            r9:r8   = c
 *            r15:rcx = d
 *            r10-r14 = a0-a4
 *            rbx     = b
 *            rdi     = r
 *            rsi     = a / t?
 */
  uint64_t tmp1, tmp2;
__asm__ __volatile__(
    "movq 24(%%rsi),%%r13\n"
    "movq 0(%%rbx),%%rax\n"
    "movq 32(%%rsi),%%r14\n"
    /* d += a3 * b0 */
    "mulq %%r13\n"
    "movq 0(%%rsi),%%r10\n"
    "movq 8(%%rsi),%%r11\n"
    "movq %%rax,%%r9\n"
    "movq 16(%%rsi),%%r12\n"
    "movq 8(%%rbx),%%rax\n"
    "movq %%rdx,%%rsi\n"
    /* d += a2 * b1 */
    "mulq %%r12\n"
    "addq %%rax,%%r9\n"
    "movq 16(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* d += a1 * b2 */
    "mulq %%r11\n"
    "addq %%rax,%%r9\n"
    "movq 24(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* d = a0 * b3 */
    "mulq %%r10\n"
    "addq %%rax,%%r9\n"
    "movq 32(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* c = a4 * b4 */
    "mulq %%r14\n"
    "movq %%rax,%%r8\n"
    "movq %%rdx,%%rcx\n"
    /* c >>= 52 (%%r8 only) */
    "shrdq $52,%%rcx,%%r8\n"
    /* d += (c & M) * R */
    "movq $0xfffffffffffff,%%r15\n"
    "movq $0x1000003d10,%%rdx\n"
    "andq %%r15,%%rax\n"
    "mulq %%rdx\n"
    "addq %%rax,%%r9\n"
    "adcq %%rdx,%%rsi\n"
    /* t3 (tmp1) = d & M */
    "movq %%r9,%%rax\n"
    "andq %%r15,%%rax\n"
    "movq %%rax,%q1\n"  
    /* d >>= 52 */
    "movq 0(%%rbx),%%rax\n"
    "shrdq $52,%%rsi,%%r9\n"
    "xor %%esi,%%esi\n"
    /* d += a4 * b0 */
    "mulq %%r14\n"
    "addq %%rax,%%r9\n"
    "movq 8(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* d += a3 * b1 */
    "mulq %%r13\n"
    "addq %%rax,%%r9\n"
    "movq 16(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* d += a2 * b2 */
    "mulq %%r12\n"
    "addq %%rax,%%r9\n"
    "movq 24(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* d += a1 * b3 */
    "mulq %%r11\n"
    "addq %%rax,%%r9\n"
    "movq 32(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* d += a0 * b4 */
    "mulq %%r10\n"
    "addq %%rax,%%r9\n"
     /* d += c * R */
    "movq %%r8,%%rax\n"
    "movq $0x1000003d10,%%r8\n"
    "adcq %%rdx,%%rsi\n"
    "mulq %%r8\n"
    "addq %%rax,%%r9\n"
    "adcq %%rdx,%%rsi\n"
    /* t4 = d & M (%%r15) */
    "movq %%r9,%%rax\n"
    "andq %%r15,%%rax\n"
    /* d >>= 52 */
    "shrdq $52,%%rsi,%%r9\n"
    "xor %%esi,%%esi\n"
    /* tx = t4 >> 48 (tmp3) */
    "movq %%rax,%%r15\n"
    "shrq $48,%%r15\n" /*Q3*/
    
    /* t4 &= (M >> 4) (tmp2) */
    "movq $0xffffffffffff,%%rdx\n"
    "andq %%rdx,%%rax\n"
    "movq %%rax,%q2\n"
    /*"movq %q2,%%r15\n" */
    "movq 0(%%rbx),%%rax\n"
    /* c = a0 * b0 */
    "mulq %%r10\n"
    "movq %%rax,%%r8\n"
    "movq 8(%%rbx),%%rax\n"
    "movq %%rdx,%%rcx\n"
    /* d += a4 * b1 */
    "mulq %%r14\n"
    "addq %%rax,%%r9\n"
    "movq 16(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* d += a3 * b2 */
    "mulq %%r13\n"
    "addq %%rax,%%r9\n"
    "movq 24(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* d += a2 * b3 */
    "mulq %%r12\n"
    "addq %%rax,%%r9\n"
    "movq 32(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* d += a1 * b4 */
    "mulq %%r11\n"
    "addq %%rax,%%r9\n"
    "adcq %%rdx,%%rsi\n"
    
    "movq %%r15,%%rax\n"  /*Q3 transfered*/
    
    /* u0 = d & M (%%r15) */
    "movq %%r9,%%rdx\n"
    "shrdq $52,%%rsi,%%r9\n"
    "movq $0xfffffffffffff,%%r15\n"
    "xor %%esi, %%esi\n"
    "andq %%r15,%%rdx\n"
    /* d >>= 52 */

    /* u0 = (u0 << 4) | tx (%%r15) */
    "shlq $4,%%rdx\n"
    "orq %%rax,%%rdx\n"
    /* c += u0 * (R >> 4) */
    "movq $0x1000003d1,%%rax\n"
    "mulq %%rdx\n"
    "addq %%rax,%%r8\n"
    "adcq %%rdx,%%rcx\n"
    /* r[0] = c & M */
    "movq %%r8,%%rax\n"
    "andq %%r15,%%rax\n"
    "movq %%rax,0(%%rdi)\n"
    /* c >>= 52 */
    "movq 0(%%rbx),%%rax\n"
    "shrdq $52,%%rcx,%%r8\n"
    "xor %%ecx,%%ecx\n"
    /* c += a1 * b0 */
    "mulq %%r11\n"
    "addq %%rax,%%r8\n"
    "movq 8(%%rbx),%%rax\n"
    "adcq %%rdx,%%rcx\n"
    /* c += a0 * b1 */
    "mulq %%r10\n"
    "addq %%rax,%%r8\n"
    "movq 16(%%rbx),%%rax\n"
    "adcq %%rdx,%%rcx\n"
    /* d += a4 * b2 */
    "mulq %%r14\n"
    "addq %%rax,%%r9\n"
    "movq 24(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* d += a3 * b3 */
    "mulq %%r13\n"
    "addq %%rax,%%r9\n"
    "movq 32(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* d += a2 * b4 */
    "mulq %%r12\n"
    "addq %%rax,%%r9\n"
    "adcq %%rdx,%%rsi\n"
    /* c += (d & M) * R */
    "movq %%r9,%%rax\n"
    "movq $0x1000003d10,%%rdx\n"
    "andq %%r15,%%rax\n"
    "mulq %%rdx\n"
    "addq %%rax,%%r8\n"
    "adcq %%rdx,%%rcx\n"
    /* d >>= 52 */
    "shrdq $52,%%rsi,%%r9\n"
    /* r[1] = c & M */
    "movq %%r8,%%rax\n"
    "andq %%r15,%%rax\n"
    "movq %%rax,8(%%rdi)\n"
    /* c >>= 52 */
    "movq 0(%%rbx),%%rax\n"
    "shrdq $52,%%rcx,%%r8\n"
    "xor %%ecx,%%ecx\n"
    /* c += a2 * b0 */
    "mulq %%r12\n"
    "addq %%rax,%%r8\n"
    "movq 8(%%rbx),%%rax\n"
    "adcq %%rdx,%%rcx\n"
    /* c += a1 * b1 */
    "mulq %%r11\n"
    "addq %%rax,%%r8\n"
    "movq 16(%%rbx),%%rax\n"
    "adcq %%rdx,%%rcx\n"
    /* c += a0 * b2 (last use of %%r10 = a0) */
    "mulq %%r10\n"
    "addq %%rax,%%r8\n"
    /* fetch t3 (%%r10, overwrites a0), t4 (%%r15) */
    "movq 24(%%rbx),%%rax\n"
    "adcq %%rdx,%%rcx\n"
    /* d += a4 * b3 */
    "mulq %%r14\n"
    "movq %q1,%%r10\n" 
    "xor %%esi, %%esi\n"
    "addq %%rax,%%r9\n"
    "movq 32(%%rbx),%%rax\n"
    "adcq %%rdx,%%rsi\n"
    /* d += a3 * b4 */
    "mulq %%r13\n"
    "addq %%rax,%%r9\n"
    "movq $0x1000003d10,%%r11\n"
    "adcq %%rdx,%%rsi\n"
    /* c += (d & M) * R */
    "movq %%r9,%%rax\n"
    "andq %%r15,%%rax\n"
    "mulq %%r11\n"
    "addq %%rax,%%r8\n"
    "adcq %%rdx,%%rcx\n"
    /* d >>= 52 (%%r9 only) */
    "shrdq $52,%%rsi,%%r9\n"
    /* r[2] = c & M */
    "movq %%r8,%%rax\n"
    "andq %%r15,%%rax\n"
    "movq %q2,%%rsi\n"
    "movq %%rax,16(%%rdi)\n"
    /* c >>= 52 */
    "shrdq $52,%%rcx,%%r8\n"
    /* c += t3 */
    "movq %%r9,%%rax\n"
    "addq %%r10,%%r8\n"
    "xor %%ecx,%%ecx\n"
    /* c += d * R */
    "mulq %%r11\n"
    "addq %%rax,%%r8\n"
    "adcq %%rdx,%%rcx\n"
    /* r[3] = c & M */
    "movq %%r8,%%rax\n"
    "andq %%r15,%%rax\n"
    "movq %%rax,24(%%rdi)\n"
    /* c >>= 52 (%%r8 only) */
    "shrdq $52,%%rcx,%%r8\n"
    /* c += t4 (%%r8 only) */
    "addq %%rsi,%%r8\n"
    /* r[4] = c */
    "movq %%r8,32(%%rdi)\n"
: "+S"(a), "=m"(tmp1), "=m"(tmp2)
: "b"(b), "D"(r)
: "%rax", "%rcx", "%rdx", "%r8", "%r9", "%r10", "%r11", "%r12", "%r13", "%r14", "%r15", "cc", "memory");
}

SECP256K1_INLINE static void secp256k1_fe_sqr_inner(uint64_t *r, const uint64_t *a) {
/**
 * Registers: rdx:rax = multiplication accumulator
 *            r9:r8   = c
 *            rcx:rbx = d
 *            r10-r14 = a0-a4
 *            r15     = M (0xfffffffffffff)
 *            rdi     = r
 *            rsi     = a / t?
 */
  uint64_t tmp1a;
__asm__ __volatile__(
    "movq 0(%%rsi),%%r10\n"
    "movq 8(%%rsi),%%r11\n"
    "movq 16(%%rsi),%%r12\n"
    "movq 24(%%rsi),%%r13\n"
    "movq 32(%%rsi),%%r14\n"
    "leaq (%%r10,%%r10,1),%%rax\n"
    "movq $0xfffffffffffff,%%r15\n"
    /* d = (a0*2) * a3 */
    "mulq %%r13\n"
    "movq %%rax,%%rbx\n"
    "leaq (%%r11,%%r11,1),%%rax\n"
    "movq %%rdx,%%rcx\n"
    /* d += (a1*2) * a2 */
    "mulq %%r12\n"
    "addq %%rax,%%rbx\n"
    "movq %%r14,%%rax\n"
    "adcq %%rdx,%%rcx\n"
    /* c = a4 * a4 */
    "mulq %%r14\n"
    "movq %%rax,%%r8\n"
    "movq %%rdx,%%r9\n"
    /* d += (c & M) * R */
    "movq $0x1000003d10,%%rdx\n"
    "andq %%r15,%%rax\n"
    "mulq %%rdx\n"
    "addq %%rax,%%rbx\n"
    "adcq %%rdx,%%rcx\n"
    /* c >>= 52 (%%r8 only) */
    "shrdq $52,%%r9,%%r8\n"
    /* t3 (tmp1) = d & M */
    "movq %%rbx,%%rsi\n"
    "andq %%r15,%%rsi\n" /*Q1 became rsi*/
    /* d >>= 52 */
    "shrdq $52,%%rcx,%%rbx\n"
    /* a4 *= 2 */
    "movq %%r10,%%rax\n"
    "addq %%r14,%%r14\n"
    /* d += a0 * a4 */
    "mulq %%r14\n"
    "xor %%ecx,%%ecx\n"
    "addq %%rax,%%rbx\n"
    "leaq (%%r11,%%r11,1),%%rax\n"
    "adcq %%rdx,%%rcx\n"
    /* d+= (a1*2) * a3 */
    "mulq %%r13\n"
    "addq %%rax,%%rbx\n"
    "movq %%r12,%%rax\n"
    "adcq %%rdx,%%rcx\n"
    /* d += a2 * a2 */
    "mulq %%r12\n"
    "addq %%rax,%%rbx\n"

    /* d += c * R */
    "movq %%r8,%%rax\n"
    "movq $0x1000003d10,%%r8\n"
    "adcq %%rdx,%%rcx\n"
    "mulq %%r8\n"
    "addq %%rax,%%rbx\n"
    "adcq %%rdx,%%rcx\n"
    /* t4 = d & M (%%rsi) */
    "movq %%rbx,%%rdx\n"
    "andq %%r15,%%rdx\n"
    /* d >>= 52 */
    "shrdq $52,%%rcx,%%rbx\n"
    "xor %%ecx,%%ecx\n"
    /* tx = t4 >> 48 (tmp3) */
    "movq %%rdx,%%r15\n"
    "shrq $48,%%r15\n" /*Q3=R15*/
    /* t4 &= (M >> 4) (tmp2) */
    "movq $0xffffffffffff,%%rax\n"
    "andq %%rax,%%rdx\n"
    "movq %%rdx,%q1\n"/*Q2 OUT - renamed to q1*/
    /* c = a0 * a0 */
    "movq %%r10,%%rax\n"
    "mulq %%r10\n"
    "movq %%rax,%%r8\n"
    "movq %%r11,%%rax\n"
    "movq %%rdx,%%r9\n"
    /* d += a1 * a4 */
    "mulq %%r14\n"
    "addq %%rax,%%rbx\n"
    "leaq (%%r12,%%r12,1),%%rax\n"
    "adcq %%rdx,%%rcx\n"
    /* d += (a2*2) * a3 */
    "mulq %%r13\n"
    "addq %%rax,%%rbx\n"
    "adcq %%rdx,%%rcx\n"
    /* u0 = d & M (%%rsi) */
    "movq %%rbx,%%rdx\n"
    "movq $0xfffffffffffff,%%rax\n"
    "andq %%rax,%%rdx\n"
    /* d >>= 52 */
    "shrdq $52,%%rcx,%%rbx\n"
    "xor %%ecx,%%ecx\n"
    /* u0 = (u0 << 4) | tx (%%rsi) */
    "shlq $4,%%rdx\n"
    "orq %%r15,%%rdx\n" /*Q3 - R15 RETURNS*/
    /* c += u0 * (R >> 4) */
    "movq $0x1000003d1,%%rax\n"
    "movq $0xfffffffffffff,%%r15\n" /*R15 back in its place*/
    "mulq %%rdx\n"
    "addq %%rax,%%r8\n"
    "adcq %%rdx,%%r9\n"    
    /* r[0] = c & M */
    "movq %%r8,%%rax\n"
    "andq %%r15,%%rax\n"
    "movq %%rax,0(%%rdi)\n"
    /* c >>= 52 */
    "shrdq $52,%%r9,%%r8\n"
    "xorq %%r9,%%r9\n"
    /* a0 *= 2 */
    "addq %%r10,%%r10\n"
    /* c += a0 * a1 */
    "movq %%r10,%%rax\n"
    "mulq %%r11\n"
    "addq %%rax,%%r8\n"
    "movq %%r12,%%rax\n"
    "adcq %%rdx,%%r9\n"
    /* d += a2 * a4 */
    "mulq %%r14\n"
    "addq %%rax,%%rbx\n"
    "movq %%r13,%%rax\n"
    "adcq %%rdx,%%rcx\n"
    /* d += a3 * a3 */
    "mulq %%r13\n"
    "addq %%rax,%%rbx\n"
    "adcq %%rdx,%%rcx\n"
    /* c += (d & M) * R */
    "movq %%rbx,%%rax\n"
    "andq %%r15,%%rax\n"
    "movq $0x1000003d10,%%rdx\n"
    "mulq %%rdx\n"
    "addq %%rax,%%r8\n"
    "adcq %%rdx,%%r9\n"
    /* d >>= 52 */
    "shrdq $52,%%rcx,%%rbx\n"
    "xor %%ecx,%%ecx\n"
    /* r[1] = c & M */
    "movq %%r8,%%rax\n"
    "andq %%r15,%%rax\n"
    "movq %%rax,8(%%rdi)\n"
    /* c >>= 52 */
    "movq %%r10,%%rax\n"
    "shrdq $52,%%r9,%%r8\n"
    "xorq %%r9,%%r9\n"
    /* c += a0 * a2 (last use of %%r10) */
    "mulq %%r12\n"
    "addq %%rax,%%r8\n"
    "movq %%r11,%%rax\n"
    "movq %q1,%%r12\n" /*Q2 RETURNS*/
    "adcq %%rdx,%%r9\n"
    /* fetch t3 (%%r10, overwrites a0),t4 (%%rsi) */
    /*"movq %q1,%%r10\n" */
    /* c += a1 * a1 */
    "mulq %%r11\n"
    "addq %%rax,%%r8\n"
    "movq %%r13,%%rax\n"
    "adcq %%rdx,%%r9\n"
    /* d += a3 * a4 */
    "mulq %%r14\n"
    "addq %%rax,%%rbx\n"
    "adcq %%rdx,%%rcx\n"
    /* c += (d & M) * R */
    "movq %%rbx,%%rax\n"
    "andq %%r15,%%rax\n"
    "movq $0x1000003d10,%%r13\n"
    "mulq %%r13\n"
    "addq %%rax,%%r8\n"
    "adcq %%rdx,%%r9\n"
    /* d >>= 52 (%%rbx only) */
    "shrdq $52,%%rcx,%%rbx\n"
    /* r[2] = c & M */
    "movq %%r8,%%rax\n"
    "andq %%r15,%%rax\n"
    "movq %%rax,16(%%rdi)\n"
    /* c >>= 52 */
    "shrdq $52,%%r9,%%r8\n"
    "xorq %%r9,%%r9\n"
    /* c += t3 */
    "movq %%rbx,%%rax\n"
    "addq %%rsi,%%r8\n" /*RSI = Q1*/
    /* c += d * R */
    "mulq %%r13\n"
    "addq %%rax,%%r8\n"
    "adcq %%rdx,%%r9\n"
    /* r[3] = c & M */
    "movq %%r8,%%rax\n"
    "andq %%r15,%%rax\n"
    "movq %%rax,24(%%rdi)\n"
    /* c >>= 52 (%%r8 only) */
    "shrdq $52,%%r9,%%r8\n"
    /* c += t4 (%%r8 only) */
    "addq %%r12, %%r8\n"
    /* r[4] = c */
    "movq %%r8,32(%%rdi)\n"
: "+S"(a), "=m"(tmp1a)
: "D"(r)
: "%rax", "%rbx", "%rcx", "%rdx", "%r8", "%r9", "%r10", "%r11", "%r12", "%r13", "%r14", "%r15", "cc", "memory");
}

#endif
