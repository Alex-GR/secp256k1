/* Written by alexgr for a i5-8250u, based on the c file for int128 implementation 

 
 -General optimizations notes - if some things don't make sense in first sight: 
 
 1) Use of a temp array for storing/restoring registers and temp variables to maximize registers. The use of mm and xmm registers as temp store
 was explored but it's somehow way slower than it should. Context switches preservation maybe? The use of temp memory variables instead of an array was
 explored, to keep rsp in use but conserve rcx. The problem was with clang (sometimes preffers rbp for its addresses). Thus storing the rbp contents
 for later rbp-restoration, on a rbp-based address won't work and the program will crash. With gcc it worked ok though - and using rsp is generally
 better. For consistency's sake, rcx was sacrificed to keep the array address and this means +2 instruction count on saving/restoring rsp contents on
 the rcx-based array. 
 
 2) shrdq $52 for c/d >> 52 to low part become shldq $12 to high part. This preserves the low part and it can be be anded with M without needing +1 move.
 
 3) After shifting right by 52, we don't need to have a zero register for the high part. We can adcq with $0 in chained adds. Saves a few xor reg/reg.
 
 4) The use of 3-operand non-destructive ANDN (with reverse mask of M) and BZHI were explored to zero bits 52+. Somehow they didn't give optimal results. 
 
 5) In one instance, it is required to zero bits 48+ instead of 52+. Now 48 is a convenient number (2 bytes precisely). If this is stored in mem, one 
 can write it directly to the array position and then issue a movw $0 with a +6 offset to overwrite bytes 7 and 8 with zero. Apparently 2 movs are much
 cheaper than bringing in the 0xFFFFFFFFFFFF mask and anding.
 
 6) If it's not saved in ram the same can be achieved with shifting left and then back right to pull zeroes. If a copy to another register is desired, 
 three commands can be turned to two. (mov regA, regB / shl $16 / shr $16) with rorx: Rorx -16 regA, regB and then shr $16 on regB to clear the last 16
 bits with zeroes pulled in. Making a copy of the register, bringing in a mask, anding with the mask is costlier - plus moving 64 bit immediates
 consumes 10 byte opcodes - and should be avoided.
 
 7) The merging of u0 and tx can be done in less moves by pulling the 4 bits from the side (shldq/shrdq) - after properly aligning the source. 
 Incidentally the alignment can also be combined with shifthing left-right to zero the upper bits from 48+. 
 
 8) All multiplications involving a and b are from cache. Somehow it ended up faster than using registers. Sqr in particular faces big stalls
 in its execution near the middle and in the end - which seems to indicate it could benefit from having registers, but still every time I wrote it
 with registers I couldn't get it up to similar speed.
 
 9) There was some effort to keep instruction count to a minimum - sometimes at the expense of the natural flow of multiplications and additions. 
 For example if a certain variable has been loaded on rdx and it can opportunistically multiply another variable with the same rdx, then the (later) mul,
 will be pulled forward in time to save that one extra load on the rdx.
 
 10) Multiplication by 2, directly from memory, is provided by rorx'ing by a value of -1 (left by 1). It is assumed variables have zeroes in the upper
 bits - so when the zero goes to the first bit =x2. This provides a lea (reg,reg,1) effect - but from memory, when loading rdx for multiplications. A
 more correct approach is SHRX but it needs a count register as it can't be fed with an immediate like "$-1".
 
 11) The code is not as it was written: It has been extesively reordered. General instruction reordering doesn't follow any kind of "rational" logic - aside
 from not violating proper sequencing of inputs/outputs. Essentially instructions are reordered based on how they bench. If, say, we issue 4 mulxs and then 
 want to add their results, sometimes it's better to first add the later muls, or have an irregular mix of sources/destinations to get big boosts in 
 performance. Apparently the internal reordering queue of the cpu is doing weird things - and some times, if the code is positioned in a certain way
 it can get way faster (or way slower). Some times it's faster to interleave adds and mulxs sometimes it's not. Sometimes it's faster to delay operations
 and expect prior calculations to finish, sometime it's not - you can issue an addition on a multiplication where the result hasn't even been finished yet 
 (based on the time it shows when profiled with perf record ./bench_verify). In any case, things are reordered based on a specific cpu architecture results, 
 and it will almost certainly not perform the same with other cpu architectures - which will need their own reordering. */
 




#ifndef SECP256K1_FIELD_INNER5X52_IMPL_H
#define SECP256K1_FIELD_INNER5X52_IMPL_H

SECP256K1_INLINE static void secp256k1_fe_mul_inner(uint64_t *r, const uint64_t *a, const uint64_t * SECP256K1_RESTRICT b) {
uint64_t aa[9];
 

__asm__ __volatile__(

"movq %%rbp, 0(%%rcx)\n" /*storing rbp*/
"movq 32(%%rsi), %%rdx\n" /*a4*/

"movq %%rsp, 8(%%rcx)\n" /*storing rsp*/   
"movq $0xFFFFFFFFFFFFF, %%rsp\n" /*M*/
"movq %%r14, 32(%%rcx)\n" /*storing r14*/

"mulx 32(%%rbx), %%r8, %%r9\n" /*a4*b4 for c*/
"shldq $12, %%r8, %%r9\n" /*c >> 52 stored with left shift 12 on high (r9), to leave low part (r8) intact for anding */
"movq 0(%%rsi), %%rdx\n" /*a0*/
"movq %%r12, 48(%%rcx)\n" /*storing r12*/

"mulx 24(%%rbx), %%r10, %%r11\n" /*a0*b3 for d*/
"movq %%r13, 40(%%rcx)\n" /*storing r13*/
"movq 8(%%rsi), %%rdx\n" /*a1*/
"movq %%rdi, 16(%%rcx)\n" /*storing rdi*/

"mulx 16(%%rbx), %%r12, %%r13\n" /*a1*b2 for d*/
"movq %%r15, 24(%%rcx)\n" /*storing r15*/
"movq 16(%%rsi), %%rdx\n" /*a2*/

"mulx 8(%%rbx), %%r14, %%r15\n" /*a2*b1 for d*/
"movq 24(%%rsi), %%rdx\n" /*a3*/

"mulx 0(%%rbx), %%rbp, %%rdi\n" /*a3*b0 for d*/
"andq %%rsp, %%r8\n" /*c & M*/
"movq $0x1000003D10, %%rdx\n" /*R*/

"mulx %%r8, %%r8, %%rax\n" /* (c & M) * R */

"addq %%r10, %%r12\n" /*first d additions */
"adcq %%r11, %%r13\n" /*first d additions */
"addq %%r14, %%r12\n" /*first d additions */
"adcq %%r15, %%r13\n" /*first d additions */
"addq %%rbp, %%r12\n" /*first d additions */
"adcq %%rdi, %%r13\n" /*first d additions */

"mulx %%r9, %%r9, %%r10\n" /* c * R */
"movq 0(%%rsi), %%rdx\n" /*a0*/

"mulx 32(%%rbx), %%r14, %%r11\n" /*a0 *b4 for 2nd d*/
"movq 8(%%rsi), %%rdx\n" /*a1*/

"mulx 24(%%rbx), %%rbp, %%r15\n" /*a1 *b3 for 2nd d*/
"addq %%r8, %%r12\n"  /*add (c & M) * R on d */
"adcq %%rax, %%r13\n" /*add (c & M) * R on d */
"movq 16(%%rsi), %%rdx\n" /*a2*/
"shldq $12, %%r12, %%r13\n" /* d >> 52 */
"andq %%rsp, %%r12\n" /* t3 = d & M */

"mulx 16(%%rbx), %%r8, %%rdi\n" /*a2 *b2 for 2nd d*/
"movq 24(%%rsi), %%rdx\n" /*a3*/
"movq %%r12, 56(%%rcx)\n" /* t3*/
"mulx 8(%%rbx), %%rax, %%r12\n" /*a3 *b1 for 2nd d*/

"movq 32(%%rsi), %%rdx\n" /*a3*/
"addq %%r14, %%rax\n" /* 2nd d additions on first d */
"adcq %%r11, %%r12\n"/* 2nd d additions on first d */
"addq %%r8, %%r9\n" /* 2nd d additions on first d */
"adcq %%rdi, %%r10\n" /* 2nd d additions on first d */
"addq %%r13, %%rax\n" /* 2nd d additions on first d */
"adcq $0, %%r12\n"   /* 2nd d additions on first d */

"mulx 0(%%rbx), %%r14, %%r11\n" /*a4 *b0 for 2nd d*/
"addq %%rbp, %%r9\n" /* 2nd d additions on first d */
"adcq %%r15, %%r10\n"/* 2nd d additions on first d */
"addq %%r14, %%r9\n" /* 2nd d additions on first d */
"adcq %%r11, %%r10\n" /* 2nd d additions on first d */
"addq %%rax, %%r9\n" /* 2nd d additions on first d */
"adcq %%r12, %%r10\n" /* 2nd d additions on first d */

"movq 24(%%rsi), %%rdx\n" /*a3*/

"mulx 16(%%rbx), %%r8, %%rdi\n" /*a3 *b2 for 3rd d*/
"movq 32(%%rsi), %%rdx\n" /*a4*/

"mulx 8(%%rbx), %%r15, %%rbp\n" /*a4 * b1 for 3rd d */
"movq %%r9, 64(%%rcx)\n" /*t4 temp store, as is - without zeroing the last 16bits (will do next)*/
"movq 16(%%rsi), %%rdx\n" /*a2*/

"mulx 24(%%rbx), %%rax, %%r12\n" /*a2 *b3 for 3rd d*/
"movq 8(%%rsi), %%rdx\n" /*a1*/
"movw $0, 70(%%rcx)\n" /*t4 temp store, zeroing last 16 bits / cheaper than ANDing with a mask and then storing it*/

"mulx 32(%%rbx), %%r14, %%r11\n" /*a1 *b4 for 3rd d*/
"shldq $12, %%r9, %%r10\n" /* d >>52 */
"movq 0(%%rsi), %%rdx\n" /*a0*/

"addq %%r15, %%r8\n" /*3rd d additions*/
"adcq %%rbp, %%rdi\n" /*3rd d additions*/
"addq %%rax, %%r8\n" /*3rd d additions*/
"adcq %%r12, %%rdi\n" /*3rd d additions*/
"addq %%r14, %%r8\n" /*3rd d additions*/
"adcq %%r11, %%rdi\n" /*3rd d additions*/
"addq %%r10, %%r8\n" /*adding 3rd d on 2nd d*/
"adcq $0, %%rdi\n" /*adding 3rd d on 2nd d*/

"rorx $-12, %%r9, %%rbp\n" /*getting u0 and tx properly aligned for merging by moving it to the 64 bit boundary - we only need 4 bits*/

"mulx 0(%%rbx), %%rax, %%r12\n" /* c = a0*b0 */
"movq 16(%%rsi), %%rdx\n" /*a2*/
"mulx 32(%%rbx), %%r13, %%r9\n" /* a2*b4 for 4th d*/
"movq 24(%%rsi), %%rdx\n" /*a3*/
"mulx 24(%%rbx), %%r14, %%r11\n" /* a3*b3 for 4th d*/
"movq 32(%%rsi), %%rdx\n" /*a4*/

"shldq $12, %%r8, %%rdi\n" /* d >> 52*/

"addq %%r14, %%r13\n" /*4th d additions*/
"adcq %%r11, %%r9\n"  /*4th d additions*/
"addq %%rdi, %%r13\n" /*4th d additions (merging with 3rd d) */
"adcq $0, %%r9\n"     /*4th d additions (merging with 3rd d  */

"andq %%rsp, %%r8\n" /* d & M = u0 */

"mulx 16(%%rbx), %%r14, %%r11\n" /* a4*b2 for 4th d*/
"shldq $4, %%rbp, %%r8\n" /*merging u0 with tx*/
"movq $0x1000003D1, %%rdx\n" /*R>>4*/

"mulx %%r8, %%r10, %%r8\n" /* u0 * (R >> 4) */
"movq 32(%%rsi), %%rdx\n" /*a4*/

"addq %%r14, %%r13\n" /*4th d additions*/
"adcq %%r11, %%r9\n"  /*4th d additions*/
"addq %%r10, %%rax\n" /* c = c + u0 * (R>>4) */
"adcq %%r8, %%r12\n" /* c = c + u0 * (R>>4) */

"mulx 24(%%rbx), %%rbp, %%r14\n" /* a4* b3 for last d*/
"movq 0(%%rsi), %%rdx\n" /*a0*/
"mulx 8(%%rbx), %%rdi, %%r8\n" /* a0*b1 for c*/
"movq 8(%%rsi), %%rdx\n" /*a1*/
"mulx 0(%%rbx), %%r10, %%r15\n" /* a1*b0 for c*/

"shldq $12, %%r13, %%r9\n"  /* d >>52 */
"shldq $12, %%rax, %%r12\n" /* c >>52 */
"andq %%rsp, %%r13\n" /* d & M */
"movq $0x1000003D10, %%rdx\n" /*R*/
"mulx %%r13, %%r13, %%r11\n" /* d & M * R */
"movq 24(%%rsi), %%rdx\n" /*a3*/

"addq %%rdi, %%r10\n" /* c additions */
"adcq %%r8, %%r15\n"  /* c additions */
"addq %%r9, %%rbp\n" /* last d additions on d*/
"adcq  $0, %%r14\n" /* last d additions on d*/
"addq %%r12, %%r10\n" /* c additions */
"adcq $0, %%r15\n"    /* c additions */

"mulx 32(%%rbx), %%rdi, %%r9\n" /*a3 * b4 for last d*/
"movq 0(%%rbx), %%rdx\n" /*b0*/

"addq %%r13, %%r10\n" /* c + (d&M)*R */
"adcq %%r11, %%r15\n" /* c + (d&M)*R */

"mulx 16(%%rsi), %%r8, %%r13\n" /*a2 * b0 for last c*/
"movq 16(%%rbx), %%rdx\n" /*b2*/

"mulx 0(%%rsi), %%r11, %%r12\n" /*a0 * b2 for last c*/
"movq 8(%%rbx), %%rdx\n" /*b1*/

"mulx 8(%%rsi), %%rbx, %%rsi\n" /*a1 * b1 for last c*/
"addq %%rdi, %%rbp\n" /* last d additions on d*/
"movq 16(%%rcx), %%rdi\n" /*restoring rdi*/
"adcq %%r9,%%r14\n" /* last d additions on d*/
"shldq $12, %%r10, %%r15\n" /* c >>52 */
"shldq $12, %%rbp, %%r14\n" /* d >>52 */
"movq $0x1000003D10, %%rdx\n" /*R*/
"andq %%rsp, %%rax\n" /* c & M = r[0]*/
"andq %%rsp, %%rbp\n" /* d & M */

"movq %%rax, 0(%%rdi)\n" /*r[0] out*/
"andq %%rsp, %%r10\n" /*c & M = r[1] */
"mulx %%rbp, %%rax, %%rbp\n" /* d&M * R */
"addq %%r8, %%r11\n" /*last 3 c-additions on top of c*/
"adcq %%r12, %%r13\n" /*last 3 c-additions on top of c*/
"mulx %%r14, %%r8, %%r14\n" /* d * R */

"addq %%rbx, %%r11\n" /*last 3 c-additions on top of c*/
"adcq %%rsi, %%r13\n" /*last 3 c-additions on top of c*/
"addq %%r11, %%rax\n" /*last 3 c-additions on top of c*/
"adcq %%r13, %%rbp\n" /*last 3 c-additions on top of c*/
"addq %%r15, %%rax\n" /*c + d&M * R*/
"adcq $0, %%rbp\n" /*c + d&M * R*/

"movq %%r10, 8(%%rdi)\n" /*r[1] out*/
"addq 56(%%rcx), %%r8\n" /*t3 added on d * R*/
"adcq $0, %%r14\n" /*t3 added on d * R*/
"movq 48(%%rcx), %%r12\n" /*restoring r12*/
"shldq $12, %%rax, %%rbp\n" /* c >> 52 */
"andq %%rsp, %%rax\n" /* c & M = r[2]*/
"addq %%rbp, %%r8\n" /* add d * R on c */
"movq 40(%%rcx), %%r13\n" /*restoring r13*/

"movq %%rax, 16(%%rdi)\n" /* r[2] out */
"adcq $0, %%r14\n" /* add d * R on c */
"andq %%r8, %%rsp\n" /* c & M = r[3]*/
"movq 24(%%rcx), %%r15\n" /*restoring r15*/

"movq %%rsp, 24(%%rdi)\n" /* r[3] out */
"movq 0(%%rcx), %%rbp\n" /*restoring rbp*/
"movq 8(%%rcx), %%rsp\n" /*restoring rsp*/
"shrdq $52, %%r14, %%r8\n" /* c >> 52 */
"addq 64(%%rcx), %%r8\n" /* c + t4  --- t4 can be prefetced to rax which is free (way ahead of time), but, somehow, it's slower */ 

"movq %%r8, 32(%%rdi)\n" /*r[4] out */
"movq 32(%%rcx), %%r14\n" /*restoring r14*/

: "+S"(a), "+b"(b)
: "D"(r),  "c"(aa)
: "%rax", "%rdx","%r8", "%r9", "%r10", "%r11", "cc", "memory" /*rsp/rbp/rdi/r12/r13/r14/r15 used also but stored and restored from memory*/
);            
  
}


SECP256K1_INLINE static void secp256k1_fe_sqr_inner(uint64_t *r, const uint64_t *a) {

uint64_t aa[9];

__asm__ __volatile__(

"movq 32(%%rsi), %%rdx\n"  /*a4*/
"movq %%rdi, 16(%%rcx)\n" /*store rdi*/
"movq %%rsp, 0(%%rcx)\n" /*store rsp*/
"movq $0xFFFFFFFFFFFFF, %%rsp\n" /*M*/

"mulx %%rdx, %%r8, %%r9\n" /*a4 * a4 for c*/
"movq %%rbp, 8(%%rcx)\n" /*store rbp*/
"shlq $1, %%rdx\n" /*a4 * 2*/
"movq %%r12, 56(%%rcx)\n" /*store r12*/

"mulx 0(%%rsi), %%r10, %%r11\n" /*a4*2 * a0 for 2nd d*/
"movq %%r13, 48(%%rcx)\n" /*store r13*/
"rorx $-1, 0(%%rsi), %%rdx\n" /*a0 *2*/
"movq %%r15, 32(%%rcx)\n" /*store r15*/

"mulx 24(%%rsi), %%r12, %%r13\n" /*a0*2  * a3 for 1st d*/
"movq %%r14, 40(%%rcx)\n" /*store r14*/
"movq 16(%%rsi), %%rdx\n" /*a2*/
"shldq $12, %%r8, %%r9\n" /* c >> 52 on r9*/

"mulx %%rdx, %%r14, %%r15\n" /* a2 * a2 for 2nd d*/
"rorx $-1, 8(%%rsi), %%rdx\n" /*a1 *2*/
"movq %%rbx, 64(%%rcx)\n" /*store rbx*/

"mulx 16(%%rsi), %%rdi, %%rbp\n" /*a1*2  * a2 for 1st d*/
"andq %%rsp, %%r8\n" /*c*M*/

"mulx 24(%%rsi), %%rax, %%rbx\n" /*a1*2  * a3 for 2nd d*/
"addq %%r10, %%r14\n" /*2nd d additions*/
"movq $0x1000003D10, %%rdx\n" /*R*/
"adcq %%r11, %%r15\n" /*2nd d additions*/

"mulx %%r8, %%r8, %%r10\n" /*c&M * R) */
"mulx %%r9, %%r9, %%r11\n" /* c * R*/
"addq %%rdi, %%r12\n" /*1st d additions*/
"adcq %%rbp, %%r13\n" /*1st d additions*/
"addq %%rax, %%r14\n" /*2nd d additions*/
"adcq %%rbx, %%r15\n" /*2nd d additions*/
"addq %%r8, %%r12\n" /*add c&M * R on d*/
"adcq %%r10, %%r13\n" /*add c&M * R on d*/
"addq %%r9, %%r14\n" /*add c* R on 2nd d*/
"adcq %%r11, %%r15\n" /*add c* R on 2nd d*/
"shldq $12, %%r12, %%r13\n" /* d >>52*/
"rorx $-1, 16(%%rsi), %%rdx\n" /*a2 *2*/
"addq %%r13, %%r14\n" /*1st d on 2nd d*/
"adcq $0, %%r15\n" /*1st d on 2nd d*/


"mulx 24(%%rsi), %%r8, %%r10\n" /*a3 * (a2*2) for 3rd d*/
"shrdq $48, %%r14, %%r13\n" /*t4 on r13 with 16 upper bits zeroed, part #1 - we first pull the useful 48 bits */
"rorx $-1, 32(%%rsi), %%rdx\n" /*a4 *2*/
"shldq $12, %%r14, %%r15\n" /*d >> 52*/

"mulx 8(%%rsi), %%rdi, %%rbp\n" /*a1 * (a4*2) for 3rd d*/
"andq %%rsp, %%r12\n" /* d & M = t3*/
"shlq $12, %%r14\n" /*different tx as it will be merged with shldq: 3 instructions become 1 (andq/shlq/orq=>shldq) */
"shrq $16, %%r13\n" /*t4 on r13 with 16 upper bits zeroed, part #2 - then we pull zeroes in to the end*/
"addq %%rdi, %%r8\n" /*3rd d additions*/
"adcq %%rbp, %%r10\n" /*3rd d additions*/
"addq %%r15, %%r8\n" /*2nd d added on 3rd d */
"adcq $0, %%r10\n" /*2nd d added on 3rd d */

"mulx 16(%%rsi), %%r15, %%rdi\n" /* a2 *a4*2 for 4th d*/
"movq 24(%%rsi), %%rdx\n" /*a3*/

"mulx %%rdx, %%r9, %%r11\n" /* a3*a3 for 4th d*/
"shldq $12, %%r8, %%r10\n" /* d >>52*/
"shldq $12, %%r14, %%r8\n" /*merges the 4 bits of tx with u0 directly - normally this is shldq $4, but we do +8 bits extra to flood the upper bits*/
"addq %%r9, %%r15\n"  /*4th d additions on 3rd d*/
"adcq %%r11, %%rdi\n" /*4th d additions on 3rd d*/
"addq %%r10, %%r15\n" /*4th d additions on 3rd d*/
"adcq $0, %%rdi\n"    /*4th d additions on 3rd d*/
"movq 0(%%rsi), %%rdx\n" /*a0*/

"mulx %%rdx, %%rax, %%rbx\n" /* c = a0*a0*/
"shlq $1, %%rdx\n" /*a0 *2*/
"shrq $8, %%r8\n" /* ...and then zero it out by shifting it right by 8 bits which pulls zeroes in - without needing an ANDing on r8 earlier. Somehow it's faster... */

"mulx 16(%%rsi), %%r14, %%rbp\n" /*a0*2  *a2 for c*/

"mulx 8(%%rsi), %%r9, %%r11\n" /* a0*2 *a1 for c*/
"movq $0x1000003D1, %%rdx\n" /*R*/

"mulx %%r8, %%r10, %%r8\n" /*u0 * (R >> 4);*/
"addq %%r10, %%rax\n"  /* u0 * R>>4 + c*/
"adcq %%r8, %%rbx\n" /* u0 * R>>4 + c*/
"movq 8(%%rsi), %%rdx\n" /*a1*/

"mulx %%rdx, %%r8, %%r10\n" /*a1 * a1 for c*/
"rorx $-1, 32(%%rsi), %%rdx\n" /*a4 *2*/
"addq %%r14, %%r8\n"  /*adding a1*1 + a0*2 *a2 for c*/
"adcq %%rbp, %%r10\n" /*adding a1*1 + a0*2 *a2 for c*/

"mulx 24(%%rsi), %%r14, %%rsi\n" /*a3 * a4*2  */
"shldq $12, %%r15, %%rdi\n" /* d >>52*/
"movq $0x1000003D10, %%rdx\n" /*R*/
"andq %%rsp, %%r15\n" /* d & M */
"addq %%rdi, %%r14\n" /*adding a3 * a4*2 on d*/

"mulx %%r15, %%r15, %%rbp\n" /* (d & M) * R */
"adcq $0, %%rsi\n" /*adding a3 * a4*2 on d*/
"shldq $12, %%rax, %%rbx\n" /* c >> 52*/
"movq 16(%%rcx), %%rdi\n" /*restore rdi*/
"andq %%rsp, %%rax\n" /* c&M = r[0]*/
"shldq $12, %%r14, %%rsi\n" /*(d >> 52)*/

"movq %%rax, 0(%%rdi)\n" /*r[0] out*/
"addq %%r15, %%r9\n"  /* c + d&M * R */
"adcq %%rbp, %%r11\n" /* c + d&M * R */
"andq %%rsp, %%r14\n" /*d & M*/
"addq %%rbx, %%r9\n" /*  c += (a0*2) * a1; */
"adcq $0, %%r11\n"   /*  c += (a0*2) * a1; */

"mulx %%r14, %%rax, %%rbp\n" /* (d & M) * R */
"shldq $12, %%r9, %%r11\n" /* c >> 52*/ 

"mulx %%rsi, %%rsi, %%rdx\n" /* d * R */
"andq %%rsp, %%r9\n" /* c & M = r[1] */
"movq 32(%%rcx), %%r15\n" /*restore r15*/
"movq 40(%%rcx), %%r14\n" /*restore r14*/

"movq %%r9, 8(%%rdi)\n" /*r[1] out*/
"movq 64(%%rcx), %%rbx\n" /*restore rbx*/
"addq %%rax, %%r8\n" /* c additions*/
"adcq %%rbp, %%r10\n" /* c additions*/
"addq %%r11, %%r8\n" /* c additions*/
"adcq $0, %%r10\n"   /* c additions*/
"shldq $12, %%r8, %%r10\n" /* c >> 52*/
"andq %%rsp, %%r8\n" /*r[2] = c&M*/
"addq %%r12, %%rsi\n" /* c + d * R + t3*/
"movq 56(%%rcx), %%r12\n" /*restore r12*/

"movq %%r8, 16(%%rdi)\n" /*r[2] out*/
"movq 8(%%rcx), %%rbp\n" /*restore rbp*/ 
"adcq $0, %%rdx\n" /* c + d * R + t3*/
"addq %%rsi, %%r10\n" /* c + d * R + t3*/
"adcq $0, %%rdx\n" /* c + d * R + t3*/
"andq %%r10, %%rsp\n" /* r[3] = c&M& */
"shldq $12, %%r10, %%rdx\n" /* c >>52*/

"movq %%rsp, 24(%%rdi)\n" /*r[3] out*/
"addq %%r13, %%rdx\n" /*add c+t4*/

"movq %%rdx, 32(%%rdi)\n"  /*r[4] out*/
"movq 0(%%rcx), %%rsp\n"  /*restore rsp*/ 
"movq 48(%%rcx), %%r13\n" /*restore r13*/


: "+S"(a)
: "D"(r), "c"(aa)
: "%rax", "%rdx", "%r8", "%r9", "%r10", "%r11", "cc", "memory" /*rsp/rbp/rdi/rbx/r12/r13/r14/r15 used also but stored and restored from memory*/
);     
   
}

#endif

