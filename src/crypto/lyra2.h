/**
 * Header file for the Lyra2 Password Hashing Scheme (PHS).
 *
 * Author: The Lyra PHC team (http://www.lyra-kdf.net/) -- 2014.
 *
 * This software is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef LYRA2_H_
#define LYRA2_H_
#   define __restrict__ __restrict
#include <stdint.h>
#include "Lyra.h"
#include "common/xmrig.h"

typedef unsigned char byte;

//struct LYRA2_ctx {
//	uint64_t* wholeMatrix;
//};

//Block length required so Blake2's Initialization Vector (IV) is not overwritten (THIS SHOULD NOT BE MODIFIED)
#define BLOCK_LEN_BLAKE2_SAFE_INT64 8                                   //512 bits (=64 bytes, =8 uint64_t)
#define BLOCK_LEN_BLAKE2_SAFE_BYTES (BLOCK_LEN_BLAKE2_SAFE_INT64 * 8)   //same as above, in bytes


#ifdef BLOCK_LEN_BITS
#define BLOCK_LEN_INT64 (BLOCK_LEN_BITS/64)      //Block length: 768 bits (=96 bytes, =12 uint64_t)
#define BLOCK_LEN_BYTES (BLOCK_LEN_BITS/8)       //Block length, in bytes
#else   //default block lenght: 768 bits
#define BLOCK_LEN_INT64 12                       //Block length: 768 bits (=96 bytes, =12 uint64_t)
#define BLOCK_LEN_BYTES (BLOCK_LEN_INT64 * 8)    //Block length, in bytes
#endif

#define NROWS 16384
#define NCOLS 4
#define TCOST 4
#define LYRA2_MEMSIZE (BLOCK_LEN_INT64 * NCOLS * 8 * NROWS)

#define memMatrix(x)  (&ctx->memory[x * BLOCK_LEN_INT64 * NCOLS])


int LYRA2(void* ctx2, void* K, int64_t kLen, const void* pwd, int32_t pwdlen);
void* LYRA2_create(void);
void LYRA2_destroy(void* c);


template<xmrig::Algo ALGO, bool SOFT_AES, xmrig::Variant VARIANT>
inline void lyra2_hash(const uint8_t* __restrict__ input, size_t size, uint8_t* __restrict__ output, Lyra2_ctx** __restrict__ ctx)
{
	LYRA2(ctx, output, 32, input, size);
}
template<xmrig::Algo ALGO, bool SOFT_AES, xmrig::Variant VARIANT>
static inline void lyra2_hash(const uint8_t* input, size_t size, uint8_t* output, void* ctx)
{
	LYRA2(ctx, output, 32, input, size);
}

#endif /* LYRA2_H_ */
