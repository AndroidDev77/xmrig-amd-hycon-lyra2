
R"===(#ifdef cl_clang_storage_class_specifiers
#pragma OPENCL EXTENSION cl_clang_storage_class_specifiers : enable
#endif

#ifndef XEVAN_CL
#define XEVAN_CL

typedef unsigned int sph_u32;
typedef int sph_s32;
typedef unsigned long sph_u64;
typedef long sph_s64;

__constant static const sph_u64 blake2b_IV[8] = {
  0x6a09e667f3bcc908UL, 0xbb67ae8584caa73bUL,
  0x3c6ef372fe94f82bUL, 0xa54ff53a5f1d36f1UL,
  0x510e527fade682d1UL, 0x9b05688c2b3e6c1fUL,
  0x1f83d9abfb41bd6bUL, 0x5be0cd19137e2179UL
};

#define ROTL64(x,n) rotate(x,(ulong)n) 
#define ROTR64(x,n)  (((0U + (x)) << (64 - (n))) | ((x) >> (n)))
#define SWAP32(x) as_ulong(as_uint2(x).s10)
#define SWAP24(x) as_ulong(as_uchar8(x).s34567012)
#define SWAP16(x) as_ulong(as_uchar8(x).s23456701)

#define G(a,b,c,d) \
  do { \
	a += b; d ^= a; d = SWAP32(d); \
	c += d; b ^= c; b = ROTR64(b,24); \
	a += b; d ^= a; d = ROTR64(d,16); \
	c += d; b ^= c; b = ROTR64(b, 63); \
\
  } while (0);

#define G_old(a,b,c,d) \
  do { \
	a += b; d ^= a; d = ROTR64(d, 32); \
	c += d; b ^= c; b = ROTR64(b, 24); \
	a += b; d ^= a; d = ROTR64(d, 16); \
	c += d; b ^= c; b = ROTR64(b, 63); \
\
  } while (0);


#define init_sponge(state) \
do { \
	state[0]  = 0; \
	state[1]  = 0; \
	state[2]  = 0; \
	state[3]  = 0; \
	state[4]  = 0; \
	state[5]  = 0; \
	state[6]  = 0; \
	state[7]  = 0; \
	state[8]  = blake2b_IV[0]; \
	state[9]  = blake2b_IV[1]; \
	state[10] = blake2b_IV[2]; \
	state[11] = blake2b_IV[3]; \
	state[12] = blake2b_IV[4]; \
	state[13] = blake2b_IV[5]; \
	state[14] = blake2b_IV[6]; \
	state[15] = blake2b_IV[7]; \
} while(0);

/*One Round of the Blake2b's compression function*/
#define round_lyra(v)  \
    G_old(v[ 0],v[ 4],v[ 8],v[12]); \
    G_old(v[ 1],v[ 5],v[ 9],v[13]); \
    G_old(v[ 2],v[ 6],v[10],v[14]); \
    G_old(v[ 3],v[ 7],v[11],v[15]); \
    G_old(v[ 0],v[ 5],v[10],v[15]); \
    G_old(v[ 1],v[ 6],v[11],v[12]); \
    G_old(v[ 2],v[ 7],v[ 8],v[13]); \
    G_old(v[ 3],v[ 4],v[ 9],v[14]); 

#define reduceDuplexRowSetup(ptrWordIn, ptrWordInOut, ptrWordOut) \
   { \
	for (int i = 0; i < 4; i++) \
	{ \
		for (int j = 0; j < 12; j++)\
		{\
		  state[j] ^= (ptrWordIn[j+(i*12)]  + ptrWordInOut[j+(i*12)]);\
		} \
		round_lyra(state);\
		for (int j = 0; j < 12; j++)\
		{\
		   ptrWordOut[j+((3-i)*12)] = ptrWordIn[j+(i*12)]  ^ state[j];\
		} \
\
		ptrWordInOut[0+(i*12)] ^= state[11]; \
		ptrWordInOut[1+(i*12)] ^= state[0]; \
		ptrWordInOut[2+(i*12)] ^= state[1]; \
		ptrWordInOut[3+(i*12)] ^= state[2]; \
		ptrWordInOut[4+(i*12)] ^= state[3]; \
		ptrWordInOut[5+(i*12)] ^= state[4]; \
		ptrWordInOut[6+(i*12)] ^= state[5]; \
		ptrWordInOut[7+(i*12)] ^= state[6]; \
		ptrWordInOut[8+(i*12)] ^= state[7]; \
		ptrWordInOut[9+(i*12)] ^= state[8]; \
		ptrWordInOut[10+(i*12)] ^= state[9]; \
		ptrWordInOut[11+(i*12)] ^= state[10]; \
	} \
 \
   } 

#define reduceDuplexRow(ptrWordIn, ptrWordInOut, ptrWordOut) \
do { \
	 for (int i = 0; i < NCOLS; i++) \
	 { \
		 for (int j = 0; j < 12; j++)\
		 { \
			 state[j] ^= (ptrWordIn[j]  + ptrWordInOut[j]);\
		 } \
		 round_lyra(state); \
		 for (int j = 0; j < 12; j++) \
		 {\
		    ptrWordOut[j] ^= state[j];\
		 } \
\
		 ptrWordInOut[0] ^= state[11]; \
		 ptrWordInOut[1] ^= state[0]; \
		 ptrWordInOut[2] ^= state[1]; \
		 ptrWordInOut[3] ^= state[2]; \
		 ptrWordInOut[4] ^= state[3]; \
		 ptrWordInOut[5] ^= state[4]; \
		 ptrWordInOut[6] ^= state[5]; \
		 ptrWordInOut[7] ^= state[6]; \
		 ptrWordInOut[8] ^= state[7]; \
		 ptrWordInOut[9] ^= state[8]; \
		 ptrWordInOut[10] ^= state[9]; \
		 ptrWordInOut[11] ^= state[10]; \
\
		 ptrWordOut += BLOCK_LEN_INT64;\
		 ptrWordInOut += BLOCK_LEN_INT64;\
		 ptrWordIn += BLOCK_LEN_INT64;\
     }\
} while(0);

#define absorbblock(in) \
do { \
	state[0] ^= in[0]; \
	state[1] ^= in[1]; \
	state[2] ^= in[2]; \
	state[3] ^= in[3]; \
	state[4] ^= in[4]; \
	state[5] ^= in[5]; \
	state[6] ^= in[6]; \
	state[7] ^= in[7]; \
	state[8] ^= in[8];\
	state[9] ^= in[9]; \
	state[10] ^= in[10]; \
	state[11] ^= in[11]; \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
  } while(0);
  
#define absorbBlockBlake2Safe(in)  { \
	state[0] ^= in[0]; \
	state[1] ^= in[1]; \
	state[2] ^= in[2]; \
	state[3] ^= in[3]; \
	state[4] ^= in[4]; \
	state[5] ^= in[5]; \
	state[6] ^= in[6]; \
	state[7] ^= in[7]; \
	for(int j = 0; j < 12; j++){round_lyra(state);}\
  }
 
#define BLOCK_LEN_BLAKE2_SAFE_INT64 8
#define BLOCK_LEN_BLAKE2_SAFE_BYTES (BLOCK_LEN_BLAKE2_SAFE_INT64 * 8)
#define NROWS 16384
#define NCOLS 4
#define TCOST 4
#define BLOCK_LEN_INT64 12
#define BLOCK_LEN_BYTES (BLOCK_LEN_INT64 * 8)
#define LYRA2_MEMSIZE (BLOCK_LEN_INT64 * NCOLS * 8 * NROWS)
#define memMatrix(x) ( &Matrix[x * BLOCK_LEN_INT64 * NCOLS + (gid * LYRA2_MEMSIZE/8 )])

/*One Round of the Blake2b's compression function*/
#define round_lyra_prefetch(v)  \
 do { \
    G_old(v[ 0],v[ 4],v[ 8],v[12]); \
    G_old(v[ 1],v[ 5],v[ 9],v[13]); \
    G_old(v[ 2],v[ 6],v[10],v[14]); \
    G_old(v[ 3],v[ 7],v[11],v[15]); \
    G_old(v[ 0],v[ 5],v[10],v[15]); \
    G_old(v[ 1],v[ 6],v[11],v[12]); \
    G_old(v[ 2],v[ 7],v[ 8],v[13]); \
    G_old(v[ 3],v[ 4],v[ 9],v[14]); \
 } while(0);

#ifndef WORKSIZE
#define WORKSIZE 1
#endif

__attribute__((reqd_work_group_size(1, 1, 1)))
__kernel void lyra2(__global unsigned long* Matrix,__global uchar* pwd,sph_s64 pwdlen, __global sph_u64* output, ulong Threads, sph_u64 nonce){
  uint gid = get_global_id(0);
  sph_s64 kLen = 32;
  int ROW_LEN_INT64 = BLOCK_LEN_INT64 * NCOLS; // 12 * 4 = 48
  int ROW_LEN_BYTES = ROW_LEN_INT64 * 8; // 48 * 8 = 384
  int BLOCK_LEN =  BLOCK_LEN_BLAKE2_SAFE_INT64; // 8
  int nBlocksInput = ((pwdlen + 6 * BLOCK_LEN_BLAKE2_SAFE_INT64) / BLOCK_LEN_BLAKE2_SAFE_BYTES) + 1; // 2
  //============================= Basic variables ============================//
  int row = 2; //index of row to be processed
  int prev = 1; //index of prev (last row ever computed/modified)
  int rowa = 0; //index of row* (a previous row, deterministically picked during Setup and randomly picked while Wandering)
  int tau; //Time Loop iterator
  int step = 1; //Visitation step (used during Setup and Wandering phases)
  int window = 2; //Visitation window (used to define which rows can be revisited during Setup)
  int gap = 1; //Modifier to the step, assuming the values 1 or -1
  //long i; //auxiliary iteration counter
  //long v64; // 64bit var for memcpy 
  
     __global unsigned long* memPtrStart = memMatrix(0);
  //==========================================================================/
  if(gid < Threads)
  {
  //========== Initializing the Memory Matrix and pointers to it =============//
  //Tries to allocate enough space for the whole memory matrix
  //printf("Matrix Size :%i", NROWS *BLOCK_LEN_INT64  * NCOLS);
  //unsigned long Matrix[NROWS *BLOCK_LEN_INT64  * NCOLS]; // 16384 * 12 * 4 = 786,432
  for(int k = 0; k < (NROWS *BLOCK_LEN_INT64  * NCOLS); k++)
  {
     Matrix[k + (gid * LYRA2_MEMSIZE/8 )] = 0;
  }

  __global unsigned char* ptrByte = (__global unsigned char*) memPtrStart;
  //Prepends the password

  for (int j = 0; j < pwdlen; j++) {
     ptrByte[j] = pwd[j]; 
  }
  //nonce
  __global unsigned long* ptrLongNonce = (__global unsigned long*)(&ptrByte[pwdlen-4]);
  ptrLongNonce[0] = nonce + gid;
  
  
  long length = pwdlen;
  for (int j = length; j < nBlocksInput * BLOCK_LEN_BLAKE2_SAFE_BYTES - pwdlen+length; j++) { 
      ptrByte[j] = 0;
  }

  //Concatenates the basil: every integer passed as parameter, in the order they are provided by the interface
  __global unsigned long* ptrLong = (__global unsigned long*)(&ptrByte[length]);
  ptrLong[0] = kLen;
  ptrLong[1] = pwdlen; // saltlen
  ptrLong[2] = (unsigned long)0;
  ptrLong[3] = TCOST;
  ptrLong[4] = NROWS;
  ptrLong[5] = NCOLS;
  
  //Now comes the padding
  ptrByte[length+48] = 0x80;
  ptrByte[nBlocksInput * BLOCK_LEN_BLAKE2_SAFE_BYTES - 1] ^= 0x01;

  //======================= Initializing the Sponge State ====================//
  //Sponge state: 16 uint64_t, BLOCK_LEN_INT64 words of them for the bitrate (b) and the remainder for the capacity (c)
  unsigned long state[16];
  init_sponge(state);

  //set up pointers
   __global unsigned long* ptrWordRowIn;
   __global unsigned long* ptrWordRowInOut;
   __global unsigned long* ptrWordRowOut;
   
   __local unsigned long tempWordRowIn[48];
   __local unsigned long tempWordRowInOut[48];
   __local unsigned long tempWordRowOut[48];
   
  //================================ Setup Phase =============================//
  //Absorbing salt, password and basil: this is the only place in which the block length is hard-coded to 512 bits
  __global unsigned long* ptrSetUp = memPtrStart;
  for(int i = 0; i < nBlocksInput; i++){
     absorbBlockBlake2Safe(ptrSetUp);
	 ptrSetUp += BLOCK_LEN;
  }
  
  //Initializes M[0] and M[1]
    
  /// reducedSqueezeRow0
  ptrWordRowIn = (((NCOLS-1) * BLOCK_LEN_INT64) + memPtrStart);
  for (int i = 0; i < NCOLS; i++)
  {
      #pragma unroll 12
	  for (int j = 0; j<12; j++) 
	  {
		  ptrWordRowIn[j] = state[j];
	  }
	  //Goes to next block (column) that will receive the squeezed data
	  ptrWordRowIn -= BLOCK_LEN_INT64;
	  round_lyra(state);
  }

  /// reducedSqueezeRow1
  ptrWordRowIn = memPtrStart;
  ptrWordRowOut = ( (NCOLS-1) * BLOCK_LEN_INT64) + memMatrix(1);
  
  for (int i = 0; i < NCOLS; i++)
  {
      #pragma unroll 12
	  for (int j = 0; j < 12; j++) 
	  {
	    state[j] ^= ptrWordRowIn[j];
	  }
	  round_lyra(state);
	  
	  #pragma unroll 12
	  for (int j = 0; j < 12; j++) 
	  {
	   ptrWordRowOut[j] = ptrWordRowIn[j] ^ state[j];
	  }
	  ptrWordRowIn += BLOCK_LEN_INT64;
	  ptrWordRowOut -= BLOCK_LEN_INT64;
	  
  }

  //M[row] = rand; 
  //M[row*] = M[row*] XOR rotW(rand)
  do{
	 ptrWordRowIn    = memMatrix(prev);
	 ptrWordRowInOut = memMatrix(rowa);
	 ptrWordRowOut   = memMatrix(row);
     
	 for(int i = 0; i < 48; i++)
	 {
	    tempWordRowIn[i]    = ptrWordRowIn[i];
	    tempWordRowInOut[i] = ptrWordRowInOut[i];
		//tempWordRowOut[i]    = ptrWordRowOut[i];
	 }

	 reduceDuplexRowSetup(tempWordRowIn, tempWordRowInOut, tempWordRowOut);
	 for(int i = 0; i < 48; i++)
	 {
	    ptrWordRowIn[i]    = tempWordRowIn[i];
		ptrWordRowInOut[i] = tempWordRowInOut[i];
		ptrWordRowOut[i]   = tempWordRowOut[i];
	 }
	 //updates the value of row* (deterministically picked during Setup))
	 rowa = (rowa + step) & (window - 1);
	 //update prev: it now points to the last row ever computed
	 prev = row;
	 //updates row: goes to the next row to be computed
	 row++;
	  //Checks if all rows in the window where visited.
	  if (rowa == 0) {
		 step = window + gap; //changes the step: approximately doubles its value
		 window *= 2; //doubles the size of the re-visitation window
		 gap = -gap; //inverts the modifier to the step
	  }
	 
  } while(row < NROWS);
  row = 0;
  for (tau = 1; tau <= TCOST; tau++) {
		//Step is approximately half the number of all rows of the memory matrix for an odd tau; otherwise, it is -1
		step = (tau % 2 == 0) ? -1 : NROWS / 2 - 1;
		do {
			rowa = state[0] & (unsigned long)(NROWS-1);  //(USE THIS IF NROWS IS A POWER OF 2)
			
			ptrWordRowIn = memMatrix(prev);
	        ptrWordRowInOut = memMatrix(rowa);
            ptrWordRowOut  = memMatrix(row);
            
			for(int i = 0; i < 48; i++)
	        {
			   //tempWordRowIn[i]    = ptrWordRowIn[i];
			   //tempWordRowInOut[i] = ptrWordRowInOut[i];
			   //tempWordRowOut[i]    = ptrWordRowOut[i];
	        }

			reduceDuplexRow(ptrWordRowIn, ptrWordRowInOut, ptrWordRowOut);

            for(int i = 0; i < 48; i++)
	        {
	           //ptrWordRowIn[i] = tempWordRowIn[i];
		       //ptrWordRowInOut[i] = tempWordRowInOut[i];
		       //ptrWordRowOut[i] = tempWordRowOut[i];
	        }
			prev = row;
			row = (row + step) & (unsigned long)(NROWS-1); //(USE THIS IF NROWS IS A POWER OF 2)
			//int nrow = (row + step) & (unsigned int)(NROWS-1); //(USE THIS IF NROWS IS A POWER OF 2)

		} while (row != 0);
  }
  
  absorbblock( memMatrix(rowa));
  
  for(int k = 0; k < 4; k++)
  {
    output[k + (gid * 4 )] = state[k];
  }
  }



}
#endif
)==="