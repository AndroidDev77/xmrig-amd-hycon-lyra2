
R"===(
#pragma OPENCL EXTENSION cl_clang_storage_class_specifiers : enable

//typedef uint sph_u32;

typedef unsigned int sph_u32;
typedef int sph_s32;
#ifndef __OPENCL_VERSION__
typedef unsigned long sph_u64;
typedef long  sph_s64;
#else
typedef unsigned long sph_u64;
typedef long sph_s64;
#endif
/*Blake2b IV Array*/
__constant static const sph_u64 blake2b_IV[8] =
{
  0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL,
  0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
  0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL,
  0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL
};


#define ROTL64(x,n) rotate(x,(ulong)n)
#define ROTR64(x,n) rotate(x,(ulong)(64-n))
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
  } while (0)

#define G_old(a,b,c,d) \
  do { \
	a += b; d ^= a; d = ROTR64(d, 32); \
	c += d; b ^= c; b = ROTR64(b, 24); \
	a += b; d ^= a; d = ROTR64(d, 16); \
	c += d; b ^= c; b = ROTR64(b, 63); \
\
  } while (0)


#define init_sponge(state)  { \
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
} while(0)

/*One Round of the Blake2b's compression function*/
#define round_lyra(v)  \
 do { \
    G(v[ 0],v[ 4],v[ 8],v[12]); \
    G(v[ 1],v[ 5],v[ 9],v[13]); \
    G(v[ 2],v[ 6],v[10],v[14]); \
    G(v[ 3],v[ 7],v[11],v[15]); \
    G(v[ 0],v[ 5],v[10],v[15]); \
    G(v[ 1],v[ 6],v[11],v[12]); \
    G(v[ 2],v[ 7],v[ 8],v[13]); \
    G(v[ 3],v[ 4],v[ 9],v[14]); \
 } while(0)


#define reduceDuplexRowSetup(rowIn, rowInOut, rowOut) \
   { \
	for (int i = 0; i < 8; i++) \
				{ \
\
		for (int j = 0; j < 12; j++) {state[j] ^= as_ulong(Matrix[12 * i + j][rowIn]) + as_ulong(Matrix[12 * i + j][rowInOut]);} \
		round_lyra(state); \
		for (int j = 0; j < 12; j++) {Matrix[j + 84 - 12 * i][rowOut] = Matrix[12 * i + j][rowIn] ^ state[j];} \
\
		Matrix[0 + 12 * i][rowInOut] ^= state[11]; \
		Matrix[1 + 12 * i][rowInOut] ^= state[0]; \
		Matrix[2 + 12 * i][rowInOut] ^= state[1]; \
		Matrix[3 + 12 * i][rowInOut] ^= state[2]; \
		Matrix[4 + 12 * i][rowInOut] ^= state[3]; \
		Matrix[5 + 12 * i][rowInOut] ^= state[4]; \
		Matrix[6 + 12 * i][rowInOut] ^= state[5]; \
		Matrix[7 + 12 * i][rowInOut] ^= state[6]; \
		Matrix[8 + 12 * i][rowInOut] ^= state[7]; \
		Matrix[9 + 12 * i][rowInOut] ^= state[8]; \
		Matrix[10 + 12 * i][rowInOut] ^= state[9]; \
		Matrix[11 + 12 * i][rowInOut] ^= state[10]; \
				} \
 \
   } 

#define reduceDuplexRow(rowIn, rowInOut, rowOut) \
do { \
	 for (int i = 0; i < 4; i++) \
	 { \
		 for (int j = 0; j < 12; j++){ \
			 state[j] ^= Matrix[12 * i + j][rowIn] + Matrix[12 * i + j][rowInOut];} \
\
		 round_lyra(state); \
		 for (int j = 0; j < 12; j++) {Matrix[j + 12 * i][rowOut] ^= state[j];} \
\
		 Matrix[0 + 12 * i][rowInOut] ^= state[11]; \
		 Matrix[1 + 12 * i][rowInOut] ^= state[0]; \
		 Matrix[2 + 12 * i][rowInOut] ^= state[1]; \
		 Matrix[3 + 12 * i][rowInOut] ^= state[2]; \
		 Matrix[4 + 12 * i][rowInOut] ^= state[3]; \
		 Matrix[5 + 12 * i][rowInOut] ^= state[4]; \
		 Matrix[6 + 12 * i][rowInOut] ^= state[5]; \
		 Matrix[7 + 12 * i][rowInOut] ^= state[6]; \
		 Matrix[8 + 12 * i][rowInOut] ^= state[7]; \
		 Matrix[9 + 12 * i][rowInOut] ^= state[8]; \
		 Matrix[10 + 12 * i][rowInOut] ^= state[9]; \
		 Matrix[11 + 12 * i][rowInOut] ^= state[10]; \
     }\
} while(0)

#define absorbblock(in) \
do { \
	state[0] ^= Matrix[0][in]; \
	state[1] ^= Matrix[1][in]; \
	state[2] ^= Matrix[2][in]; \
	state[3] ^= Matrix[3][in]; \
	state[4] ^= Matrix[4][in]; \
	state[5] ^= Matrix[5][in]; \
	state[6] ^= Matrix[6][in]; \
	state[7] ^= Matrix[7][in]; \
	state[8] ^= Matrix[8][in]; \
	state[9] ^= Matrix[9][in]; \
	state[10] ^= Matrix[10][in]; \
	state[11] ^= Matrix[11][in]; \
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
  } while(0)
  
#define absorbBlockBlake2Safe(in)  { \
	state[0] ^= Matrix[0][in]; \
	state[1] ^= Matrix[1][in]; \
	state[2] ^= Matrix[2][in]; \
	state[3] ^= Matrix[3][in]; \
	state[4] ^= Matrix[4][in]; \
	state[5] ^= Matrix[5][in]; \
	state[6] ^= Matrix[6][in]; \
	state[7] ^= Matrix[7][in]; \
	round_lyra(state);\
  }
//lyra2 algo 

#define BLOCK_LEN_BLAKE2_SAFE_INT64 8                           
#define BLOCK_LEN_BLAKE2_SAFE_BYTES (BLOCK_LEN_BLAKE2_SAFE_INT64 * 8)
#define NROWS 16384
#define NCOLS 4
#define TCOST 4
#define BLOCK_LEN_INT64 12
#define BLOCK_LEN_BYTES (BLOCK_LEN_INT64 * 8)
#define LYRA2_MEMSIZE (BLOCK_LEN_INT64 * NCOLS * 8 * NROWS)

__attribute__((reqd_work_group_size(8, 1, 1)))
__kernel void lyra2(__global uchar* pwd, sph_u64 pwdlen, __global sph_u64* output )
{
  uint gid = get_global_id(0);
  sph_u64 kLen = 32;
  //__global hash_t *hash = &(hashes[gid-get_global_offset(0)]);
  long ROW_LEN_INT64 = BLOCK_LEN_INT64 * NCOLS;
  long ROW_LEN_BYTES = ROW_LEN_INT64 * 8;
  long BLOCK_LEN =  BLOCK_LEN_BLAKE2_SAFE_INT64;
  long nBlocksInput = ((pwdlen + 6 * 8) / BLOCK_LEN_BLAKE2_SAFE_BYTES) + 1;
  //============================= Basic variables ============================//
  long row = 2; //index of row to be processed
  long prev = 1; //index of prev (last row ever computed/modified)
  long rowa = 0; //index of row* (a previous row, deterministically picked during Setup and randomly picked while Wandering)
  long tau; //Time Loop iterator
  long step = 1; //Visitation step (used during Setup and Wandering phases)
  long window = 2; //Visitation window (used to define which rows can be revisited during Setup)
  long gap = 1; //Modifier to the step, assuming the values 1 or -1
  long i; //auxiliary iteration counter
  long v64; // 64bit var for memcpy
  //==========================================================================/

  ulong Matrix[16384][4]; // very uncool
  
  __global ulong* password = (__global ulong*)pwd;
  long  zeroPadding = nBlocksInput * BLOCK_LEN_BLAKE2_SAFE_BYTES - pwdlen;
  
  //__global uint2 Matrix[16384][4]; // very uncool
  
  int matIdx = 0;
  for (int i = 0; i < pwdlen/8; i++) { 
     Matrix[matIdx][0] = (*password+i);
	 matIdx = matIdx + 1; //password;
  }
  for (int i = 0; i < zeroPadding; i++) { 
     Matrix[matIdx][0] = 0;
	 matIdx = matIdx + 1;
	 } //0s
  
  //Concatenates the basil: every integer passed as parameter, in the order they are provided by the interface

  Matrix[matIdx][0] = kLen;
  matIdx += 1;
  Matrix[matIdx][0] = pwdlen;
  matIdx += 1;
  Matrix[matIdx][0] = 0;
  matIdx += 1;;
  Matrix[matIdx][0] = TCOST;
  matIdx += 1;
  Matrix[matIdx][0] = NROWS;
  matIdx += 1;
  Matrix[matIdx][0] = NCOLS;
  matIdx += 1;
  
  //Now comes the padding
  unsigned char* paddingByte = (unsigned char *)(&Matrix[matIdx][0]);
  *paddingByte = 0x80;
  unsigned char* paddingByte2 = ((unsigned char *)(&Matrix[0][0])) + (nBlocksInput * BLOCK_LEN_BLAKE2_SAFE_BYTES - 1);
  *paddingByte2 ^= 0x01;

  // Init State
  ulong state[16];
  init_sponge(state);
  
  
  int blockIndex = 0;
  
  //================================ Setup Phase =============================//
  //Absorbing salt, password and basil: this is the only place in which the block length is hard-coded to 512 bits
  for(int i = 0; i < nBlocksInput; i++){
     absorbBlockBlake2Safe(blockIndex);
	 blockIndex += 1;
  }
  
  //Initializes M[0] and M[1]
  
  //blake2blyra x2 

  for (int i = 0; i < 24; i++) {round_lyra(state);} //because 12 is not enough

  
  /// reducedSqueezeRow0

  for (int i = 0; i < 4; i++)
  {
	  for (int j = 0; j<12; j++) {Matrix[j + 84 - 12 * i][0] = state[j];}
	  round_lyra(state);
  }

  /// reducedSqueezeRow1

  for (int i = 0; i < 4; i++)
  {
	  for (int j = 0; j < 12; j++) {state[j] ^= Matrix[j + 12 * i][0];}
	  round_lyra(state);
	  for (int j = 0; j < 12; j++) {Matrix[j + 84 - 12 * i][1] = Matrix[j + 12 * i][0] ^ state[j];}
  }
  //M[row] = rand; 
  //M[row*] = M[row*] XOR rotW(rand)
  
  do{
         reduceDuplexRowSetup(prev, rowa, row);
		 
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
  
  
  //reduceDuplexRowSetup(1, 0, 2);
  //reduceDuplexRowSetup(2, 1, 3);
  //reduceDuplexRowSetup(3, 0, 4);
  //reduceDuplexRowSetup(4, 3, 5);
  //reduceDuplexRowSetup(5, 2, 6);
  //reduceDuplexRowSetup(6, 1, 7);

  /*sph_u32 rowa;
  rowa = state[0].x & 7;

  reduceDuplexRow(7, rowa, 0);
  rowa = state[0].x & 7;
  reduceDuplexRow(0, rowa, 3);
  rowa = state[0].x & 7;
  reduceDuplexRow(3, rowa, 6);
  rowa = state[0].x & 7;
  reduceDuplexRow(6, rowa, 1);
  rowa = state[0].x & 7;
  reduceDuplexRow(1, rowa, 4);
  rowa = state[0].x & 7;
  reduceDuplexRow(4, rowa, 7);
  rowa = state[0].x & 7;
  reduceDuplexRow(7, rowa, 2);
  rowa = state[0].x & 7;
  reduceDuplexRow(2, rowa, 5);*/
  row = 0;
  for (tau = 1; tau <= TCOST; tau++) {
		//Step is approximately half the number of all rows of the memory matrix for an odd tau; otherwise, it is -1
		step = (tau % 2 == 0) ? -1 : NROWS / 2 - 1;
		do {
			//Selects a pseudorandom index row*
			//------------------------------------------------------------------------------------------
			rowa = state[0] & (uint)(NROWS-1);  //(USE THIS IF NROWS IS A POWER OF 2)
			//rowa = state[0] % NROWS; //(USE THIS FOR THE "GENERIC" CASE)
			//------------------------------------------------------------------------------------------
			//__builtin_prefetch((ulong*)(memMatrix(rowa))+0);
			//__builtin_prefetch((ulong*)(memMatrix(rowa))+4);
			//__builtin_prefetch((ulong*)(memMatrix(rowa))+8);

			//Performs a reduced-round duplexing operation over M[row*] XOR M[prev], updating both M[row*] and M[row]
			reduceDuplexRow(prev, rowa, row);

			//update prev: it now points to the last row ever computed
			prev = row;

			//updates row: goes to the next row to be computed
			//------------------------------------------------------------------------------------------
			row = (row + step) & (uint)(NROWS-1); //(USE THIS IF NROWS IS A POWER OF 2)
			//row = (row + step) % NROWS; //(USE THIS FOR THE "GENERIC" CASE)
			//------------------------------------------------------------------------------------------
			int nrow = (row + step) & (unsigned int)(NROWS-1); //(USE THIS IF NROWS IS A POWER OF 2)
			//__builtin_prefetch((uint64_t*)(memMatrix(nrow))+0);
			//__builtin_prefetch((uint64_t*)(memMatrix(nrow))+4);
			//__builtin_prefetch((uint64_t*)(memMatrix(nrow))+8);

		} while (row != 0);
	}

  absorbblock(rowa);
  
  for (int i = 0; i < 4; i++)
  {
     output[i] = state[i];
  } 
  barrier(CLK_LOCAL_MEM_FENCE);

}
}
)==="