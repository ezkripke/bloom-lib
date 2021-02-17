#ifndef _BLOOM_H_
#define _BLOOM_H_

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SEED 42
void MurmurHash3_x64_128(const void *key, const int len, const uint32_t seed,
                         void *out);

typedef struct {
  double bpe;
  uint32_t num_hashes;
  uint32_t len_fp;
  uint8_t *fp;
} BloomFilter;

void bloom_init(BloomFilter *bf, uint64_t num_entries, double fpr) {
  bf->bpe = -log2(fpr) / log(2);
  bf->num_hashes = (uint32_t)ceil(bf->bpe * log(2));
  bf->len_fp = (uint32_t)ceil(bf->bpe * num_entries);
  bf->fp = (uint8_t *)malloc((bf->len_fp / 8 + (bf->len_fp % 8 != 0)) *
                             sizeof(uint8_t));
}

void bloom_done(BloomFilter *bf) { free(bf->fp); }

void bloom_update(BloomFilter *bf, int key) {
  for (int i = 0; i < bf->num_hashes; i++) {
    uint64_t hash_val[2];
    MurmurHash3_x64_128(&key, sizeof(int), SEED, &hash_val);
    int n = (hash_val[0] + i * hash_val[1]) % bf->len_fp;
    bf->fp[n / 8] |= 1UL << (7 - (n % 8));
  }
}

int bloom_contains(BloomFilter *bf, int key) {
  for (int i = 0; i < bf->num_hashes; i++) {
    uint64_t hash_val[2];
    MurmurHash3_x64_128(&key, sizeof(int), SEED, &hash_val);
    int n = (hash_val[0] + i * hash_val[1]) % bf->len_fp;
    if (!(bf->fp[n / 8] & 1UL << (7 - (n % 8)))) {
      return 0;
    }
  }
  return 1;
}

//-----------------------------------------------------------------------------
// The following code is taken from Peter Scott's port of MurmurHash to C.
// (https://github.com/PeterScott/murmur3)

#ifdef __GNUC__
#define FORCE_INLINE __attribute__((always_inline)) inline
#else
#define FORCE_INLINE inline
#endif

static FORCE_INLINE uint64_t rotl64(uint64_t x, int8_t r) {
  return (x << r) | (x >> (64 - r));
}

#define ROTL64(x, y) rotl64(x, y)

#define BIG_CONSTANT(x) (x##LLU)

//-----------------------------------------------------------------------------
// Block read - if your platform needs to do endian-swapping or can only
// handle aligned reads, do the conversion here

#define getblock(p, i) (p[i])

//-----------------------------------------------------------------------------
// Finalization mix - force all bits of a hash block to avalanche
static FORCE_INLINE uint64_t fmix64(uint64_t k) {
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xff51afd7ed558ccd);
  k ^= k >> 33;
  k *= BIG_CONSTANT(0xc4ceb9fe1a85ec53);
  k ^= k >> 33;

  return k;
}
//-----
//
void MurmurHash3_x64_128(const void *key, const int len, const uint32_t seed,
                         void *out) {
  const uint8_t *data = (const uint8_t *)key;
  const int nblocks = len / 16;
  int i;

  uint64_t h1 = seed;
  uint64_t h2 = seed;

  uint64_t c1 = BIG_CONSTANT(0x87c37b91114253d5);
  uint64_t c2 = BIG_CONSTANT(0x4cf5ad432745937f);

  //----------
  // body

  const uint64_t *blocks = (const uint64_t *)(data);

  for (i = 0; i < nblocks; i++) {
    uint64_t k1 = getblock(blocks, i * 2 + 0);
    uint64_t k2 = getblock(blocks, i * 2 + 1);

    k1 *= c1;
    k1 = ROTL64(k1, 31);
    k1 *= c2;
    h1 ^= k1;

    h1 = ROTL64(h1, 27);
    h1 += h2;
    h1 = h1 * 5 + 0x52dce729;

    k2 *= c2;
    k2 = ROTL64(k2, 33);
    k2 *= c1;
    h2 ^= k2;

    h2 = ROTL64(h2, 31);
    h2 += h1;
    h2 = h2 * 5 + 0x38495ab5;
  }

  //----------
  // tail

  const uint8_t *tail = (const uint8_t *)(data + nblocks * 16);

  uint64_t k1 = 0;
  uint64_t k2 = 0;

  switch (len & 15) {
  case 15:
    k2 ^= (uint64_t)(tail[14]) << 48;
  case 14:
    k2 ^= (uint64_t)(tail[13]) << 40;
  case 13:
    k2 ^= (uint64_t)(tail[12]) << 32;
  case 12:
    k2 ^= (uint64_t)(tail[11]) << 24;
  case 11:
    k2 ^= (uint64_t)(tail[10]) << 16;
  case 10:
    k2 ^= (uint64_t)(tail[9]) << 8;
  case 9:
    k2 ^= (uint64_t)(tail[8]) << 0;
    k2 *= c2;
    k2 = ROTL64(k2, 33);
    k2 *= c1;
    h2 ^= k2;

  case 8:
    k1 ^= (uint64_t)(tail[7]) << 56;
  case 7:
    k1 ^= (uint64_t)(tail[6]) << 48;
  case 6:
    k1 ^= (uint64_t)(tail[5]) << 40;
  case 5:
    k1 ^= (uint64_t)(tail[4]) << 32;
  case 4:
    k1 ^= (uint64_t)(tail[3]) << 24;
  case 3:
    k1 ^= (uint64_t)(tail[2]) << 16;
  case 2:
    k1 ^= (uint64_t)(tail[1]) << 8;
  case 1:
    k1 ^= (uint64_t)(tail[0]) << 0;
    k1 *= c1;
    k1 = ROTL64(k1, 31);
    k1 *= c2;
    h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len;
  h2 ^= len;

  h1 += h2;
  h2 += h1;

  h1 = fmix64(h1);
  h2 = fmix64(h2);

  h1 += h2;
  h2 += h1;

  ((uint64_t *)out)[0] = h1;
  ((uint64_t *)out)[1] = h2;
}

#ifdef __cplusplus
}
#endif

#endif // _BLOOM_H_
