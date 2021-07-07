#ifndef PSEUDO_RANDOM_INC
#define PSEUDO_RANDOM_INC
/*
 Pseudo-random number generator
 */

#include <stdint.h>

static uint64_t pseudo_random_seeds_[2] = { 11ULL, 1181783497276652981ULL };

static inline uint64_t xorshift128plus(uint64_t s[static 2])
{
	uint64_t x, y;
	x = s[0], y = s[1];
	s[0] = y;
	x ^= x << 23;
	s[1] = x ^ y ^ (x >> 17) ^ (y >> 26);
	y += s[1];
	return y;
}

static inline double pseudo_rand_double(void)
{
	return (double)(xorshift128plus(pseudo_random_seeds_) >> 11) * (1.0 / 9007199254740992.0);
}

static inline float pseudo_rand_float(void)
{
	return (float)(xorshift128plus(pseudo_random_seeds_) >> 11) * (1.0 / 9007199254740992.0);
}

#endif // PSEUDO_RANDOM_INC
