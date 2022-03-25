/*
Copyright 2022, Yves Gallot

primegen is free source code, under the MIT license (see LICENSE). You can redistribute, use and/or modify it.
Please give feedback to the authors if improvement is realized. It is distributed in the hope that it will be useful.
*/

#include <cstdint>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>

constexpr static uint64_t add_mod(const uint64_t x, const uint64_t y, const uint64_t m)
{
	const uint64_t c = (x >= m - y) ? m : 0;
	return x + y - c;
}

const size_t M = 2 * 3 * 5 * 7;

// Create tables.
inline size_t gcd(const size_t x, const size_t y)
{
	size_t a = x, b = y;
	while (b != 0)
	{
		const size_t t = b;
		b = a % b;
		a = t;
	}
	return a;
}

inline void gen_tables()
{
	size_t r_size = 0;
	uint32_t r[M];

	std::cout << "r: ";
	for (size_t i = 1; i < M; ++i)
	{
		if (gcd(i, M) == 1)
		{
			r[r_size++] = i;
			std::cout << i << ", ";
		}
	}
	std::cout << std::endl;

	std::cout << "wheel: ";
	for (size_t i = 0; i < r_size - 1; ++i)
	{
		std::cout << r[i + 1] - r[i] << ", ";
	}
	std::cout << r[0] + M - r[r_size - 1] << std::endl;

	std::cout << "r_inv: ";
	for (size_t i = 0; i < M; ++i)
	{
		bool found = false;
		for (size_t j = 0; j < r_size; ++j)
		{
			if (i == r[j])
			{
				std::cout << j << ", ";
				found = true;
				break;
			}
		}
		if (!found) std::cout << "-1, ";
	}
	std::cout << std::endl;
}

// If packed the size of the sieve is 2 * 3 * 5 * ... * P_wheel * sieve_blk_size otherwise it is (2 - 1) * (3 - 1) * (5 - 1) * ... * (P_wheel - 1) * sieve_blk_size.
// Packed sieve is slower.
// #define PACKED_SIEVE	true

// Code is faster if sieve_blk_size is a constant expression 
template <uint32_t sieve_blk_size>
static uint64_t sieve_evalT(const uint64_t p_max)
{
	const size_t r_size = (2 - 1) * (3 - 1) * (5 - 1) * (7 - 1);
	const uint8_t r[r_size] = { 1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
								113, 121, 127, 131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181, 187, 191, 193, 197, 199, 209 };
	// wheel[i] = r[i + 1] - r[i];
	const uint8_t wheel[r_size] = { 10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4,
									2, 4, 8, 6, 4, 6, 2, 4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2 };
#ifdef PACKED_SIEVE
	// if i < r_size then r_inv[r[i]] = i else r_inv = -1
	const int8_t r_inv[M] = { -1, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, 2, -1, -1, -1, 3, -1, 4, -1, -1, -1, 5, -1, -1, -1, -1, -1, 6, -1, 7,
							  -1, -1, -1, -1, -1, 8, -1, -1, -1, 9, -1, 10, -1, -1, -1, 11, -1, -1, -1, -1, -1, 12, -1, -1, -1, -1, -1, 13, -1, 14,
							  -1, -1, -1, -1, -1, 15, -1, -1, -1, 16, -1, 17, -1, -1, -1, -1, -1, 18, -1, -1, -1, 19, -1, -1, -1, -1, -1, 20, -1, -1,
							  -1, -1, -1, -1, -1, 21, -1, -1, -1, 22, -1, 23, -1, -1, -1, 24, -1, 25, -1, -1, -1, 26, -1, -1, -1, -1, -1, -1, -1, 27,
							  -1, -1, -1, -1, -1, 28, -1, -1, -1, 29, -1, -1, -1, -1, -1, 30, -1, 31, -1, -1, -1, 32, -1, -1, -1, -1, -1, 33, -1, 34,
							  -1, -1, -1, -1, -1, 35, -1, -1, -1, -1, -1, 36, -1, -1, -1, 37, -1, 38, -1, -1, -1, 39, -1, -1, -1, -1, -1, 40, -1, 41,
							  -1, -1, -1, -1, -1, 42, -1, -1, -1, 43, -1, 44, -1, -1, -1, 45, -1, 46, -1, -1, -1, -1, -1, -1, -1, -1, -1, 47 };

#endif

#ifdef PACKED_SIEVE
	#define sieve_size(blk_size)	(blk_size * r_size)
	#define sieve_index(u, v)		(u * r_size + v)
	#define sieve_rindex(n)			((n / M) * r_size + r_inv[n % M])
#else
	#define sieve_size(blk_size)	(blk_size * M)
	#define sieve_index(u, v)		(u * M + r[v])
	#define sieve_rindex(n)			(n)
#endif

	// Create sieve bitmap
	std::vector<bool> sieve(sieve_size(sieve_blk_size), true);
	sieve[sieve_index(0, 0)] = false;

	// A vector containing: the prime p (p^2 < p_max), the current index m and the value of the wheel c = r[m] * p
	struct sp { uint32_t p; uint64_t c; size_t m; sp(const uint32_t p, const uint64_t c, const size_t m) : p(p), c(c), m(m) {} };
	std::vector<sp> v_pcm;

	// initialize v_pmc with the wheel factorized sieve of Eratosthenes
	bool end = false;
	for (size_t j = 0; j < sieve_blk_size; ++j)
	{
		for (size_t i = 0; i < r_size; ++i)
		{
			if (sieve[sieve_index(j, i)])
			{
				const uint32_t p = j * M + r[i];
				if (p * uint64_t(p) > p_max) { end = true; break; }

				size_t m = 1;
				uint64_t c = r[m] * uint64_t(p);

				v_pcm.push_back(sp(p, c, m));

				// Fill the sieve
				while (c < M * sieve_blk_size)
				{
					sieve[sieve_rindex(c)] = false;
					c += wheel[m] * uint64_t(p);
					m = add_mod(m, 1, r_size);
				}
			}
		}
		if (end) break;
	}

	std::cout << "p_max = " << double(p_max) << ", sieve_blk_size = " << sieve_blk_size << ", bitmap_size = " << sieve.size() / 8 / 1024 << " kB"
			  << ", prm_count = " << v_pcm.size() << ", prm_max = " << v_pcm.back().p << ", prm_table_size = " << v_pcm.size() * sizeof(sp) / 1024 << " kB";

	// Compute the prime-counting function as an example 
	uint64_t prm_pi = 4;	// 2, 3, 5, 7

	// Start sieving, k is segment index
	end = false;
	for (size_t k = 0; !end; ++k)
	{
		// Erase the sieve
		sieve.assign(sieve_size(sieve_blk_size), true);
		if (k == 0) sieve[sieve_index(0, 0)] = false;

		for (sp & pcm : v_pcm)
		{
			uint64_t c = pcm.c;
			size_t m = pcm.m;

			const uint32_t p = pcm.p;

			// Fill the sieve for multiples of p
			while (c < M * sieve_blk_size)
			{
				sieve[sieve_rindex(c)] = false;
				c += wheel[m] * uint64_t(p);
				m = add_mod(m, 1, r_size);
			}

			pcm.c = c - M * sieve_blk_size;
			pcm.m = m;
		}

		// Read the sieve and count primes
		for (size_t j = 0; j < sieve_blk_size; ++j)
		{
			for (size_t i = 0; i < r_size; ++i)
			{
				const uint64_t p = M * (uint64_t(sieve_blk_size) * k + j) + r[i];
				if (p <= p_max)
				{
					prm_pi += (sieve[sieve_index(j, i)]) ? 1 : 0;
				}
				else end = true;
			}
			if (end) break;
		}
	}

	return prm_pi;
}

static uint64_t sieve_eval(const uint64_t p_max)
{
	uint64_t prm_pi = 0;
	
	// We must have p_max <= (sieve_blk_size * M)^2
	const uint32_t sieve_blk_size = std::lrint(std::sqrt(double(p_max)) / M) + 1;

	if      (sieve_blk_size <= (uint32_t(1) << 13)) prm_pi = sieve_evalT<uint32_t(1) << 13>(p_max);
	else if (sieve_blk_size <= (uint32_t(1) << 14)) prm_pi = sieve_evalT<uint32_t(1) << 14>(p_max);
	else if (sieve_blk_size <= (uint32_t(1) << 15)) prm_pi = sieve_evalT<uint32_t(1) << 15>(p_max);
	else if (sieve_blk_size <= (uint32_t(1) << 16)) prm_pi = sieve_evalT<uint32_t(1) << 16>(p_max);
	else if (sieve_blk_size <= (uint32_t(1) << 17)) prm_pi = sieve_evalT<uint32_t(1) << 17>(p_max);
	else if (sieve_blk_size <= (uint32_t(1) << 18)) prm_pi = sieve_evalT<uint32_t(1) << 18>(p_max);
	else if (sieve_blk_size <= (uint32_t(1) << 19)) prm_pi = sieve_evalT<uint32_t(1) << 19>(p_max);
	else if (sieve_blk_size <= (uint32_t(1) << 20)) prm_pi = sieve_evalT<uint32_t(1) << 20>(p_max);
	else if (sieve_blk_size <= (uint32_t(1) << 21)) prm_pi = sieve_evalT<uint32_t(1) << 21>(p_max);
	else if (sieve_blk_size <= (uint32_t(1) << 22)) prm_pi = sieve_evalT<uint32_t(1) << 22>(p_max);
	else if (sieve_blk_size <= (uint32_t(1) << 23)) prm_pi = sieve_evalT<uint32_t(1) << 23>(p_max);
	else std::cerr << "p_max = " << double(p_max) << " is too large." << std::endl;

	return prm_pi;
}

int main()
{
	std::cerr << "primegen: fast prime number list generator" << std::endl;
	std::cerr << " Copyright (c) 2022, Yves Gallot" << std::endl;
	std::cerr << " primegen is free source code, under the MIT license." << std::endl << std::endl;

	// Default value is M = 2 * 3 * 5 * 7 but gen_tables can be used to compute another set of parameters.
	// gen_tables(); return EXIT_SUCCESS;

	// Iterate p_max to check speed
	const uint64_t e9 = uint64_t(1000000000);
	for (uint64_t p_max = e9; p_max <= 10 * e9 * e9; p_max *= 10)
	{
		const auto t0 = std::chrono::steady_clock::now();

		// Compute the prime-counting function as an example 
		const uint64_t prm_pi = sieve_eval(p_max);

		if (prm_pi != 0)
		{
			const double duration = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
			std::cout << ", number of primes: " << prm_pi << ", " << round(duration * 100) / 100 << " sec." << std::endl;
		}
		else break;
	}

	return EXIT_SUCCESS;
}
