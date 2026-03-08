#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

typedef unsigned __int128 uint128_t;

bool quiet;
bool verbose;
bool skip_zeros;

//
// Fast Modular Arithmetic
//

uint64_t 
add_mod (
    uint64_t a, 
    uint64_t b, 
    uint64_t p
    ) 

{
    uint64_t r = a + b;
    return (r >= p) ? r - p : r;
}

uint64_t 
sub_mod (
    uint64_t a, 
    uint64_t b, 
    uint64_t p
    ) 

{
    return (a >= b) ? a - b : a + p - b;
}

uint64_t 
mul_mod (
    uint64_t a, 
    uint64_t b, 
    uint64_t p
    ) 

{
    return (uint64_t)(((uint128_t) a * b) % p);
}

uint64_t 
pow_mod (
    uint64_t b, 
    uint64_t e, 
    uint64_t p
    ) 

{
    uint64_t r = 1;
    b %= p;
    while (e > 0) {
        if (e % 2 == 1) r = mul_mod(r, b, p);
        b = mul_mod(b, b, p);
        e /= 2;
    }
    return r;
}

uint64_t 
inv_mod (
    uint64_t a, 
    uint64_t p
    ) 

{
    return pow_mod(a, p - 2, p);
}

//
// Primality Testing (Miller-Rabin)
//

bool 
miller_rabin_test (
    uint64_t d, 
    uint64_t p, 
    uint64_t a
    ) 

{
    uint64_t x = pow_mod(a, d, p);
    if (x == 1 || x == p - 1) return true;
    while (d != p - 1) {
        x = mul_mod(x, x, p);
        d *= 2;
        if (x == 1) return false; 
        if (x == p - 1) return true;
    }
    return false;
}

bool 
is_prime (
    uint64_t p
    ) 

{
    if (p <= 1 || p == 4) return false;
    if (p <= 3) return true;

    uint64_t d = p - 1;
    while (d % 2 == 0) d /= 2;
    
    static const uint64_t bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int i = 0; i < 12; i++) {
        if (p <= bases[i]) break;
        if (!miller_rabin_test(d, p, bases[i])) return false;
    }
    return true;
}

//
// Number Theory Helpers
//

// Find k-th root of 1 modulo p. In other words, find element of order k modulo p.
uint64_t 
find_primitive_root (
    uint64_t k, 
    uint64_t p
    ) 

{
    uint64_t factors[64];
    int num_factors = 0;
    uint64_t t = k;
    
    for (uint64_t i = 2; i * i <= t; i++) {
        if (t % i == 0) {
            factors[num_factors++] = i;
            while (t % i == 0) t /= i;
        }
    }
    if (t > 1) factors[num_factors++] = t;
    
    uint64_t exponent = (p - 1) / k;
    for (uint64_t x = 2; x < p; x++) {
        uint64_t z = pow_mod(x, exponent, p);
        if (z == 1) continue;
        
        bool is_primitive = true;
        for (int i = 0; i < num_factors; i++) {
            if (pow_mod(z, k / factors[i], p) == 1) {
                is_primitive = false;
                break;
            }
        }
        if (is_primitive) return z;
    }
    return 0; 
}

//
// The O(K log N) 1D Kasteleyn Engine
//

// Returns Lucas number U(n+1) mod p, where U(0) = 0, U(1) = 1, and U(n) = a * U(n-1) + U(n-2) for n > 1.
uint64_t 
lucas_number (
    uint64_t n, 
    uint64_t a,
    uint64_t p
    ) 

{
    if (n == 0) return 1;

    // Initialize U(k) and U(k+1) for k = 1
    uint64_t uk = 1;
    uint64_t uk1 = a;

    for (int i = 62 - __builtin_clzll(n); i >= 0; i--) {
        // Doubling step: compute U(2k) and U(2k+1)
        uint64_t u2k  = ((uint128_t) uk * (((uint128_t) (2 * uk1) + ((uint128_t) (p - a) * uk)) % p)) % p;
        uint64_t u2k1 = ((uint128_t) uk * uk + (uint128_t) uk1 * uk1) % p;

        // Jump from k to either 2k or 2k+1 depending on the current bit
        if ((n >> i) & 1) {
            uint64_t u2k2 = ((uint128_t) a * u2k1 + u2k) % p;
            uk  = u2k1;
            uk1 = u2k2;
        }
        else {
            uk  = u2k;
            uk1 = u2k1;
        }
    }

    return uk1;
}

uint64_t 
evaluate_kasteleyn_mod_p (
    uint64_t n, 
    uint64_t m, 
    uint64_t p
    ) 

{
    // Impossible to perfectly tile an odd-by-odd grid.
    if (n % 2 == 1 && m % 2 == 1) return 0;

    // Isolate dimensions for the 1D optimization
    uint64_t k = (n < m) ? n : m;
    uint64_t n_ = (n < m) ? m : n;

    // We need the 2(K+1) root of unity to generate 2a_j
    uint64_t t = 2 * (k + 1);
    uint64_t w = find_primitive_root(t, p);
    
    uint64_t inv_w = inv_mod(w, p);
    uint64_t v1 = add_mod(w, inv_w, p); // v1 = w + w^-1
    uint64_t v_prev = 2;                // v0 = 2
    uint64_t v_curr = v1;
    
    uint64_t r = 1;
    uint64_t k_half = k / 2;
    
    for (uint64_t j = 1; j <= k_half; j++) {

        uint64_t factor = lucas_number(n_, v_curr, p);
        r = mul_mod(r, factor, p);

        // Advance the Chebyshev sequence to prep for the next loop (V_j = 2 cos(j \theta))
        uint64_t v_next = sub_mod(mul_mod(v1, v_curr, p), v_prev, p);
        v_prev = v_curr;
        v_curr = v_next;
    }

    return r;
}

//
// Divide-and-Conquer CRT using GMP
//

void 
crt_combine (
    mpz_t r0, 
    mpz_t p0, 
    mpz_t r1, 
    mpz_t p1
    ) 

{
    mpz_t inv, t, P;
    mpz_inits(inv, t, P, NULL);
    
    mpz_invert(inv, p0, p1); 
    mpz_sub(t, r1, r0);
    mpz_mul(t, t, inv);
    mpz_mod(t, t, p1); 
    
    mpz_mul(P, p0, p1);
    mpz_mul(t, t, p0);
    mpz_add(r0, r0, t);
    mpz_mod(r0, r0, P); 
    
    mpz_set(p0, P);
    mpz_clears(inv, t, P, NULL);
}

void 
crt_tree (
    mpz_t* r, 
    mpz_t* p, 
    int start, 
    int end
    ) 

{
    if (start == end) return;
    if (start + 1 == end) {
        crt_combine(r[start], p[start], r[end], p[end]);
        return;
    }
    int mid = start + (end - start) / 2;
    crt_tree(r, p, start, mid);
    crt_tree(r, p, mid + 1, end);
    crt_combine(r[start], p[start], r[mid + 1], p[mid + 1]);
}

//
// O(K) Target Bits Calculation (Prevents Float Overflow)
//

double 
calculate_target_bits (
    uint64_t n, 
    uint64_t m
    ) 

{
    double exact_bits = 0.0;
    uint64_t k = (n < m) ? n : m;
    uint64_t n_ = (n < m) ? m : n;
    uint64_t k_half = k / 2;
    
    for (uint64_t j = 1; j <= k_half; j++) {
        double a_j = cos(M_PI * j / (k + 1.0));
        double c_j = a_j + sqrt(1.0 + a_j * a_j);
        double c_bar_j = a_j - sqrt(1.0 + a_j * a_j);
        
        if (n_ < 500) {
            // Safe to compute natively
            double u = (pow(c_j, n_ + 1) - pow(c_bar_j, n_ + 1)) / (c_j - c_bar_j);
            exact_bits += log2(u);
        } 
        else {
            // For large N, bar_c_j approaches 0. Using logarithmic arithmetic 
            // to bypass the double-precision boundary of 10^308.
            exact_bits += (n_ + 1) * log2(c_j) - log2(c_j - c_bar_j);
        }
    }
    
    return exact_bits + 64.0;
}

//
// Main
//

int 
main (
    int argc, 
    char** argv
    ) 

{
    uint64_t m_min = 0;
    uint64_t m_max = 0;
    uint64_t n = 0;
    char* output_file_name = NULL;

    for (int i = 1; i < argc; i++) {
        char* arg = argv[i];
        if (arg[0] == '-') {
            if (i < argc - 1 && strcmp(arg, "--min") == 0) {
                m_min = strtoull(argv[++i], NULL, 10);
                continue;
            }
            else if (i < argc - 1 && (strcmp(arg, "-m") == 0 || strcmp(arg, "--max") == 0)) {
                m_max = strtoull(argv[++i], NULL, 10);
                continue;
            }
            else if (i < argc - 1 && (strcmp(arg, "-o") == 0 || strcmp(arg, "--output") == 0)) {
                output_file_name = argv[++i];
                continue;
            }
            else if (strcmp(arg, "-q") == 0 || strcmp(arg, "--quiet") == 0) {
                quiet = 1;
                continue;
            }
            else if (strcmp(arg, "-v") == 0 || strcmp(arg, "--verbose") == 0) {
                verbose = 1;
                continue;
            }
            else if (strcmp(arg, "--skip-zeros") == 0) {
                skip_zeros = 1;
                continue;
            }
            printf("ERROR: Unexpected option '%s'.\n", arg);
            return 1;
        }
        else if (n == 0) {
            n = strtoull(arg, NULL, 10);
            continue;
        }
        else if (m_min == 0) {
            m_max = m_min = strtoull(arg, NULL, 10);
            continue;
        }
    }

    if (n == 0) {
        printf("Usage: %s <N> <M> [--max <M_max>] [--verbose] [--skip-zeros] [--output <file name>]\n", argv[0]);
    }

    if (m_min == 0) m_min = 1;
    if (m_max == 0) m_max = 1000;

    FILE* output = NULL;
    if (output_file_name != NULL) {
        output = fopen(output_file_name, "w");
    }

    double total_start_time;
    if (verbose) total_start_time = omp_get_wtime();

    int max_threads = omp_get_max_threads();
    int max_primes = (int)(calculate_target_bits(n, m_max) / 62.0) + max_threads + 10;

    mpz_t* remainders = malloc(max_primes * sizeof(mpz_t));
    mpz_t* moduli = malloc(max_primes * sizeof(mpz_t));
    for (int i = 0; i < max_primes; i++) {
        mpz_inits(remainders[i], moduli[i], NULL);
    }

    for (uint64_t m = m_min; m <= m_max; m++) {
        if (n % 2 == 1 && m % 2 == 1) {
            if (!skip_zeros) {
                if (m_min != m_max) printf("%llu\t", m);
                printf("0\n");
                if (output != NULL) fprintf(output, "0\n");
            }
            continue;
        }

        double start_time;
        if (verbose) start_time = omp_get_wtime();

        // 1D logic only requires roots based on the smallest dimension!
        uint64_t k = (n < m) ? n : m;
        uint64_t l = 2 * (k + 1); 
        
        double section_time;
        if (verbose) {
            printf("Finding primes...");
            fflush(stdout);
            section_time = omp_get_wtime();
        }

        int primes_found = 0;
        double target_bits = calculate_target_bits(n, m);
        double collected_bits = 0.0;

        uint64_t base_candidate = (1ull << 62) + (1ull << 61);
        base_candidate = base_candidate + (l - (base_candidate % l)) + 1; 

        uint64_t candidate_offset = 0;

        #pragma omp parallel
        {
            while (1) {
                // Check if we hit our target bits. 
                // Using an atomic read prevents memory cache staleness across cores.
                double current_bits;
                #pragma omp atomic read
                current_bits = collected_bits;
                
                if (current_bits >= target_bits) break;

                // Atomically claim the next candidate offset to test
                uint64_t my_offset;
                #pragma omp atomic capture
                my_offset = candidate_offset++;
                
                uint64_t p_candidate = base_candidate + my_offset * l;
                
                // This is the heavy lifting. All threads do this concurrently!
                if (is_prime(p_candidate)) {
                    
                    // Atomically claim our index in the output array
                    int my_index;
                    #pragma omp atomic capture
                    my_index = primes_found++;
                    
                    mpz_set_ui(moduli[my_index], p_candidate);
                    
                    // Atomically add the bits we just discovered
                    double my_bits = log2((double)p_candidate);
                    #pragma omp atomic
                    collected_bits += my_bits;
                }
            }
        }
        
        int num_primes = primes_found;

        if (verbose) {
            printf("\rFound %d primes. %.3f sec.%10s\n", num_primes, omp_get_wtime() - section_time, "");
            section_time = omp_get_wtime();
        }

        int completed = 0;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < num_primes; i++) {
            uint64_t p_val = mpz_get_ui(moduli[i]); 
            uint64_t r = evaluate_kasteleyn_mod_p(n, m, p_val);
            mpz_set_ui(remainders[i], r);

            int current_completed;
            #pragma omp atomic capture
            current_completed = ++completed;
            
            if (verbose && current_completed % 100 == 0) {
                printf("Evaluating across %d threads. %d / %d primes...\r", max_threads, current_completed, num_primes);
                fflush(stdout);
            }
        }
        
        if (verbose) {
            printf("\rEvaluated across %d threads. %d / %d primes. %.3f sec.%10s\n", max_threads, completed, num_primes, omp_get_wtime() - section_time, "");
            printf("Folding CRT tree...\r");
            section_time = omp_get_wtime();
        }
        
        crt_tree(remainders, moduli, 0, num_primes - 1);
        
        if (verbose) {
            printf("Folded CRT tree. %.3f sec.%10s\n", omp_get_wtime() - section_time, "");
            printf("Time: %.3f sec.%10s\n", omp_get_wtime() - start_time, "");
        }

        if (!quiet) {
            if (m_min != m_max) printf("%llu\t", m);
            mpz_out_str(stdout, 10, remainders[0]);
            printf("\n");
        }

        if (output != NULL) {
            mpz_out_str(output, 10, remainders[0]);
            fprintf(output, "\n");
        }
    }

    for (int i = 0; i < max_primes; i++) {
        mpz_clears(remainders[i], moduli[i], NULL);
    }
    free(remainders);
    free(moduli);

    if (verbose) printf("Total time: %.3f sec.\n", omp_get_wtime() - total_start_time);
    if (output != NULL) fclose(output);

    return 0;
}
