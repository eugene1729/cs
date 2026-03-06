#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

typedef unsigned __int128 uint128_t;

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
        if (e % 2 == 1) {
            r = mul_mod(r, b, p);
        }
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
    if (x == 1 || x == p - 1) {
        return true;
    }
    while (d != p - 1) {
        x = mul_mod(x, x, p);
        d *= 2;
        if (x == 1) { 
            return false; 
        }
        if (x == p - 1) {
            return true;
        }
    }
    return false;
}

bool 
is_prime (
    uint64_t p
    ) 

{
    if (p <= 1 || p == 4) {
        return false;
    }
    if (p <= 3) {
        return true;
    }

    uint64_t d = p - 1;
    while (d % 2 == 0) {
        d /= 2;
    }
    
    // Deterministic bases for 64-bit integers.
    static const uint64_t bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (int i = 0; i < 12; i++) {
        if (p <= bases[i]) {
            break;
        }
        if (!miller_rabin_test(d, p, bases[i])) {
            return false;
        }
    }
    return true;
}

//
// Number Theory Helpers
//

uint64_t 
gcd (
    uint64_t a, 
    uint64_t b
    ) 

{
    while (b != 0) {
        uint64_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

uint64_t 
lcm (
    uint64_t a, 
    uint64_t b
    ) 

{
    return (a / gcd(a, b)) * b;
}

// Finds k-th root of unity modulo p.

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
        if (z == 1) {
            continue;
        }
        
        bool is_primitive = true;
        for (int i = 0; i < num_factors; i++) {
            if (pow_mod(z, k / factors[i], p) == 1) {
                is_primitive = false;
                break;
            }
        }
        if (is_primitive) {
            return z;
        }
    }
    return 0; 
}

uint64_t 
evaluate_kasteleyn_mod_p (
    uint64_t n, 
    uint64_t m, 
    uint64_t p
    ) 

{
    uint64_t r;

    // Impossible to perfectly tile an odd-by-odd grid.
    if (n % 2 == 1 && m % 2 == 1) {
        r = 0;
    }

    // Fast path for square grids.
    else if (n == m) {
        uint64_t t_n = 2 * (n + 1);
        uint64_t w_n = find_primitive_root(t_n, p);
        uint64_t n_half = (n + 1) / 2;
        
        uint64_t* a = malloc((n_half + 1) * sizeof(uint64_t));
        
        // Setup Chebyshev Recurrence: W_j = w^(2j) + w^(-2j)
        uint64_t w2 = mul_mod(w_n, w_n, p);
        uint64_t inv_w2 = inv_mod(w2, p);
        
        uint64_t w1 = add_mod(w2, inv_w2, p);
        uint64_t w_prev = 2;       // W_0 = w^0 + w^0 = 2
        uint64_t w_curr = w1;      // W_1 = w^2 + w^-2
        
        a[1] = add_mod(w1, 2, p);
        
        // Generate the 1D array using exactly 1 MUL, 1 SUB, 1 ADD per step
        for (uint64_t j = 2; j <= n_half; j++) {
            uint64_t w_next = sub_mod(mul_mod(w1, w_curr, p), w_prev, p);
            a[j] = add_mod(w_next, 2, p);
            w_prev = w_curr;
            w_curr = w_next;
        }
        
        uint64_t diag = 1;
        uint64_t off_diag = 1;
        
        // Matrix Symmetry Evaluation (Upper Triangle Only)
        for (uint64_t j = 1; j <= n_half; j++) {
            
            // 1. The Diagonal (j == k)
            uint64_t term_diag = add_mod(a[j], a[j], p);
            diag = mul_mod(diag, term_diag, p);
            
            // 2. The Upper Triangle (k > j)
            for (uint64_t k = j + 1; k <= n_half; k++) {
                uint64_t term = add_mod(a[j], a[k], p);
                off_diag = mul_mod(off_diag, term, p);
            }
        }

        free(a);
        
        // Square the off-diagonal product to account for the lower triangle
        off_diag = mul_mod(off_diag, off_diag, p);
        r = mul_mod(diag, off_diag, p);
    }

    // Regular path for non-square grids.
    else {
        uint64_t t_m = 2 * (m + 1);
        uint64_t t_n = 2 * (n + 1);
        
        uint64_t w_m = find_primitive_root(t_m, p);
        uint64_t w_n = find_primitive_root(t_n, p);
        
        uint64_t m_half = (m + 1) / 2;
        uint64_t n_half = (n + 1) / 2;
        
        uint64_t* a = malloc((m_half + 1) * sizeof(uint64_t));
        uint64_t* b = malloc((n_half + 1) * sizeof(uint64_t));
        
        // Setup Recurrence for a (Dimension M)
        uint64_t w2_m = mul_mod(w_m, w_m, p);
        uint64_t inv_w2_m = inv_mod(w2_m, p);
        uint64_t w1_m = add_mod(w2_m, inv_w2_m, p);
        uint64_t w_prev_m = 2, w_curr_m = w1_m;
        
        a[1] = add_mod(w1_m, 2, p);
        for (uint64_t j = 2; j <= m_half; j++) {
            uint64_t w_next_m = sub_mod(mul_mod(w1_m, w_curr_m, p), w_prev_m, p);
            a[j] = add_mod(w_next_m, 2, p);
            w_prev_m = w_curr_m;
            w_curr_m = w_next_m;
        }
        
        // Setup Recurrence for b (Dimension N)
        uint64_t w2_n = mul_mod(w_n, w_n, p);
        uint64_t inv_w2_n = inv_mod(w2_n, p);
        uint64_t w1_n = add_mod(w2_n, inv_w2_n, p);
        uint64_t w_prev_n = 2, w_curr_n = w1_n;
        
        b[1] = add_mod(w1_n, 2, p);
        for (uint64_t k = 2; k <= n_half; k++) {
            uint64_t w_next_n = sub_mod(mul_mod(w1_n, w_curr_n, p), w_prev_n, p);
            b[k] = add_mod(w_next_n, 2, p);
            w_prev_n = w_curr_n;
            w_curr_n = w_next_n;
        }
        
        r = 1;
        // O(N*M) Double Loop
        for (uint64_t j = 1; j <= m_half; j++) {
            for (uint64_t k = 1; k <= n_half; k++) {
                uint64_t term = add_mod(a[j], b[k], p);
                r = mul_mod(r, term, p);
            }
        }
        
        free(a);
        free(b);
    }

    return r;
}

//
// Divide-and-Conquer CRT using GMP
//

// Combines (r0 mod p0) and (r1 mod p1) in-place into r0 and p0.
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
    
    // inv = p0^-1 mod p1
    mpz_invert(inv, p0, p1); 
    
    // t = (r1 - r0) * inv mod p1
    mpz_sub(t, r1, r0);
    mpz_mul(t, t, inv);
    mpz_mod(t, t, p1); 
    
    // P = p0 * p1
    mpz_mul(P, p0, p1);
    
    // r0 = r0 + p0 * t
    mpz_mul(t, t, p0);
    mpz_add(r0, r0, t);
    mpz_mod(r0, r0, P); 
    
    // p0 = P
    mpz_set(p0, P);
    
    mpz_clears(inv, t, P, NULL);
}

// Recursively builds the CRT tree.
void 
crt_tree (
    mpz_t* r, 
    mpz_t* p, 
    int start, 
    int end
    ) 

{
    if (start == end) {
        return;
    }
    if (start + 1 == end) {
        crt_combine(r[start], p[start], r[end], p[end]);
        return;
    }
    int mid = start + (end - start) / 2;
    crt_tree(r, p, start, mid);
    crt_tree(r, p, mid + 1, end);
    crt_combine(r[start], p[start], r[mid + 1], p[mid + 1]);
}

double
calculate_target_bits (
    uint64_t n,
    uint64_t m
    )

{
    // Calculate the EXACT bit bound using floating-point Kasteleyn.
    // This runs in O(N*M) but uses fast floats, taking < 5ms even for massive grids.
    double exact_bits = 0.0;
    uint64_t m_half = m / 2; 
    uint64_t n_half = n / 2;
    for (uint64_t j = 1; j <= m_half; j++) {
        double cos_j = cos(M_PI * j / (m + 1.0));
        double term_j = 4.0 * cos_j * cos_j;
        for (uint64_t k = 1; k <= n_half; k++) {
            double cos_k = cos(M_PI * k / (n + 1.0));
            double term_k = 4.0 * cos_k * cos_k;
            exact_bits += log2(term_j + term_k);
        }
    }
    
    // Add 64 bits to absorb any floating-point ulp drift and provide a bulletproof CRT cushion.
    return exact_bits + 64.0;
}

//
// Main
//

int main (
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
                if (m_min < 1) {
                    printf("ERROR: Minimum M is invalid.\n");
                    return 1;
                }
                continue;
            }
            else if (i < argc - 1 && (strcmp(arg, "-m") == 0 || strcmp(arg, "--max") == 0)) {
                m_max = strtoull(argv[++i], NULL, 10);
                if (m_max < 1) {
                    printf("ERROR: Maximum M is invalid.\n");
                    return 1;
                }
                continue;
            }
            else if (i < argc - 1 && (strcmp(arg, "-o") == 0 || strcmp(arg, "--output") == 0)) {
                output_file_name = argv[++i];
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
            if (n < 1) {
                printf("ERROR: N is invalid.\n");
                return 1;
            }
            continue;
        }
        else if (m_min == 0) {
            m_max = m_min = strtoull(arg, NULL, 10);
            if (m_max < 1) {
                printf("ERROR: M is invalid.\n");
                return 1;
            }
            continue;
        }

        printf("ERROR: Unexpected argument '%s'.\n", arg);
        return 1;
    }

    if (n == 0) {
        printf("Usage: %s <N> <M> [--max <M_max>] [--verbose] [--skip-zeros] [--output <file name>]\n", argv[0]);
    }

    if (m_min == 0) {
        m_min = 1;
    }
    if (m_max == 0) {
        m_max = 1000;
    }

    FILE* output = NULL;
    if (output_file_name != NULL) {
        output = fopen(output_file_name, "w");
    }

    double total_start_time;
    if (verbose) {
        total_start_time = omp_get_wtime();
    }

    // Calculate Required Primes based on max possible answer size.
    // The number of tilings is roughly e^(G/pi * n * m), bounded by 2^(n*m/2).
    int max_primes = (int)(calculate_target_bits(n, m_max) / 62.0) + 2;

    mpz_t* remainders = malloc(max_primes * sizeof(mpz_t));
    mpz_t* moduli = malloc(max_primes * sizeof(mpz_t));
    for (int i = 0; i < max_primes; i++) {
        mpz_inits(remainders[i], moduli[i], NULL);
    }

    for (uint64_t m = m_min; m <= m_max; m++) {

        // Handle special case of odd N and odd M.
        if (n % 2 == 1 && m % 2 == 1) {
            if (!skip_zeros) {
                if (m_min != m_max) {
                    printf("%llu\t", m);
                }
                printf("0\n");
                if (output != NULL) {
                    fprintf(output, "0\n");
                }
            }
            continue;
        }

        double start_time;
        if (verbose) {
            start_time = omp_get_wtime();
        }

        // Calculate Required Primes.

        // Find the LCM of 2(n+1) and 2(m+1) to ensure roots exist.
        uint64_t l = lcm(2 * (n + 1), 2 * (m + 1));
        
        // Prime Sieve.
        double section_time;
        if (verbose) {
            printf("Finding primes...");
            fflush(stdout);
            section_time = omp_get_wtime();
        }

        int primes_found = 0;
        double target_bits = calculate_target_bits(n, m);
        double collected_bits = 0.0;

        // Start high in the 63-bit range to pack maximum bits per prime.
        // We start at 2^62 + 2^61 (giving ~62.6 bits per prime) 
        // which safely avoids the uint64_t addition overflow limit of 2^63 - 1.
        uint64_t p_candidate = (1ull << 62) + (1ull << 61);
        p_candidate = p_candidate + (l - (p_candidate % l)) + 1; 
        
        while (collected_bits < target_bits) {
            if (is_prime(p_candidate)) {
                mpz_set_ui(moduli[primes_found++], p_candidate);
                // Track the EXACT bit capacity of this specific prime
                collected_bits += log2((double)p_candidate); 
            }
            p_candidate += l;
        }
        
        int num_primes = primes_found; // Lock in the final count for your OpenMP loop

        if (verbose) {
            printf("\rFound %d primes. %.3f sec.%10s\n", num_primes, omp_get_wtime() - section_time, "");
        }

        // Kasteleyn Evaluation (Parallel - the heavy lifting)
        if (verbose) {
            section_time = omp_get_wtime();
        }
        int completed = 0;

        // This single line spins up the thread pool and divides the loop evenly.
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < num_primes; i++) {
            // Grab the prime for this specific thread's index
            uint64_t p_val = mpz_get_ui(moduli[i]); 
            
            // Do the heavy 2D loop
            uint64_t r = evaluate_kasteleyn_mod_p(n, m, p_val);
            
            // Store the result. GMP's mpz_set_ui is thread-safe here 
            // because each thread writes to a strictly distinct array index (remainders[i]).
            mpz_set_ui(remainders[i], r);

            // Safely increment and capture the exact count in one atomic CPU instruction
            int current_completed;
            #pragma omp atomic capture
            current_completed = ++completed;
            
            // Any thread that crosses a 100-mark milestone updates the screen
            if (verbose && current_completed % 100 == 0) {
                printf("Evaluating across %d threads. %d / %d primes...\r", 
                       omp_get_max_threads(), current_completed, num_primes);
                fflush(stdout);
            }
        }
        if (verbose) {
            printf("\rEvaluated across %d threads. %d / %d primes. %.3f sec.%10s\n", omp_get_max_threads(), completed, num_primes, omp_get_wtime() - section_time, "");
        }

        // Reconstruct using Divide-and-Conquer CRT.
        if (verbose) {
            printf("Folding CRT tree...\r");
            section_time = omp_get_wtime();
        }
        crt_tree(remainders, moduli, 0, num_primes - 1);
        if (verbose) {
            printf("Folded CRT tree. %.3f sec.%10s\n", omp_get_wtime() - section_time, "");
        }

        // Output the exact answer.
        if (verbose) {
            printf("Time: %.3f sec.%10s\n", omp_get_wtime() - start_time, "");
        }

        if (m_min != m_max) {
            printf("%llu\t", m);
        }
        mpz_out_str(stdout, 10, remainders[0]);
        printf("\n");
        if (output != NULL) {
            mpz_out_str(output, 10, remainders[0]);
            fprintf(output, "\n");
        }
    }

    // Cleanup
    for (int i = 0; i < max_primes; i++) {
        mpz_clears(remainders[i], moduli[i], NULL);
    }
    free(remainders);
    free(moduli);

    if (verbose) {
        printf("Total time: %.3f sec.\n", omp_get_wtime() - total_start_time);
    }

    if (output != NULL) {
        fclose(output);
    }

    return 0;
}
