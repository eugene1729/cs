#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

bool verbose = false;
bool skip_zeros = false;

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
            else if (strcmp(arg, "-v") == 0 || strcmp(arg, "--verbose") == 0) {
                verbose = true;
                continue;
            }
            else if (strcmp(arg, "--skip-zeros") == 0) {
                skip_zeros = true;
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
        return 1;
    }

    if (m_min == 0) m_min = 1;
    if (m_max == 0) m_max = 1000;

    FILE* output = NULL;
    if (output_file_name != NULL) {
        output = fopen(output_file_name, "w");
    }

    // --- DP Setup ---
    int max_mask = 1 << n;
    mpz_t *dp_cur = malloc(max_mask * sizeof(mpz_t));
    mpz_t *dp_next = malloc(max_mask * sizeof(mpz_t));
    for (int i = 0; i < max_mask; i++) {
        mpz_init(dp_cur[i]);
        mpz_init(dp_next[i]);
    }
    mpz_set_ui(dp_cur[0], 1);

    // --- Berlekamp-Massey Setup ---
    int max_seq = 50000; 
    mpq_t *C = malloc(max_seq * sizeof(mpq_t));
    mpq_t *B = malloc(max_seq * sizeof(mpq_t));
    mpq_t *T = malloc(max_seq * sizeof(mpq_t));
    
    for (int i = 0; i < max_seq; i++) {
        mpq_init(C[i]); mpq_set_ui(C[i], 0, 1);
        mpq_init(B[i]); mpq_set_ui(B[i], 0, 1);
        mpq_init(T[i]); mpq_set_ui(T[i], 0, 1);
    }
    mpq_set_ui(C[0], 1, 1);
    mpq_set_ui(B[0], 1, 1);

    int L = 0;
    int m_var = 1;
    int L_B = 0;
    int consecutive_zeros = 0;
    bool locked = false;

    mpq_t b, d, c_val, temp;
    mpq_inits(b, d, c_val, temp, NULL);
    mpq_set_ui(b, 1, 1);

    mpz_t *coeff = NULL;
    mpz_t *x = malloc(m_max * sizeof(mpz_t));
    for (uint64_t i = 0; i < m_max; i++) {
        mpz_init(x[i]);
    }

    // --- Core Generation Loop ---
    for (uint64_t step = 0; step < m_max; step++) {
        
        if (!locked) {
            // Generate next column using square-by-square DP
            for (int r = 0; r < n; r++) {
                for (int mask = 0; mask < max_mask; mask++) {
                    mpz_set_ui(dp_next[mask], 0);
                }

                for (int mask = 0; mask < max_mask; mask++) {
                    if (mpz_sgn(dp_cur[mask]) == 0) continue;

                    if ((mask & 1) != 0) {
                        // Square occupied
                        int nmask = mask >> 1;
                        mpz_add(dp_next[nmask], dp_next[nmask], dp_cur[mask]);
                    } 
                    else {
                        // Square empty - Option 1: Horizontal Domino
                        int nmask_h = (mask >> 1) | (1 << (n - 1));
                        mpz_add(dp_next[nmask_h], dp_next[nmask_h], dp_cur[mask]);

                        // Square empty - Option 2: Vertical Domino
                        if (r < n - 1 && (mask & 2) == 0) {
                            int nmask_v = (mask >> 1) | 1;
                            mpz_add(dp_next[nmask_v], dp_next[nmask_v], dp_cur[mask]);
                        }
                    }
                }
                
                mpz_t *t_ptr = dp_cur; 
                dp_cur = dp_next; 
                dp_next = t_ptr;
            }

            mpz_set(x[step], dp_cur[0]);

            // Feed new value into BM Algorithm
            if (step < max_seq) {
                mpq_set_z(d, x[step]);
                for (int i = 1; i <= L; i++) {
                    mpq_set_z(temp, x[step - i]);
                    mpq_mul(temp, temp, C[i]);
                    mpq_add(d, d, temp);
                }

                if (mpq_sgn(d) == 0) {
                    m_var++;
                    consecutive_zeros++;
                    
                    // If prediction is perfect for 15 straight iterations, check if we've found the integer recurrence
                    if (consecutive_zeros >= 15 && L > 0) {
                        bool integer_coeffs = true;
                        for (int i = 1; i <= L; i++) {
                            if (mpz_cmp_ui(mpq_denref(C[i]), 1) != 0) {
                                integer_coeffs = false;
                                break;
                            }
                        }
                        
                        if (integer_coeffs) {
                            locked = true;
                            coeff = malloc((L + 1) * sizeof(mpz_t));
                            for (int i = 1; i <= L; i++) {
                                mpz_init(coeff[i]);
                                // Store exactly as x_n = SUM(coeff[i] * x_{n-i})
                                mpz_set(coeff[i], mpq_numref(C[i]));
                                mpz_neg(coeff[i], coeff[i]);
                            }
                            if (verbose) {
                                printf("Berlekamp-Massey target acquired! Linear recurrence of degree %d found at M = %llu. Switching to overdrive.\n", L, step + 1);
                            }
                        }
                    }
                } 
                else {
                    consecutive_zeros = 0;
                    for (int i = 0; i <= L; i++) mpq_set(T[i], C[i]);
                    int T_len = L;

                    mpq_div(c_val, d, b);

                    int new_L_bound = L > (L_B + m_var) ? L : (L_B + m_var);
                    for (int i = L + 1; i <= new_L_bound; i++) {
                        mpq_set_ui(C[i], 0, 1);
                    }

                    for (int i = 0; i <= L_B; i++) {
                        mpq_mul(temp, c_val, B[i]);
                        mpq_sub(C[i + m_var], C[i + m_var], temp);
                    }

                    if (2 * L <= step) {
                        L = step + 1 - L;
                        for (int i = 0; i <= T_len; i++) mpq_set(B[i], T[i]);
                        L_B = T_len;
                        mpq_set(b, d);
                        m_var = 1;
                    } 
                    else {
                        m_var++;
                    }
                }
            }
        } 
        else {
            // BM Locked: Fast sequence generation bypassing DP completely
            mpz_set_ui(x[step], 0);
            for (int i = 1; i <= L; i++) {
                mpz_addmul(x[step], coeff[i], x[step - i]);
            }
        }

        // Print Logic
        uint64_t current_m = step + 1;
        if (!skip_zeros || mpz_sgn(x[step]) != 0) {
            if (current_m >= m_min) {
                if (m_min != m_max) printf("%llu\t", current_m);
                mpz_out_str(stdout, 10, x[step]);
                printf("\n");
                
                if (output != NULL) {
                    mpz_out_str(output, 10, x[step]);
                    fprintf(output, "\n");
                }
            }
        }
    }

    // Cleanup
    for (int i = 0; i < max_mask; i++) {
        mpz_clears(dp_cur[i], dp_next[i], NULL);
    }
    free(dp_cur);
    free(dp_next);

    for (int i = 0; i < max_seq; i++) {
        mpq_clears(C[i], B[i], T[i], NULL);
    }
    free(C); 
    free(B); 
    free(T);
    
    mpq_clears(b, d, c_val, temp, NULL);
    
    if (coeff != NULL) {
        for (int i = 1; i <= L; i++) mpz_clear(coeff[i]);
        free(coeff);
    }

    for (uint64_t i = 0; i < m_max; i++) {
        mpz_clear(x[i]);
    }
    free(x);
    
    if (output != NULL) fclose(output);

    return 0;
}
