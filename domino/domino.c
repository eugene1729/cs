#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>

#define VERBOSE     1

unsigned int n;
unsigned int m;
unsigned int t;
unsigned int ti;

int** cyclotomics;

#if VERBOSE
unsigned int num_polynomials_multiplied;
unsigned int num_polynomials_to_multiply;
clock_t start_time;
clock_t time_;
#endif

unsigned int
gcd (
    unsigned int n,
    unsigned int m
    )

{
    while (n != 0) {
        unsigned int n_ = m % n;
        m = n;
        n = n_;
    }
    return m;
}

int*
cyclotomic (
    unsigned int n
    )

{
    int* p = cyclotomics[n];
    if (p != 0) return p;

    p = calloc(n + 1, sizeof(int));
    p[0] = -1; 
    p[n] = 1; 
    int pl = n;
    int* q = calloc(n + 1, sizeof(int));

    for (unsigned int i = 1; i < n; i++) {

        if (n % i == 0) {

            int* c = cyclotomic(i);
            int cl;
            for (cl = i; c[cl] == 0; cl--);

            int ql = pl - cl;
            while (pl >= cl) {
                if (p[pl] != 0) {
                    int a = p[pl] / c[cl];
                    q[pl - cl] = a;
                    for (int j = 0; j <= cl; j++) {
                        p[pl - cl + j] -= c[j] * a;
                    }
                }
                pl--;
            }
            // Check that remainder is 0.
            for (int i = 0; i <= pl; i++) assert(p[i] == 0);

            int* z = p;
            p = q;
            q = z;
            pl = ql;
        }
    }

    free(q);

    return cyclotomics[n] = p;
}

typedef struct {
    mpz_t* p;
    mpz_t* p_;
    unsigned int power_of_2;
} POLYNOMIAL_PRODUCT;

#if VERBOSE
void
update_multiplying_polynomials_status ()

{
    num_polynomials_multiplied++;
    clock_t now = clock();
    if (now - time_ > CLOCKS_PER_SEC * 0.1) {
        printf("\rMultiplying polynomials... %u / %u", num_polynomials_multiplied, num_polynomials_to_multiply);
        fflush(stdout);
        time_ = now;
    }
}
#endif

void
swap_polynomials (
    POLYNOMIAL_PRODUCT* pp
    )

{
    mpz_t* t = pp->p;
    pp->p = pp->p_;
    pp->p_ = t;
}

void
multiply_by_polynomial (
    POLYNOMIAL_PRODUCT* pp,
    unsigned int j,
    unsigned int k
    )

{
    unsigned int p_exp;

    // When i == j, polynomial has form [4, ..., 2..., 2, ...], then factor out 2 and multiply by power of 2 at the end.
    if (j == k) {
        pp->power_of_2++;
        p_exp = 1;
    }
    else {
        p_exp = 2;
    }

    // Zero out target polynomial while multiplying by 2^p_exp at the same time.
    for (unsigned int i = 0; i < t; i++) {
        if (mpz_sgn(pp->p[i]) != 0) {
            mpz_mul_2exp(pp->p_[i], pp->p[i], p_exp);
        }
        else if (mpz_sgn(pp->p_[i]) != 0) {
            mpz_set_ui(pp->p_[i], 0);
        }
    }

    for (unsigned int i = 0; i < t; i++) {
        if (mpz_sgn(pp->p[i]) != 0) {
            unsigned int i_;
            i_ = (i + j) % t;
            mpz_add(pp->p_[i_], pp->p_[i_], pp->p[i]);
            i_ = (t + i - j) % t;
            mpz_add(pp->p_[i_], pp->p_[i_], pp->p[i]);
            if (j != k) {
                i_ = (i + k) % t;
                mpz_add(pp->p_[i_], pp->p_[i_], pp->p[i]);
                i_ = (t + i - k) % t;
                mpz_add(pp->p_[i_], pp->p_[i_], pp->p[i]);
            }
        }
    }

    swap_polynomials(pp);

#if VERBOSE
    update_multiplying_polynomials_status();
#endif
}

void
square_polynomial (
    POLYNOMIAL_PRODUCT* pp
    )

{
    for (unsigned int i = 0; i < t; i++) {
        if (mpz_sgn(pp->p_[i]) != 0) {
            mpz_set_ui(pp->p_[i], 0);
        }
    }

    mpz_t p;
    mpz_init(p);
    for (unsigned int i = 1; i < t; i++) {
        for (unsigned int j = 0; j < i; j++) {
            mpz_mul(p, pp->p[i], pp->p[j]);
            mpz_addmul_ui(pp->p_[(i + j) % t], p, 2);
        }
    }
    mpz_clear(p);

    for (unsigned int i = 0; i < t; i++) {
        mpz_addmul(pp->p_[(i * 2) % t], pp->p[i], pp->p[i]);
    }

    swap_polynomials(pp);

#if VERBOSE
    update_multiplying_polynomials_status();
#endif
}

int
main (
    int argc,
    char** argv
    )

{
    // Initialize from command line parameters.
    if (argc != 2 && argc != 3) {
        printf("ERROR: Required parameters are missing.\n");
        return 1;
    }
    n = atoi(argv[1]);
    if (n < 1) {
        printf("ERROR: N is required.\n");
        return 1;
    }

    unsigned int m_min;
    unsigned int m_max;

    if (argc == 3) {
        m_max = m_min = atoi(argv[2]);
        if (m_max < 1) {
            printf("ERROR: Invalid M.\n");
            return 1;
        }
    }
    else {
        m_min = 1;
        m_max = 65536;
    }

    unsigned int max_t = (n + 1) * (m_max + 1);
    cyclotomics = calloc(max_t + 1, sizeof(int*));

    for (m = m_min; m <= m_max; m++) {

        // Handle special case of odd N and odd M.
        if (n % 2 == 1 && m % 2 == 1) {
            if (m_min != m_max) {
                printf("%u\t", m);
            }
            printf("0\n");
            continue;
        }

#if VERBOSE
        time_ = start_time = clock();
#endif

        /*
            Kasteleyn formula is 
            \prod_{j=1}^{\lceil\frac{m}{2}\rceil} \prod_{k=1}^{\lceil\frac{n}{2}\rceil} ( 4 \cos^2 \frac{\pi j}{m+1} + 4 \cos^2 \frac{\pi k}{n+1} )
            Each element in the product is 
            4 \cos^2 \frac{\pi j}{m+1} + 4 \cos^2 \frac{\pi k}{n+1}
            4 + e^{\frac{2 \pi i j}{m+1}} + e^{-\frac{2 \pi i j}{m+1}} + e^{\frac{2 \pi i k}{n+1}} + e^{-\frac{2 \pi i k}{n+1}}
            Let t = (m+1)(n+1)
            Let z = e^{\frac{2 \pi i}{t}} -- t-th root of unity
            4 + z^{j (n+1)} + z^{-j (n+1)} + z^{k (m+1)} + z^{-k (m+1)}
            So each element in the product can be represented by polynomial in z with integer coefficients.
        */

        int square = n == m;
        ti = gcd(n + 1, m + 1);
        t = (n + 1) * (m + 1) / ti;

#if VERBOSE
        printf("Polynomial degree: %u.\n", t);
        num_polynomials_multiplied = 0;
#endif

        mpz_t* p1 = malloc(t * sizeof(mpz_t));
        mpz_t* p2 = malloc(t * sizeof(mpz_t));
        for (unsigned int i = 0; i < t; i++) {
            mpz_init(p1[i]);
            mpz_init(p2[i]);
        }

        POLYNOMIAL_PRODUCT pp;
        pp.power_of_2 = 0;
        pp.p = p1;
        pp.p_ = p2;
        mpz_set_ui(pp.p[0], 1);

        if (square) {
            num_polynomials_to_multiply = ((n + 1) / 2 - 1) * ((n + 1) / 2) / 2 + (n + 1) / 2 + 1;
            // Lower triangle.
            for (unsigned int j = 2; j <= t / 2; j++) {
                for (unsigned int k = 1; k < j; k++) {
                    multiply_by_polynomial(&pp, j, k);
                }
            }
            // Square to get lower and upper triangle (upper is equal to lower).
            square_polynomial(&pp);
            // Diagonal.
            for (unsigned int j = 1; j <= t / 2; j++) {
                multiply_by_polynomial(&pp, j, j);
            }
        }

        else {
            num_polynomials_to_multiply = ((m + 1) / 2) * ((n + 1) / 2); // Not the same as t / 4.
            unsigned int ji = (n + 1) / ti;
            unsigned int ki = (m + 1) / ti;
            for (unsigned int j = ji; j <= t / 2; j += ji) {
                for (unsigned int k = ki; k <= t / 2; k += ki) {
                    multiply_by_polynomial(&pp, j, k);
                }
            }
        }

#if VERBOSE
        printf("\rMultiplying polynomials... %.3f sec.%10s\n", (double) (clock() - start_time) / CLOCKS_PER_SEC, "");
        printf("Calculating cyclotomic polynomial...");
        fflush(stdout);
        start_time = clock();
#endif

        /*
            We now have a polynomial P(x). r = P(z). r is an integer and is the answer we are looking for.
            By Little Bezout's theorem, r is the remainder from dividing P(x) by x-z. P(x) = Q(x)(x-z) + r.
            In fact P(z^j) = P(z) = r for any z^j where j is coprime with t. 
            This is because of the structure of polynomials we multiplied in the first step.
            Product of all z^i where i is in [1 .. t-1] and i is coprime with t is 
            t-th cyclotomic polynomial \Phi_{t}(x) that has all integer coefficients. 
            So the answer is the remainder from division of P(x) by \Phi_{t}(x). 
            Two polynomials with integer coefficients.
        */

        // Calculate cyclotomic polynomial.

        unsigned int t_ = (m + 1) * (n + 1);
        int* c_ = cyclotomic(t_);
        int* c;
        unsigned int cl;
        for (cl = t_; c_[cl] == 0; cl--);
        assert(c_[cl] == 1);
        if (ti != 1) {
            assert(cl % ti == 0);
            cl /= ti;
            c = calloc(cl + 1, sizeof(int));
            for (unsigned int i = 0; i <= cl; i++) {
                c[i] = c_[i * ti];
            }
        }
        else {
            c = c_;
        }

#if VERBOSE
        printf(" %.3f sec.\n", (double) (clock() - start_time) / CLOCKS_PER_SEC);
        printf("Cyclotomic polynomial degree: %u.\n", cl);
        printf("Dividing by cyclotomic polynomial...");
        fflush(stdout);
        time_ = start_time = clock();
#endif

        // Divide by cyclotomic polynomial. The answer is the remainder from that division.

        mpz_t* p = pp.p;
        mpz_t a;
        mpz_init(a);
        int pl = t - 1;
        while (pl >= cl) {
#if VERBOSE
            clock_t now = clock();
            if (now - time_ > CLOCKS_PER_SEC * 0.1) {
                printf("\rDividing by cyclotomic polynomial... %u / %u", t - 1 - pl, t - cl);
                fflush(stdout);
                time_ = now;
            }
#endif
            if (mpz_sgn(p[pl]) != 0) {
                for (unsigned int j = 0; j <= cl; j++) {
                    if (c[j] == 1) {
                        mpz_sub(p[pl - cl + j], p[pl - cl + j], p[pl]);
                    }
                    else if (c[j] == -1) {
                        mpz_add(p[pl - cl + j], p[pl - cl + j], p[pl]);
                    }
                    else if (c[j] > 0) {
                        mpz_submul_ui(p[pl - cl + j], p[pl], c[j]);
                    }
                    else if (c[j] < 0) {
                        mpz_addmul_ui(p[pl - cl + j], p[pl], -c[j]);
                    }
                }
            }
            pl--;
        }
        mpz_clear(a);

        if (c != c_) {
            free(c);
        }

#if VERBOSE
        printf("\rDividing by cyclotomic polynomial... %.3f sec.%10s\n", (double) (clock() - start_time) / CLOCKS_PER_SEC, "");
#endif

        if (m_min != m_max) {
            printf("%u\t", m);
        }
        if (pp.power_of_2 != 0) {
            mpz_mul_2exp(p[0], p[0], pp.power_of_2);
        }
        mpz_out_str(stdout, 10, p[0]);
        printf("\n");

        for (unsigned int i = 0; i < t; i++) {
            mpz_clear(p1[i]);
            mpz_clear(p2[i]);
        }
        free(p1);
        free(p2);
    }

    for (int i = 0; i <= max_t; i++) free(cyclotomics[i]);
    free(cyclotomics);

    return 0;
}
