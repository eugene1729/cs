#include <fstream>
#include <iostream>
#include <map>
#include <gmpxx.h>

using namespace std;

#define max_w           24
#define log_max_h       10

typedef uint64_t cache_key_t;
typedef map<cache_key_t, mpz_class> cache_t;

uint w;
cache_t cache;

uint
reverse_bits (
    uint m
    )

{
    uint m_ = 0;
    for (int i = 0; i < w; i++) {
        m_ |= (m & 1) ? (1 << (w - 1 - i)) : 0;
        m >>= 1;
    }
    return m_;
}

mpz_class
t (
    int h,                  // remaining height
    int m = 0,              // mask
    int i = 0,              // offset
    int m_ = 0              // new mask
    )

{
    if (i == w) {
        if (h == 1) {
            return m_ == 0 ? 1 : 0;
        }
        uint m_r = reverse_bits(m_);
        if (m_r < m_) {
            m_ = m_r;
        }
        uint64_t k = (((cache_key_t) m_) << log_max_h) | (h - 1);
        map<cache_key_t, mpz_class>::iterator it = cache.find(k);
        if (it != cache.end()) {
            return it->second;
        }
        return cache[k] = t(h - 1, m_);
    }

    mpz_class sum = t(h, m, i + 1, m_ | (m & (1 << i) ? 0 : (1 << i)));

    if (i < w - 1 && (m & (3 << i)) == 0) {
        sum = sum + t(h, m, i + 2, m_);
    }

    return sum;
}

int
main (
    int argc,
    char** argv
    )

{
    if (argc != 2) {
        printf("ERROR: Required parameter is missing.\n");
        return 1;
    }
    w = atoi(argv[1]);
    if (w < 1 || w > max_w) {
        printf("ERROR: Width is required.\n");
        return 1;
    }

    for (int h = 1; h <= 1 << log_max_h; h++) {
        clock_t time = clock();
        mpz_class r = t(h);
        time = clock() - time;
        cout << "D(" << setfill(' ') << setw(2) << w << ", " << setw(2) << h << ") = " << r << " [" << std::fixed << setprecision(3) << (double) time / CLOCKS_PER_SEC << "s]" << endl;
    }
}
