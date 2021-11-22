#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef unsigned long long ull;
typedef unsigned int uint;

const ull SENTINEL = (ull) -1;

uint n;
uint k;
uint mask_size_in_bits;
ull all_ones;

typedef struct {
    ull m;
    ull mv[6];
    ull l;
    ull i;
    ull c;
} MASK;

MASK* masks;
uint max_masks;
uint num_masks;
MASK** heap;

void
search (
    int i,
    ull p,
    ull m
    )

{
    int x = i % n;
    int y = i / n;

    if (x == n - 1 && y == k - 1) {
        if (num_masks == max_masks) {
            max_masks *= 2;
            masks = realloc(masks, sizeof(MASK) * max_masks);
        }
        masks[num_masks++].m = all_ones ^ m;
    }
    else {
        if (x > 0) {
            int j = i - 1;
            if ((p & (1ull << j)) == 0) search(j, p | (1ull << j), m | (1ull << (y * (n - 1) + x - 1)));
        }
        if (x < n - 1) {
            int j = i + 1;
            if ((p & (1ull << j)) == 0) search(j, p | (1ull << j), m | (1ull << (y * (n - 1) + x)));
        }
        if (y > 0) {
            int j = i - n;
            if ((p & (1ull << j)) == 0) search(j, p | (1ull << j), m | (1ull << (k * (n - 1) + (y - 1) * n + x)));
        }
        if (y < k - 1) {
            int j = i + n;
            if ((p & (1ull << j)) == 0) search(j, p | (1ull << j), m | (1ull << (k * (n - 1) + y * n + x)));
        }
    }
}

int
main (
    int argc,
    char* argv[]
    ) 

{
    if (argc != 3) {
        printf("ERROR: Required parameters are missing.\n");
        exit(1);
    }
    n = atoi(argv[1]);
    if (n < 1 || n > 19) {
        printf("ERROR: N needs to be between 1 and 19.\n");
        exit(1);
    }
    k = atoi(argv[2]);
    if (k < 1 || k > 19) {
        printf("ERROR: K needs to be between 1 and 19.\n");
        exit(1);
    }

    mask_size_in_bits = n * (k - 1) + k * (n - 1);
    all_ones = (1ull << mask_size_in_bits) - 1;

    printf("%ux%u\n", n, k);
    printf("mask_size_in_bits = %u\n", mask_size_in_bits);
    if (mask_size_in_bits > 64) {
        printf("ERROR: Mask is larger than 64 bits.\n");
        exit(1);
    }

    max_masks = 16 * 1024;
    masks = malloc(sizeof(MASK) * max_masks);
    search(0, 1, 0);
    printf("num_masks = %u\n", num_masks);

    masks = realloc(masks, sizeof(MASK) * num_masks);
    heap = malloc(sizeof(MASK*) * num_masks);
    ull t = 0;
    for (int i = 0; i < num_masks; i++) {
        MASK* mask = masks + i;
        ull m = mask->m;
        mask->l = 1ull << __builtin_popcountll(m);
        t += mask->l;
        ull mk = ~m << 1;
        for (int j = 0; j < 6; j++) {
            ull mp = mk;
            mp ^= mp <<  1;
            mp ^= mp <<  2;
            mp ^= mp <<  4;
            mp ^= mp <<  8;
            mp ^= mp << 16;
            mp ^= mp << 32;
            ull mv_ = mp & m;
            mask->mv[j] = mv_;
            m = (m ^ mv_) | (mv_ >> (1ull << j));
            mk &= ~mp;
        }
        heap[i] = mask;
    }

    time_t start_time = time(NULL);
    ull y_ = 0;
    ull y = 0;
    ull z = 0;
    ull last = SENTINEL;
    for (;;) {

        MASK* m = heap[0];
        ull min = m->c;
        if (min == SENTINEL) {
            break;
        }
        y++;
        if (min != last) {
            last = min;
            z++;
        }
        if (y >= y_) {
            time_t cur_time = time(NULL);
            time_t elapsed_time = cur_time - start_time;
            uint h = elapsed_time / 3600;
            uint m = (elapsed_time % 3600) / 60;
            uint s = elapsed_time % 60;
            time_t eta = start_time + (double) elapsed_time * t / y;
            char* eta_string = ctime(&eta);
            eta_string[24] = ' ';
            printf("%0*llx | %16llu | %6.2f%% | Elapsed: %4u:%02u:%02u | ETA: %s\r", (mask_size_in_bits + 3) / 4, min, z, (double) y / t * 100, h, m, s, eta_string);
            fflush(stdout);
            y_ += 50000000;
        }

        m->i++;
        if (m->i >= m->l) {
            m->c = SENTINEL;
        }
        else {
            ull x = m->i;
            ull* a = m->mv;
            x = (x & ~a[5]) | ((x << 32) & a[5]);
            x = (x & ~a[4]) | ((x << 16) & a[4]);
            x = (x & ~a[3]) | ((x <<  8) & a[3]);
            x = (x & ~a[2]) | ((x <<  4) & a[2]);
            x = (x & ~a[1]) | ((x <<  2) & a[1]);
            x = (x & ~a[0]) | ((x <<  1) & a[0]);
            m->c = x & m->m;
        }

        int i = 0;
        for (;;) {
            int l = i * 2 + 1;
            int r = l + 1;
            int min = i;
            if (l < num_masks && heap[l]->c < heap[min]->c) {
                min = l;
            }
            if (r < num_masks && heap[r]->c < heap[min]->c) {
                min = r;
            }
            if (min == i) {
                break;
            }
            MASK* t = heap[min];
            heap[min] = heap[i];
            heap[i] = t;
            i = min;
        }
    }

    time_t cur_time = time(NULL);
    time_t elapsed_time = cur_time - start_time;
    uint h = elapsed_time / 3600;
    uint m = (elapsed_time % 3600) / 60;
    uint s = elapsed_time % 60;
    printf("%-100s\rElapsed: %4u:%02u:%02u\n", "", h, m, s);
    printf("%ux%u\t%llu\n", n, k, z);

    return 0;
}
