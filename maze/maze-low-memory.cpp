#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <atomic>
#include <cstdlib>
#include <ctime>
#include <omp.h>
#include <gmpxx.h>

typedef mpz_class bigint;
typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned __int128 ulong128; // 5-bit component IDs

const uchar MAX_N = 25; 

uint n;
uint max_k;

// --- FAST 128-BIT HASHING ---
struct Hash128 {
    size_t operator()(const ulong128& x) const {
        uint64_t hi = (uint64_t)(x >> 64);
        uint64_t lo = (uint64_t)x;
        return hi ^ (lo * 11400714819323198485llu);
    }
};

// --- SHARDED MAP FOR ZERO-CONTENTION STATE DISCOVERY ---
const int NUM_SHARDS = 1024;
std::unordered_map<ulong128, uint, Hash128> xs_index[NUM_SHARDS];
omp_lock_t shard_locks[NUM_SHARDS];
std::atomic<uint> global_xsi(2); // ID 1 is the starting state

ulong128* xs_reverse_index;

uint get_or_create_xsi(ulong128 xs) {
    uint shard = Hash128()(xs) % NUM_SHARDS;
    
    omp_set_lock(&shard_locks[shard]);
    auto it = xs_index[shard].find(xs);
    if (it != xs_index[shard].end()) {
        uint id = it->second;
        omp_unset_lock(&shard_locks[shard]);
        return id;
    }
    
    uint new_id = global_xsi.fetch_add(1);
    xs_index[shard][xs] = new_id;
    xs_reverse_index[new_id] = xs;
    omp_unset_lock(&shard_locks[shard]);
    return new_id;
}

// --- THREAD ISOLATION ---
struct ThreadContext {
    ulong128 xs;
    uchar counts[MAX_N + 1];
    ulong128* row;
    ulong128* row_end;

    ThreadContext(uint n_val) {
        row = new ulong128[1 << n_val]; // 2^n max branches
        for (int i = 0; i <= MAX_N; i++) counts[i] = 0;
    }
    ~ThreadContext() {
        delete[] row;
    }
};

// --- CORE RECURSIVE TRANSITIONS ---
void next_v_rec(ThreadContext& ctx, uint i = 0, uint q = 0) {
    if (i == n) {
        if (ctx.counts[0] == 0) return;
        uchar map[MAX_N];
        for (uint j = 0; j < n; j++) map[j] = MAX_N;

        uchar x_ = 1;
        ulong128 xs__ = 0;
        for (uint j = 0; j < n; j++) {
            if ((q & (1 << j)) != 0) {
                uchar x = (ctx.xs >> (j * 5)) & 0x1f; 
                if (x != 0) {
                    if (ctx.counts[x] == 1) x = x_++;
                    else {
                        if (map[x] == MAX_N) map[x] = x_++;
                        x = map[x];
                    }
                    xs__ |= ((ulong128)x) << (j * 5);
                }
            } 
            else {
                xs__ |= ((ulong128)x_) << (j * 5);
                x_++;
            }
        }
        *ctx.row_end++ = xs__;
        return;
    }

    next_v_rec(ctx, i + 1, q);
    uint x = (ctx.xs >> (i * 5)) & 0x1f;
    ctx.counts[x]++;
    next_v_rec(ctx, i + 1, q | (1 << i));
    ctx.counts[x]--;
}

void next_h_rec(ThreadContext& ctx, uint i = 0, uint q = 0) {
    if (i == n - 1) {
        uchar xs_[MAX_N];
        for (uint j = 0; j < n; j++) {
            xs_[j] = (ctx.xs >> (j * 5)) & 0x1f;
            ctx.counts[j] = 0;
        }

        uchar map[MAX_N];
        for (;;) {
            for (uint j = 0; j < n; j++) map[j] = MAX_N;
            bool changed = false;

            for (uint j = 0; j < n - 1; j++) {
                if ((q & (1 << j)) != 0) {
                    uchar x = xs_[j];
                    uchar x_ = xs_[j + 1];
                    if (x < x_) {
                        if (map[x_] > x) { map[x_] = x; changed = true; }
                    } 
                    else if (x_ < x) {
                        if (map[x] > x_) { map[x] = x_; changed = true; }
                    }
                }
            }
            if (!changed) break;
            
            for (uint j = 0; j < n; j++) {
                uchar x = xs_[j];
                if (map[x] != MAX_N) {
                    do { x = map[x]; } while (map[x] != MAX_N);
                    xs_[j] = x;
                }
            }
        }

        for (uint j = 0; j < n; j++) ctx.counts[xs_[j]]++;
        if (ctx.counts[0] == 0) return;

        for (uint j = 0; j < n; j++) map[j] = MAX_N;

        uchar x_ = 1;
        ulong128 xs__ = 0;
        for (uint j = 0; j < n; j++) {
            uchar x = xs_[j];
            if (x != 0) {
                if (ctx.counts[x] == 1) x = x_++;
                else {
                    if (map[x] == MAX_N) map[x] = x_++;
                    x = map[x];
                }
                xs__ |= ((ulong128)x) << (j * 5);
            }
        }
        for (uint j = 0; j < n; j++) ctx.counts[j] = 0;

        *ctx.row_end++ = xs__;
        return;
    }

    next_h_rec(ctx, i + 1, q);
    next_h_rec(ctx, i + 1, q | (1 << i));
}

// OEIS A000245
size_t 
dimension (
    uint i
    ) 

{
    ulong r_val = 3;
    for (uint j = 1; j <= i - 1; j++) {
        r_val = (r_val * (2 * i - j + 1)) / j;
    }
    r_val /= i + 2;
    return r_val;
}

int 
main (
    int argc, 
    char* argv[]
    ) 

{
    max_k = 1000;
    char* output_filename = NULL;

    for (int i = 1; i < argc; i++) {
        char* arg = argv[i];
        if (arg[0] == '-') {
            if (i < argc - 1 && strcmp(arg, "-o") == 0) output_filename = argv[++i];
            else return 1;
        } 
        else if (n == 0) {
            n = atoi(arg);
            if (n < 1 || n > MAX_N) return 1;
        } 
        else {
            max_k = atoi(arg);
        }
    }
    if (n == 0) return 1;

    size_t dim = dimension(n);
    std::cout << "Allocating memory for up to " << dim << " dynamic states..." << std::endl;

    xs_reverse_index = new ulong128[dim + 1];
    bigint* r = new bigint[dim + 1];
    bigint* r_ = new bigint[dim + 1];

    omp_lock_t* r_locks = new omp_lock_t[dim + 1];
    for (size_t i = 0; i <= dim; i++) omp_init_lock(&r_locks[i]);
    for (int i = 0; i < NUM_SHARDS; i++) omp_init_lock(&shard_locks[i]);

    // Initial state injection
    ulong128 xs_init = 0;
    for (uint j = 0; j < n; j++) xs_init |= ((ulong128)j) << (j * 5);
    
    xs_index[Hash128()(xs_init) % NUM_SHARDS][xs_init] = 1;
    xs_reverse_index[1] = xs_init;
    r[1] = 1;

    std::ofstream output_file;
    if (output_filename != NULL) output_file.open(output_filename, std::ofstream::out);

    int max_threads = omp_get_max_threads();
    std::vector<ThreadContext*> ctx_arr(max_threads);
    for (int t = 0; t < max_threads; t++) ctx_arr[t] = new ThreadContext(n);

    clock_t time_ = clock();

    // --- ON THE FLY DYNAMIC PROGRAMMING ---
    for (uint current_k = 1; current_k <= max_k; current_k++) {
        
        // Only loop over states we have discovered so far!
        uint current_max_xsi = global_xsi.load() - 1;

        // --- HORIZONTAL PASS ---
        #pragma omp parallel for schedule(dynamic, 32)
        for (uint xsi = 1; xsi <= current_max_xsi; xsi++) {
            bigint r_xsi = r[xsi];
            if (r_xsi == 0) continue;

            ThreadContext& ctx = *ctx_arr[omp_get_thread_num()];
            ctx.xs = xs_reverse_index[xsi];
            ctx.row_end = ctx.row;
            memset(ctx.counts, 0, sizeof(ctx.counts));
            
            next_h_rec(ctx);
            std::sort(ctx.row, ctx.row_end);

            ulong128 prev = ~0;
            uint weight = 0;
            for (ulong128* it = ctx.row; it <= ctx.row_end; it++) {
                bool end = (it == ctx.row_end);
                if (!end && *it == prev) {
                    weight++;
                } 
                else {
                    if (weight > 0) {
                        uint target_xsi = get_or_create_xsi(prev);
                        
                        // Extremely fine-grained GMP locking
                        omp_set_lock(&r_locks[target_xsi]);
                        r_[target_xsi] += (weight == 1) ? r_xsi : r_xsi * weight;
                        omp_unset_lock(&r_locks[target_xsi]);
                    }
                    if (!end) { prev = *it; weight = 1; }
                }
            }
            r[xsi] = 0; // Clear source

            if (xsi % 1000 == 0) {
                #pragma omp critical(print_lock)
                {
                    clock_t time_now = clock();
                    if (time_now > time_ + CLOCKS_PER_SEC / 5) {
                        std::cout << "Applying horizontal transition... " << xsi << " / " << current_max_xsi << "        \r" << std::flush;
                        time_ = time_now;
                    }
                }
            }
        }

        // We update current_max_xsi because next_h likely discovered new states
        current_max_xsi = global_xsi.load() - 1;
        
        // Output evaluation
        std::cout << std::setw(80) << "" << "\r" << std::flush;
        bigint v = 0;
        
        // Instead of keeping a separate list, just check the 5-bit condition natively
        #pragma omp parallel for reduction(+:v) schedule(static)
        for (uint xsi = 1; xsi <= current_max_xsi; xsi++) {
            if (((xs_reverse_index[xsi] >> (5 * (n - 1))) & 0x1f) == 0) {
                v += r_[xsi];
            }
        }
        
        std::cout << current_k << "\t" << v << "  [Active States: " << current_max_xsi << "]" << std::endl;
        if (output_file.is_open()) output_file << v << std::endl;

        if (current_k == max_k) break;

        // --- VERTICAL PASS ---
        #pragma omp parallel for schedule(dynamic, 32)
        for (uint xsi = 1; xsi <= current_max_xsi; xsi++) {
            bigint r_xsi = r_[xsi];
            if (r_xsi == 0) continue;

            ThreadContext& ctx = *ctx_arr[omp_get_thread_num()];
            ctx.xs = xs_reverse_index[xsi];
            ctx.row_end = ctx.row;
            memset(ctx.counts, 0, sizeof(ctx.counts));
            
            next_v_rec(ctx);
            std::sort(ctx.row, ctx.row_end);

            ulong128 prev = ~0;
            uint weight = 0;
            for (ulong128* it = ctx.row; it <= ctx.row_end; it++) {
                bool end = (it == ctx.row_end);
                if (!end && *it == prev) {
                    weight++;
                } 
                else {
                    if (weight > 0) {
                        uint target_xsi = get_or_create_xsi(prev);
                        
                        omp_set_lock(&r_locks[target_xsi]);
                        r[target_xsi] += (weight == 1) ? r_xsi : r_xsi * weight;
                        omp_unset_lock(&r_locks[target_xsi]);
                    }
                    if (!end) { prev = *it; weight = 1; }
                }
            }
            r_[xsi] = 0; 

            if (xsi % 1000 == 0) {
                #pragma omp critical(print_lock)
                {
                    clock_t time_now = clock();
                    if (time_now > time_ + CLOCKS_PER_SEC / 5) {
                        std::cout << "Applying vertical transition... " << xsi << " / " << current_max_xsi << "        \r" << std::flush;
                        time_ = time_now;
                    }
                }
            }
        }
    }

    // --- MEMORY CLEANUP ---
    for (size_t i = 0; i <= dim; i++) omp_destroy_lock(&r_locks[i]);
    for (int i = 0; i < NUM_SHARDS; i++) omp_destroy_lock(&shard_locks[i]);
    delete[] r_locks;
    delete[] r;
    delete[] r_;
    delete[] xs_reverse_index;
    for (int t = 0; t < max_threads; t++) delete ctx_arr[t];

    return 0;
}
