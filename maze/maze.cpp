#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <queue>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <algorithm>
#include <omp.h>
#include <gmpxx.h>

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long ulong;

const uchar MAX_LAYER_SIZE = 16;
const uchar MAX = (uchar) -1;

const uint HORIZONTAL = 0;
const uint VERTICAL = 1;

typedef uint XSI;
typedef std::map<ulong, XSI> XS_INDEX;

struct __attribute__((packed)) PAIR {
    XSI xsi;
    ushort weight;
};

struct EDGE {
    uchar i;
    uchar j;
};

struct RawTransition {
    ulong target_xs;
    ushort weight;
};

// Thread-local context
struct ThreadContext {
    ulong xs;
    uchar counts[MAX_LAYER_SIZE + 1];
    ulong* row_ulong; 
    ulong* row_ulong_end;
    
    std::vector<ulong> local_new_states;
    std::vector<std::pair<XSI, std::vector<RawTransition>>> local_h_trans;
    std::vector<std::pair<XSI, std::vector<RawTransition>>> local_v_trans;
};

std::vector<XSI> xsi_final;
std::vector<PAIR*> pairs[2];
XS_INDEX xs_index;
uint layer_size;
EDGE* h_edge_map;
uint edges_per_h_layer;

uint dimensions[3];
uint k;

// Simplified: just collects the raw state targets for the thread
void 
state_transition (
    ThreadContext& ctx, 
    uint q, 
    ulong xs_
    ) 

{
    *ctx.row_ulong_end++ = xs_;
}

void 
next_v (
    ThreadContext& ctx, 
    uint i = 0, 
    uint q = 0
    ) 

{
    if (i == layer_size) {
        if (ctx.counts[0] == 0) return;
        uchar map[MAX_LAYER_SIZE];
        for (uint j = 0; j < layer_size; j++) map[j] = MAX;

        uchar x_ = 1;
        ulong xs__ = 0;
        for (uint j = 0; j < layer_size; j++) {
            if ((q & 1 << j) != 0) {
                uchar x = (ctx.xs >> j * 4) & 0xf;
                if (x != 0) {
                    if (ctx.counts[x] == 1) x = x_++;
                    else {
                        if (map[x] == MAX) map[x] = x_++;
                        x = map[x];
                    }
                    xs__ |= (ulong) x << j * 4;
                }
            } 
            else {
                xs__ |= (ulong) x_ << j * 4;
                x_++;
            }
        }
        state_transition(ctx, q, xs__);
        return;
    }
    next_v(ctx, i + 1, q);
    uint x = (ctx.xs >> i * 4) & 0xf;
    ctx.counts[x]++;
    next_v(ctx, i + 1, q | (1 << i));
    ctx.counts[x]--;
}

void 
next_h (
    ThreadContext& ctx, 
    uint i = 0, 
    uint q = 0
    ) 

{
    if (i == edges_per_h_layer) {
        uchar xs_[MAX_LAYER_SIZE];
        for (uint j = 0; j < layer_size; j++) {
            xs_[j] = (ctx.xs >> j * 4) & 0xf;
            ctx.counts[j] = 0;
        }

        uchar map[MAX_LAYER_SIZE];
        for (;;) {
            for (uint j = 0; j < layer_size; j++) map[j] = MAX;
            bool changed = false;
            for (uint j = 0; j < edges_per_h_layer; j++) {
                if ((q & 1 << j) != 0) {
                    EDGE* e = h_edge_map + j;
                    uchar x = xs_[e->i];
                    uchar x_ = xs_[e->j];
                    if (x < x_) {
                        if (map[x_] > x) { map[x_] = x; changed = true; }
                    } 
                    else if (x_ < x) {
                        if (map[x] > x_) { map[x] = x_; changed = true; }
                    }
                }
            }
            if (!changed) break;
            for (uint j = 0; j < layer_size; j++) {
                uchar x = xs_[j];
                if (map[x] != MAX) {
                    do { x = map[x]; } while (map[x] != MAX);
                    xs_[j] = x;
                }
            }
        }

        for (uint j = 0; j < layer_size; j++) ctx.counts[xs_[j]]++;
        if (ctx.counts[0] == 0) return;
        for (uint j = 0; j < layer_size; j++) map[j] = MAX;

        uchar x_ = 1;
        ulong xs__ = 0;
        for (uint j = 0; j < layer_size; j++) {
            uchar x = xs_[j];
            if (x != 0) {
                if (ctx.counts[x] == 1) x = x_++;
                else {
                    if (map[x] == MAX) map[x] = x_++;
                    x = map[x];
                }
                xs__ |= (ulong) x << j * 4;
            }
        }
        for (uint j = 0; j < layer_size; j++) ctx.counts[j] = 0;
        state_transition(ctx, q, xs__);
        return;
    }
    next_h(ctx, i + 1, q);
    next_h(ctx, i + 1, q | (1 << i));
}

inline 
uint 
index (
    uint i, 
    uint j, 
    uint k
    ) 

{
    return i * dimensions[1] * dimensions[2] + j * dimensions[2] + k;
}

int 
main (
    int argc, 
    char* argv[]
    ) 

{
    uint max_k = 0;
    char* output_filename = NULL;

    for (int i = 1; i < argc; i++) {
        char* arg = argv[i];
        if (arg[0] == '-') {
            if (i < argc - 1 && strcmp(arg, "-o") == 0) {
                output_filename = argv[++i];
            } 
            else {
                std::cout << "ERROR: Unexpected argument '" << arg << "'." << std::endl;
                return 1;
            }
        } 
        else if (dimensions[0] == 0) {
            for (uint j = 0; j < 3; j++) dimensions[j] = 1;
            for (uint j = 0; j <= 3; j++) {
                char* token = strtok(j == 0 ? arg : NULL, "_.,:/xX");
                if (token == NULL) break;
                if (j == 3) { std::cout << "ERROR: Too many dimensions." << std::endl; return 1; }
                dimensions[j] = atoi(token);
                if (dimensions[j] < 1 || dimensions[j] > MAX_LAYER_SIZE) {
                    std::cout << "ERROR: Dimension is out of range." << std::endl; return 1;
                }
            }
        } 
        else if (max_k == 0) {
            max_k = atoi(arg);
        } 
        else {
            std::cout << "ERROR: Unexpected argument '" << arg << "'." << std::endl;
            return 1;
        }
    }

    if (dimensions[0] == 0) { std::cout << "ERROR: Dimensions are required." << std::endl; return 1; }
    if (max_k == 0) max_k = 1000;

    layer_size = 1;
    for (uint i = 0; i < 3; i++) layer_size *= dimensions[i];
    if (layer_size > MAX_LAYER_SIZE) { std::cout << "ERROR: Dimensions out of range." << std::endl; return 1; }

    bool all_ones = true;
    std::cout << "Calculating for ";
    for (uint i = 0; i < 3; i++) {
        if (dimensions[i] != 1) {
            std::cout << dimensions[i] << " × ";
            all_ones = false;
        }
    }
    if (all_ones) {
        std::cout << "1";
    }
    std::cout << "K with K ∈ [1 .. " << max_k << "]." << std::endl;

    h_edge_map = new EDGE[3 * layer_size];
    EDGE* e = h_edge_map + 0;
    for (uint i = 0; i < dimensions[0]; i++) {
        for (uint j = 0; j < dimensions[1]; j++) {
            for (uint k = 0; k < dimensions[2]; k++) {
                if (i < dimensions[0] - 1) { e->i = index(i, j, k); e->j = index(i + 1, j, k); e++; }
                if (j < dimensions[1] - 1) { e->i = index(i, j, k); e->j = index(i, j + 1, k); e++; }
                if (k < dimensions[2] - 1) { e->i = index(i, j, k); e->j = index(i, j, k + 1); e++; }
            }
        }
    }
    edges_per_h_layer = e - h_edge_map;

    for (uint type = HORIZONTAL; type <= VERTICAL; type++) {
        pairs[type].push_back(NULL); pairs[type].push_back(NULL);
        pairs[type].reserve(1 << 24); 
    }
    xsi_final.reserve(1 << 24);
    if (layer_size == 1) xsi_final.push_back(1);
    
    ulong initial_xs = 0xfedcba9876543210ul & ((layer_size < 16 ? 1ul << layer_size * 4 : 0) - 1);
    xs_index[initial_xs] = 1;
    
    std::vector<ulong> current_frontier;
    std::vector<ulong> next_frontier;
    current_frontier.push_back(initial_xs);
    
    ulong total_num_pairs = 0;
    clock_t time_ = 0;

    // --- PHASE 1: LOCK-FREE BFS ---
    while (!current_frontier.empty()) {
        next_frontier.clear();
        ulong local_total_pairs = 0;

        #pragma omp parallel
        {
            ThreadContext ctx;
            ctx.row_ulong = new ulong[(1 << (edges_per_h_layer > layer_size ? edges_per_h_layer : layer_size)) + 1];
            memset(ctx.counts, 0, sizeof(ctx.counts));

            // PASS A: Process chunk, store raw targets privately
            #pragma omp for schedule(dynamic)
            for (size_t i = 0; i < current_frontier.size(); i++) {
                ctx.xs = current_frontier[i];
                XSI source_xsi = xs_index[ctx.xs]; // Safe read

                for (uint type = HORIZONTAL; type <= VERTICAL; type++) {
                    ctx.row_ulong_end = ctx.row_ulong;
                    type == HORIZONTAL ? next_h(ctx) : next_v(ctx);

                    std::sort(ctx.row_ulong, ctx.row_ulong_end);
                    
                    std::vector<RawTransition> trans;
                    ulong current_xs_ = ~0ul; // Guarantees mismatch on first element
                    ushort weight = 0;
                    
                    for (ulong* it = ctx.row_ulong; it < ctx.row_ulong_end; it++) {
                        if (*it != current_xs_) {
                            if (weight != 0) {
                                trans.push_back({current_xs_, weight});
                                if (xs_index.find(current_xs_) == xs_index.end()) {
                                    ctx.local_new_states.push_back(current_xs_);
                                }
                            }
                            current_xs_ = *it;
                            weight = 0;
                        }
                        weight++;
                    }
                    if (weight != 0) {
                        trans.push_back({current_xs_, weight});
                        if (xs_index.find(current_xs_) == xs_index.end()) {
                            ctx.local_new_states.push_back(current_xs_);
                        }
                    }

                    if (type == HORIZONTAL) ctx.local_h_trans.push_back({source_xsi, trans});
                    else ctx.local_v_trans.push_back({source_xsi, trans});
                }
            }

            // PASS B: Only lock for a fraction of a millisecond to issue official IDs
            #pragma omp critical(state_management)
            {
                for (ulong new_xs : ctx.local_new_states) {
                    if (xs_index.find(new_xs) == xs_index.end()) {
                        XSI xsi = xs_index.size() + 1;
                        xs_index[new_xs] = xsi;
                        next_frontier.push_back(new_xs);
                        
                        if (((new_xs >> 4 * (layer_size - 1)) & 0xf) == 0) {
                            xsi_final.push_back(xsi);
                        }
                        pairs[HORIZONTAL].push_back(NULL);
                        pairs[VERTICAL].push_back(NULL);
                    }
                }
            }

            // Wait for all IDs to be assigned
            #pragma omp barrier 

            // PASS C: Completely lock-free array building
            ulong thread_pair_count = 0;
            
            for (auto& item : ctx.local_h_trans) {
                XSI src = item.first;
                PAIR* p = new PAIR[item.second.size() + 1];
                for (size_t j = 0; j < item.second.size(); j++) {
                    p[j].xsi = xs_index[item.second[j].target_xs]; // safe read
                    p[j].weight = item.second[j].weight;
                }
                p[item.second.size()].xsi = 0;
                pairs[HORIZONTAL][src] = p; // safe write, src is unique
                thread_pair_count += item.second.size();
            }
            
            for (auto& item : ctx.local_v_trans) {
                XSI src = item.first;
                PAIR* p = new PAIR[item.second.size() + 1];
                for (size_t j = 0; j < item.second.size(); j++) {
                    p[j].xsi = xs_index[item.second[j].target_xs]; // safe read
                    p[j].weight = item.second[j].weight;
                }
                p[item.second.size()].xsi = 0;
                pairs[VERTICAL][src] = p;
                thread_pair_count += item.second.size();
            }

            #pragma omp atomic
            local_total_pairs += thread_pair_count;

            delete[] ctx.row_ulong;
        }

        current_frontier = next_frontier;
        total_num_pairs += local_total_pairs;

        clock_t time = clock();
        if (time > time_ + CLOCKS_PER_SEC / 5) {
            std::cout << "Building transfer matrices. Processed depth layer of size " << current_frontier.size() << "." << std::setw(80) << "" << "\r" << std::flush;
            time_ = time;
        }
    }

    pairs[HORIZONTAL].shrink_to_fit();
    pairs[VERTICAL].shrink_to_fit();
    xsi_final.shrink_to_fit();
    uint dim = xs_index.size();
    std::cout << "Built " << dim << " × " << dim << " transfer matrices. Nonzeros: " << total_num_pairs << "." << std::setw(80) << "" << std::endl;

    // --- PHASE 2: GRAPH TRANSPOSITION ---
    // We reverse the edges so threads can pull data lock-free
    std::vector<PAIR*> in_pairs_h(dim + 1, nullptr);
    std::vector<PAIR*> in_pairs_v(dim + 1, nullptr);
    std::vector<uint> in_degree_h(dim + 1, 0);
    std::vector<uint> in_degree_v(dim + 1, 0);

    for (XSI xsi = 1; xsi <= dim; xsi++) {
        if (pairs[HORIZONTAL][xsi]) for (PAIR* p = pairs[HORIZONTAL][xsi]; p->xsi != 0; p++) in_degree_h[p->xsi]++;
        if (pairs[VERTICAL][xsi]) for (PAIR* p = pairs[VERTICAL][xsi]; p->xsi != 0; p++) in_degree_v[p->xsi]++;
    }

    for (XSI target = 1; target <= dim; target++) {
        in_pairs_h[target] = new PAIR[in_degree_h[target] + 1];
        in_pairs_v[target] = new PAIR[in_degree_v[target] + 1];
    }

    std::vector<uint> current_idx_h(dim + 1, 0);
    std::vector<uint> current_idx_v(dim + 1, 0);

    for (XSI xsi = 1; xsi <= dim; xsi++) {
        if (pairs[HORIZONTAL][xsi]) {
            for (PAIR* p = pairs[HORIZONTAL][xsi]; p->xsi != 0; p++) {
                XSI target = p->xsi;
                uint idx = current_idx_h[target]++;
                in_pairs_h[target][idx] = {xsi, p->weight};
            }
        }
        if (pairs[VERTICAL][xsi]) {
            for (PAIR* p = pairs[VERTICAL][xsi]; p->xsi != 0; p++) {
                XSI target = p->xsi;
                uint idx = current_idx_v[target]++;
                in_pairs_v[target][idx] = {xsi, p->weight};
            }
        }
    }

    for (XSI target = 1; target <= dim; target++) {
        in_pairs_h[target][in_degree_h[target]].xsi = 0;
        in_pairs_v[target][in_degree_v[target]].xsi = 0;
    }

    // --- PHASE 3: LOCK-FREE DYNAMIC PROGRAMMING ---
    std::ofstream output_file;
    if (output_filename != NULL) output_file.open(output_filename, std::ofstream::out);

    mpz_class* r = new mpz_class[dim + 1];
    mpz_class* r_ = new mpz_class[dim + 1];
    for (XSI i = 0; i <= dim; i++) { r[i] = 0; r_[i] = 0; }
    r[1] = 1;

    for (k = 1; ; k++) {
        
        // Target nodes PULL from source nodes. Perfectly safe to parallelize!
        #pragma omp parallel for schedule(static)
        for (XSI target = 1; target <= dim; target++) {
            mpz_class sum = 0;
            for (PAIR* p = in_pairs_h[target]; p->xsi != 0; p++) {
                sum += p->weight != 1 ? p->weight * r[p->xsi] : r[p->xsi];
            }
            r_[target] = sum;
        }

        // Clear r for the next vertical pull
        #pragma omp parallel for schedule(static)
        for (XSI target = 1; target <= dim; target++) r[target] = 0;

        mpz_class v = 0;
        for (XSI final_xsi : xsi_final) v += r_[final_xsi]; // Sequential reduction
        
        std::cout << k << "\t" << v << std::endl;
        if (output_file.is_open()) output_file << v << std::endl;
        if (k == max_k) break;

        #pragma omp parallel for schedule(static)
        for (XSI target = 1; target <= dim; target++) {
            mpz_class sum = 0;
            for (PAIR* p = in_pairs_v[target]; p->xsi != 0; p++) {
                sum += p->weight != 1 ? p->weight * r_[p->xsi] : r_[p->xsi];
            }
            r[target] = sum;
        }

        #pragma omp parallel for schedule(static)
        for (XSI target = 1; target <= dim; target++) r_[target] = 0;
    }

    // --- MEMORY CLEANUP ---
    delete[] h_edge_map;
    delete[] r;
    delete[] r_;

    for (uint type = HORIZONTAL; type <= VERTICAL; type++) {
        for (size_t i = 0; i < pairs[type].size(); i++) {
            if (pairs[type][i] != NULL) delete[] pairs[type][i];
        }
    }
    for (XSI target = 1; target <= dim; target++) {
        delete[] in_pairs_h[target];
        delete[] in_pairs_v[target];
    }

    return 0;
}
