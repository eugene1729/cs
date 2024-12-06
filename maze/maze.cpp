#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <queue>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <gmpxx.h>

typedef mpz_class bigint;
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

std::vector<XSI> xsi_final;
std::vector<PAIR*> pairs[2];
XS_INDEX xs_index;
std::queue<ulong> queue;
ulong xs;
uchar counts[MAX_LAYER_SIZE + 1];
XSI* row;
XSI* row_end;
uint layer_size;
EDGE* h_edge_map;
uint edges_per_h_layer;

uint dimensions[3];
uint k;

void
state_transition (
    uint q,
    ulong xs_
    )

{
    XS_INDEX::iterator it = xs_index.find(xs_);
    XSI xsi;
    if (it == xs_index.end()) {
        queue.push(xs_);
        xsi = xs_index[xs_] = xs_index.size() + 1;
        if (((xs_ >> 4 * (layer_size - 1)) & 0xf) == 0) {
            xsi_final.push_back(xsi);
        }
        pairs[HORIZONTAL].push_back(NULL);
        pairs[VERTICAL].push_back(NULL);
    }
    else {
        xsi = it->second;
    }
    *row_end++ = xsi;
}

void
next_v (
    uint i = 0,
    uint q = 0
    )

{
    if (i == layer_size) {

        if (counts[0] == 0) {
            return;
        }

        uchar map[layer_size];
        for (uint j = 0; j < layer_size; j++) {
            map[j] = MAX;
        }

        uchar x_ = 1;
        ulong xs__ = 0;
        for (uint j = 0; j < layer_size; j++) {
            if ((q & 1 << j) != 0) {
                uchar x = (xs >> j * 4) & 0xf;
                if (x != 0) {
                    if (counts[x] == 1) {
                        x = x_++;
                    }
                    else {
                        if (map[x] == MAX) {
                            map[x] = x_++;
                        }
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

        state_transition(q, xs__);

        return;
    }

    next_v(i + 1, q);
    uint x = (xs >> i * 4) & 0xf;
    counts[x]++;
    next_v(i + 1, q | (1 << i));
    counts[x]--;
}

void
next_h (
    uint i = 0,
    uint q = 0
    )

{
    if (i == edges_per_h_layer) {

        uchar xs_[layer_size];
        for (uint j = 0; j < layer_size; j++) {
            xs_[j] = (xs >> j * 4) & 0xf;
            counts[j] = 0;
        }

        uchar map[layer_size];

        for (;;) {

            for (uint j = 0; j < layer_size; j++) {
                map[j] = MAX;
            }

            bool changed = false;

            for (uint j = 0; j < edges_per_h_layer; j++) {
                if ((q & 1 << j) != 0) {
                    EDGE* e = h_edge_map + j;
                    uchar x = xs_[e->i];
                    uchar x_ = xs_[e->j];
                    if (x < x_) {
                        if (map[x_] > x) {
                            map[x_] = x;
                            changed = true;
                        }
                    }
                    else if (x_ < x) {
                        if (map[x] > x_) {
                            map[x] = x_;
                            changed = true;
                        }
                    }
                }
            }

            if (!changed) {
                break;
            }

            for (uint j = 0; j < layer_size; j++) {
                uchar x = xs_[j];
                if (map[x] != MAX) {
                    do {
                        x = map[x];
                    } while (map[x] != MAX);
                    xs_[j] = x;
                }
            }
        }

        for (uint j = 0; j < layer_size; j++) {
            counts[xs_[j]]++;
        }

        if (counts[0] == 0) {
            return;
        }

        for (uint j = 0; j < layer_size; j++) {
            map[j] = MAX;
        }

        uchar x_ = 1;
        ulong xs__ = 0;
        for (uint j = 0; j < layer_size; j++) {
            uchar x = xs_[j];
            if (x != 0) {
                if (counts[x] == 1) {
                    x = x_++;
                }
                else {
                    if (map[x] == MAX) {
                        map[x] = x_++;
                    }
                    x = map[x];
                }
                xs__ |= (ulong) x << j * 4;
            }
        }

        for (uint j = 0; j < layer_size; j++) {
            counts[j] = 0;
        }

        state_transition(q, xs__);

        return;
    }

    next_h(i + 1, q);
    next_h(i + 1, q | (1 << i));
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
            for (uint i = 0; i < 3; i++) {
                dimensions[i] = 1;
            }
            for (uint i = 0; i < 3; i++) {
                char* token = strtok(i == 0 ? arg : NULL, "_.,:/xX");
                if (token == NULL) {
                    break;
                }
                dimensions[i] = atoi(token);
                if (dimensions[i] < 1 || dimensions[i] > MAX_LAYER_SIZE) {
                    std::cout << "ERROR: Dimension is out of range." << std::endl;
                    return 1;
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

    if (max_k == 0) {
        max_k = 1000;
    }

    layer_size = 1;
    for (uint i = 0; i < 3; i++) {
        layer_size *= dimensions[i];
    }
    if (layer_size > MAX_LAYER_SIZE) {
        std::cout << "ERROR: Dimensions are out of range." << std::endl;
        return 1;
    }

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
                if (i < dimensions[0] - 1) {
                    e->i = index(i, j, k);
                    e->j = index(i + 1, j, k);
                    e++;
                }
                if (j < dimensions[1] - 1) {
                    e->i = index(i, j, k);
                    e->j = index(i, j + 1, k);
                    e++;
                }
                if (k < dimensions[2] - 1) {
                    e->i = index(i, j, k);
                    e->j = index(i, j, k + 1);
                    e++;
                }
            }
        }
    }
    edges_per_h_layer = e - h_edge_map;

    for (uint type = HORIZONTAL; type <= VERTICAL; type++) {
        pairs[type].push_back(NULL);
        pairs[type].push_back(NULL);
        pairs[type].reserve(1 << 20);
    }
    xsi_final.reserve(1 << 20);
    if (layer_size == 1) {
        xsi_final.push_back(1);
    }
    row = new XSI[(1 << (edges_per_h_layer > layer_size ? edges_per_h_layer : layer_size)) + 1];
    xs = 0xfedcba9876543210ul & ((layer_size < 16 ? 1ul << layer_size * 4 : 0) - 1);
    xs_index[xs] = 1;
    queue.push(xs);
    ulong total_num_pairs = 0;

    clock_t time_ = 0;
    while (!queue.empty()) {

        xs = queue.front();
        queue.pop();
        XSI xsi = xs_index[xs];

        for (uint type = HORIZONTAL; type <= VERTICAL; type++) {

            row_end = row;
            type == HORIZONTAL ? next_h() : next_v();

            std::sort(row, row_end);
            uint num_pairs = 0;
            XSI xsi_ = 0;
            for (XSI* it = row; it < row_end; it++) {
                if (*it != xsi_) {
                    num_pairs++;
                    xsi_ = *it;
                }
            }
            *row_end++ = 0;
            total_num_pairs += num_pairs;

            PAIR* p = pairs[type][xsi] = new PAIR[num_pairs + 1];
            xsi_ = 0;
            ushort weight = 0;
            for (XSI* it = row; it < row_end; it++) {
                if (*it != xsi_) {
                    if (weight != 0) {
                        p->xsi = xsi_;
                        p->weight = weight;
                        p++;
                    }
                    xsi_ = *it;
                    weight = 0;
                }
                weight++;
            }
            p->xsi = 0;

            clock_t time = clock();
            if (time > time_ + CLOCKS_PER_SEC / 5) {
                std::cout << "Building transition matrices... " << xs_index.size() - queue.size() << "\r" << std::flush;
                time_ = time;
            }
        }
    }
    pairs[HORIZONTAL].shrink_to_fit();
    pairs[VERTICAL].shrink_to_fit();
    xsi_final.shrink_to_fit();
    uint dim = xs_index.size();
    std::cout << "Built " << dim << " × " << dim << " transition matrices. "
              << "Nonzero elements: " << total_num_pairs << "." << std::endl;

    std::ofstream output_file;
    if (output_filename != NULL) {
        output_file.open(output_filename, std::ofstream::out);
        if (!output_file.is_open()) {
            std::cout << "ERROR: Failed creating file '" << output_filename << "'." << std::endl;
            return 1;
        }
    }

    bigint* r = new bigint[dim + 1];
    bigint* r_ = new bigint[dim + 1];

    r[1] = 1;

    for (k = 1; ; k++) {

        for (XSI xsi = 1; xsi <= dim; xsi++) {
            bigint r_xsi = r[xsi];
            if (r_xsi != 0) {
                for (PAIR* p = pairs[HORIZONTAL][xsi]; p->xsi != 0; p++) {
                    r_[p->xsi] += p->weight != 1 ? p->weight * r_xsi : r_xsi;
                }
                r[xsi] = 0;
            }
        }

        bigint v = 0;
        for (std::vector<XSI>::iterator it = xsi_final.begin(); it < xsi_final.end(); it++) {
            v += r_[*it];
        }
        std::cout << k << "\t" << v << std::endl;
        if (output_file.is_open()) {
            output_file << v << std::endl;
        }
        if (k == max_k) {
            break;
        }

        for (XSI xsi = 1; xsi <= dim; xsi++) {
            bigint r_xsi = r_[xsi];
            if (r_xsi != 0) {
                for (PAIR* p = pairs[VERTICAL][xsi]; p->xsi != 0; p++) {
                    r[p->xsi] += p->weight != 1 ? p->weight * r_xsi : r_xsi;
                }
                r_[xsi] = 0;
            }
        }
    }

    return 0;
}
