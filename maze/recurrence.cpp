#include <iostream>
#include <fstream>
#include <iomanip>
#include <gmpxx.h>

typedef mpz_class bigint;

void
print_m (
    int n,
    bigint** m
    )

{
    for (int j = 0; j < n; j++) {
        for (int i = 0; i <= n; i++) {
            std::cout << m[j][i] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int
main (
    int argc,
    char** argv
    )

{
    if (argc < 3) {
        std::cout << "ERROR: Required arguments are missing." << std::endl;
        return 1;
    }

    int n = atoi(argv[1]);
    char* input_filename = argv[2];
    char* output_filename = NULL;

    for (int i = 3; i < argc; i++) {
        char* arg = argv[i];
        if (arg[0] == '-') {
            if (i < argc - 1 && strcmp(arg, "-o") == 0) {
                output_filename = argv[++i];
            }
            else {
                std::cout << "ERROR: Unexpected argument." << std::endl;
                return 1;
            }
        }
        else if (n == 0) {
            std::cout << "ERROR: Unexpected argument." << std::endl;
            return 1;
        }
    }

    std::ifstream input_file(input_filename);
    if (!input_file.is_open()) {
        std::cout << "ERROR: Failed opening file '" << input_filename << "'." << std::endl;
        return 1;
    }

    bigint** m = new bigint* [n];
    for (int i = 0; i < n; i++) {
        m[i] = new bigint [n + 1];
    }

    for (int i = 0; i < 2 * n; i++) {
        bigint x;
        if (!(input_file >> x)) {
            std::cout << "ERROR: File '" << input_filename << "' does not have enough data." << std::endl;
            return 1;
        }
        for (int j = i - n > 0 ? i - n : 0; j < (n < i + 1 ? n : i + 1); j++) {
            m[j][i - j] = x;
        }
    }

    input_file.close();

    int n_ = 0;
    for (int i = 0; i < n - 1; i++) {
        std::cout << "Transforming into row echelon form... " << i << "\r" << std::flush;
        bool all_zeros = false;
        for (int j = i + 1; j < n; j++) {
            if (m[i][i] != 0) {
                bigint q1 = m[j][i];
                bigint q2 = m[i][i];
                bigint g = 0;
                for (int k = i; k <= n; k++) {
                    m[j][k] = m[j][k] * q2 - q1 * m[i][k];
                    g = gcd(g, m[j][k]);
                }
                if (g == 0) {
                    all_zeros = true;
                    break;
                }
                if (g != 1) {
                    for (int k = 0; k <= n; k++) {
                        m[j][k] /= g;
                    }
                }
            }
        }
        if (all_zeros) {
            n_ = i + 1;
            break;
        }
    }
    std::cout << std::setw(80) << "" << "\r" << std::flush;

    if (n_ == 0) {
        std::cout << "Could not find a recurrence." << std::endl;
        return 0;
    }
    n = n_;

    std::cout << "Found recurrence of order " << n << "." << std::endl;
    for (int i = n - 1; i >= 0; i--) {
        std::cout << "Substituting... " << i << "\r" << std::flush;
        bigint q = m[i][i];
        for (int k = i; k <= n; k++) {
            m[i][k] /= q;
        }
        for (int j = i - 1; j >= 0; j--) {
            bigint q1 = m[j][i];
            bigint q2 = m[i][i];
            for (int k = j; k <= n; k++) {
                m[j][k] = m[j][k] * q2 - q1 * m[i][k];
            }
        }
    }
    std::cout << std::setw(80) << "" << "\r" << std::flush;

    std::ofstream output_file;
    if (output_filename != NULL) {
        output_file.open(output_filename, std::ofstream::out);
        if (!output_file.is_open()) {
            std::cout << "ERROR: Failed creating file '" << output_filename << "'." << std::endl;
            return 1;
        }
    }

    std::cout << "a(n) = ";
    if (output_file.is_open()) {
        output_file << "a(n) = ";
    }
    for (int j = 0; j < n; j++) {
        bigint q = m[n - j - 1][n];
        if (j != 0) {
            std::cout << std::setw(7) << (q < 0 ? "- " : "+ ");
            if (output_file.is_open()) {
                output_file << std::setw(7) << (q < 0 ? "- " : "+ ");
            }
            if (q < 0) {
                q = -q;
            }
        }
        std::cout << q << " * a(n-" << (j + 1) << ")" << std::endl;
        if (output_file.is_open()) {
            output_file << q << " * a(n-" << (j + 1) << ")" << std::endl;
        }
    }

    return 0;
}
