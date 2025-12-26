#include <iostream>
#include <string>
#include <cstdlib>   // std::strtoul
#include <cerrno>    // errno
#include <climits>   // ULONG_MAX

#include "gspan.h"  

static void usage() {
    std::cerr
        << "Usage: gspan_cli [options] < input.data > output.txt\n"
        << "Options:\n"
        << "  -m <minsup>    Minimum support (same as -s)\n"
        << "  -s <minsup>    Minimum support (same as -m)\n"
        << "  -n <minnodes>  Minimum number of nodes in a pattern\n"
        << "  -L <maxpat>    Maximum pattern size (edges) / limit (as in original)\n"
        << "  -d             Enable encoding (same as -e)\n"
        << "  -e             Enable encoding (same as -d)\n"
        << "  -w             Output where (occurrences)\n"
        << "  -D             Directed graphs\n"
        << "  -h             Show this help\n";
}

// Strict-ish unsigned int parser with basic validation.
static bool parse_uint(const char* s, unsigned int& out) {
    if (!s || !*s) return false;
    errno = 0;
    char* end = nullptr;
    unsigned long v = std::strtoul(s, &end, 10);
    if (errno != 0 || end == s || *end != '\0') return false;
    if (v > static_cast<unsigned long>(UINT_MAX)) return false;
    out = static_cast<unsigned int>(v);
    return true;
}

int main(int argc, char** argv) {
    unsigned int minsup   = 1;
    unsigned int maxpat   = 0xffffffffu;
    unsigned int minnodes = 0;
    bool where    = false;
    bool enc      = false;
    bool directed = false;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];

        auto need_value = [&](const char* flag) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value after " << flag << "\n";
                usage();
                std::exit(-1);
            }
            return argv[++i];
        };

        if (a == "-h") {
            usage();
            return 0;
        } else if (a == "-m" || a == "-s") {
            unsigned int v = 0;
            const char* s = need_value(a.c_str());
            if (!parse_uint(s, v)) {
                std::cerr << "Invalid value for " << a << ": " << s << "\n";
                usage();
                return -1;
            }
            minsup = v;
        } else if (a == "-n") {
            unsigned int v = 0;
            const char* s = need_value("-n");
            if (!parse_uint(s, v)) {
                std::cerr << "Invalid value for -n: " << s << "\n";
                usage();
                return -1;
            }
            minnodes = v;
        } else if (a == "-L") {
            unsigned int v = 0;
            const char* s = need_value("-L");
            if (!parse_uint(s, v)) {
                std::cerr << "Invalid value for -L: " << s << "\n";
                usage();
                return -1;
            }
            maxpat = v;
        } else if (a == "-d" || a == "-e") {
            enc = true;
        } else if (a == "-w") {
            where = true;
        } else if (a == "-D") {
            directed = true;
        } else {
            std::cerr << "Unknown option: " << a << "\n";
            usage();
            return -1;
        }
    }

    // Optional I/O speedups (safe and portable)
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    GSPAN::gSpan gspan;
    gspan.run(std::cin, std::cout, minsup,minnodes, maxpat, enc, where, directed);
    return 0;
}
