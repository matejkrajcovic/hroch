#include"constants.h"
#include <iomanip>

string datapath("data/");

int LOWER_RANGE = 100;
int UPPER_RANGE = 199;
string TEST_CASE = "F4";

int num_bases = 4;
char bases[] = {'A', 'C', 'G', 'T'};
int base_id[256];
char base_inv[256];
int debuging = 0;
int stats = 0;
int do_cheeryness = 1;
int strict_compare = 0;
int fail_on_error = 1;
int error_happened = 0;
double epsilon = 1e-10;

void setup_constants() {
    random_init();
    For(i, num_bases) {
        base_id[(int)bases[i]] = i;
        base_inv[(int)bases[i]] = bases[num_bases-1-i];
    }
    cout << fixed << setprecision(6);
}

set<string> parse_arguments(int argc, char **argv) {
    do_cheeryness = 1;
    set<string> res;
    For(i, argc) {
        string arg = argv[i];
        if (arg == "--stats" || arg == "-s") stats = 1;
        if (arg == "--debug" || arg == "-d") debuging = 1;
        if (arg == "--no-cherry") do_cheeryness = 0;
        if (arg == "--strict") strict_compare = 1;
        if (arg == "--special-strict") strict_compare = SPECIAL_TRAINING;
        res.insert(arg);
    }
    return res;
}
