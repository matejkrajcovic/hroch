#include <iostream>
#include <iomanip>

#include "constants.h"
#include "random.h"
#include "utils.h"

using namespace std;

string datapath("data/");

int LOWER_RANGE = 100;
int UPPER_RANGE = 199;
string TEST_CASE = "F4";

int num_bases = 4;
char bases[] = {'A', 'C', 'G', 'T'};
int debugging = 0;
int stats = 0;
int do_cherryness = 1;
int strict_compare = 0;
int fail_on_error = 1;
int error_happened = 0;
double epsilon = 1e-10;

void setup_constants() {
    random_init();
    cout << fixed << setprecision(6);
}

set<string> parse_arguments(int argc, char **argv) {
    do_cherryness = 1;
    set<string> res;
    For(i, argc) {
        string arg = argv[i];
        if (arg == "--stats" || arg == "-s") stats = 1;
        if (arg == "--debug" || arg == "-d") debugging = 1;
        if (arg == "--no-cherry") do_cherryness = 0;
        if (arg == "--strict") strict_compare = 1;
        if (arg == "--special-strict") strict_compare = SPECIAL_TRAINING;
        res.insert(arg);
    }
    return res;
}
