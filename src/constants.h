#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>
#include <vector>
#include <set>
#include <tuple>

typedef std::pair<int, int> pii;

typedef std::vector<int> vi;
typedef std::vector<double> vdo;
typedef std::vector<vdo> vvdo;

extern std::string datapath;

extern int LOWER_RANGE;
extern int UPPER_RANGE;
extern std::string TEST_CASE;

extern int error_happened;
extern int debuging;
extern int stats;
extern int do_cheeryness;
#define SPECIAL_TRAINING 47
extern int strict_compare;
extern int num_bases;
extern char bases[];
extern int base_id[256];
extern char base_inv[256];
extern double epsilon;

void setup_constants();
std::set<std::string> parse_arguments(int argc, char **argv);

#endif
