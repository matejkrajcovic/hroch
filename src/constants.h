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
extern int debugging;
extern int stats;
extern int do_cherryness;
#define SPECIAL_TRAINING 47
extern int strict_compare;
extern int num_bases;
extern char bases[];
extern double epsilon;

extern int strategy;
extern int reconstructions_count;
extern std::string atoms_file;
extern std::string trees_dir;
extern std::string output_file_suffix;

extern bool no_annealing;
extern double starting_temperature;
extern int annealing_steps;
extern double prob_previously_used_event;

extern int gen_count;
extern std::string gen_prefix;
extern double gen_time;

enum class annealing_schedule_enum {
  simple,
  advanced,
  baseline_advanced,
};

extern annealing_schedule_enum annealing_schedule;

enum class operation_mode {
  solve,
  gen,
  gen_all,
  gen_test,
  train,
  test_c,
  test_s,
  rec,
};

operation_mode parse_arguments(int argc, char **argv);

#endif
