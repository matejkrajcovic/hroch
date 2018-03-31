#include <iostream>
#include <iomanip>

#include "constants.h"
#include "random.h"
#include "utils.h"
#include "history.h"
#include "../lib/cxxopts.hpp"

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

int strategy;
int reconstructions_count = -1;
string atoms_file = "";
string trees_dir = "";
string output_file_suffix = "";

bool no_annealing;
double starting_temperature;
int annealing_steps;
double prob_previously_used_event;

operation_mode parse_arguments(int argc, char **argv) {
    random_init();
    cout << fixed << setprecision(6);

    try {
        cxxopts::Options options(argv[0], "HROCH: Heuristic reconstruction of cluster histories.");
        options.add_options()
            ("help", "")
            ("s,stats", "")
            ("d,debug", "")
            ("no-cherry", "")
            ("strict", "")
            ("special-strict", "")
            ;

        options.add_options("Modes of operation")
            ("solve", "")
            ("gen-all", "generate test and train data")
            ("gen-test", "generate test data")
            ("train", "produce lr-train file")
            ("test-c", "make statistics of candidate proposal algorithms")
            ("test-s", "make statistics of scoring algorithms")
            ("rec", "make statistics of history reconstruction algorithms")
            ;

        options.add_options("Solve")
            ("atoms_file", "atom file", cxxopts::value<string>())
            ("trees_dir", "trees directory", cxxopts::value<string>())
            ("count", "count of reconstructed histories", cxxopts::value<int>(), "N")
            ("strategy", "reconstruction strategy", cxxopts::value<int>()->default_value(to_string(SCORE_LR)))
            ("output_file_suffix", "", cxxopts::value<string>()->default_value(""))
            ;

        options.add_options("Simulated annealing")
            ("no_annealing", "Disable simulated annealing", cxxopts::value<bool>()->default_value("false"))
            ("starting_temperature", "Starting temperature", cxxopts::value<double>()->default_value("0.2"))
            ("annealing_steps", "Annealing steps", cxxopts::value<int>()->default_value("10"))
            ("prob_previously_used_event", "Minimum probability of using an event from previous reconstruction", cxxopts::value<double>()->default_value("0"))
            ;

        auto results = options.parse(argc, argv);

        if (results.count("stats")) {
            stats = 1;
        }
        if (results.count("debug")) {
            debugging = 1;
        }
        if (results.count("no-cherry")) {
            do_cherryness = 0;
        }
        if (results.count("strict")) {
            strict_compare = 1;
        }
        if (results.count("special-strict")) {
            strict_compare = SPECIAL_TRAINING;
        }

        if (results.count("atoms_file")) {
            atoms_file = results["atoms_file"].as<string>();
        }
        if (results.count("trees_dir")) {
            trees_dir = results["trees_dir"].as<string>();
        }
        if (results.count("count")) {
            reconstructions_count = results["count"].as<int>();
        }
        strategy = results["strategy"].as<int>();
        output_file_suffix = results["output_file_suffix"].as<string>();

        no_annealing = results["no_annealing"].as<bool>();
        starting_temperature = results["starting_temperature"].as<double>();
        annealing_steps = results["annealing_steps"].as<int>();
        prob_previously_used_event = results["prob_previously_used_event"].as<double>();

        if (results.count("help")) {
            cout << options.help({"", "Modes of operation", "Solve", "Simulated annealing"}) << endl;
            exit(0);
        } else if (results.count("solve")) {
            if (atoms_file.empty() || trees_dir.empty() || reconstructions_count < 0) {
                cout << "Missing arguments: atoms_file, trees_dir or count." << endl;
                exit(0);
            }
            return operation_mode::solve;
        } else if (results.count("gen-all")) {
            return operation_mode::gen_all;
        } else if (results.count("gen-test")) {
            return operation_mode::gen_test;
        } else if (results.count("train")) {
            return operation_mode::train;
        } else if (results.count("test-c")) {
            return operation_mode::test_c;
        } else if (results.count("test-s")) {
            return operation_mode::test_s;
        } else {
            cout << "no mode of operation, see help" << endl;
            exit(1);
        }
    } catch (const cxxopts::OptionException& e) {
        cout << "error parsing options: " << e.what() << endl;
        exit(1);
    }
}
