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
annealing_schedule_enum annealing_schedule;

int gen_count;
string gen_prefix;
double gen_time;

string reconstructions_file;

neighbor_selection_enum neighbor_selection;

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
            ("gen", "generate custom training data")
            ("gen-all", "generate test and train data")
            ("gen-test", "generate test data")
            ("train", "produce lr-train file")
            ("test-c", "make statistics of candidate proposal algorithms")
            ("test-s", "make statistics of scoring algorithms")
            ("rec", "make statistics of history reconstruction algorithms")
            ("evaluate-reconstructions", "compare correct history and reconstructions")
            ;

        options.add_options("Solve")
            ("atoms_file", "atom file", cxxopts::value<string>())
            ("trees_dir", "trees directory", cxxopts::value<string>())
            ("count", "count of reconstructed histories", cxxopts::value<int>(), "N")
            ("strategy", "reconstruction strategy", cxxopts::value<int>()->default_value(to_string(SCORE_LR)))
            ("output_file_suffix", "", cxxopts::value<string>()->default_value(""))
            ("reconstructions_file", "file with reconstructions to compare", cxxopts::value<string>())
            ;

        options.add_options("Simulated annealing")
            ("no_annealing", "Disable simulated annealing", cxxopts::value<bool>()->default_value("false"))
            ("starting_temperature", "Starting temperature", cxxopts::value<double>()->default_value("0.2"))
            ("annealing_steps", "Annealing steps", cxxopts::value<int>()->default_value("10"))
            ("prob_previously_used_event", "Minimum probability of using an event from previous reconstruction", cxxopts::value<double>()->default_value("0"))
            ("annealing_schedule", "", cxxopts::value<string>()->default_value("advanced"))
            ("neighbor_selection", "method to select neighbors", cxxopts::value<string>()->default_value("none"))
            ;

        options.add_options("Generate histories")
            ("gen_count", "Count of generated histories", cxxopts::value<int>()->default_value("100"))
            ("gen_prefix", "Prefix of generated histories", cxxopts::value<string>()->default_value("F2"))
            ("gen_time", "Time to simulate events", cxxopts::value<double>()->default_value("0.04"))
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
        if (results.count("reconstructions_file")) {
            reconstructions_file = results["reconstructions_file"].as<string>();
        }
        strategy = results["strategy"].as<int>();
        output_file_suffix = results["output_file_suffix"].as<string>();

        no_annealing = results["no_annealing"].as<bool>();
        starting_temperature = results["starting_temperature"].as<double>();
        annealing_steps = results["annealing_steps"].as<int>();
        prob_previously_used_event = results["prob_previously_used_event"].as<double>();
        string schedule_string = results["annealing_schedule"].as<string>();
        string neighbor_selection_string = results["neighbor_selection"].as<string>();

        if (schedule_string == "simple") {
            annealing_schedule = annealing_schedule_enum::simple;
        } else if (schedule_string == "advanced") {
            annealing_schedule = annealing_schedule_enum::advanced;
        } else if (schedule_string == "baseline_advanced") {
            annealing_schedule = annealing_schedule_enum::baseline_advanced;
        } else {
            cerr << "Wrong argument to parameter \"annealing_schedule\" " << schedule_string << endl;
            exit(1);
        }

        if (neighbor_selection_string == "none") {
            neighbor_selection = neighbor_selection_enum::none;
        } else if (neighbor_selection_string == "prioritize_used_events") {
            neighbor_selection = neighbor_selection_enum::prioritize_used_events;
        } else if (neighbor_selection_string == "change_event_simple") {
            neighbor_selection = neighbor_selection_enum::change_event_simple;
        } else {
            cerr << "Wrong argument to parameter \"neighbor_selection\" " << neighbor_selection_string << endl;
            exit(1);
        }

        if (results.count("gen_count")) {
            gen_count = results["gen_count"].as<int>();
        }
        if (results.count("gen_prefix")) {
            gen_prefix = results["gen_prefix"].as<string>();
        }
        if (results.count("gen_time")) {
            gen_time = results["gen_time"].as<double>();
        }

        if (results.count("help")) {
            cout << options.help({"", "Modes of operation", "Solve", "Simulated annealing", "Generate histories"}) << endl;
            exit(0);
        } else if (results.count("solve")) {
            if (atoms_file.empty() || trees_dir.empty() || reconstructions_count < 0) {
                cout << "Missing arguments: atoms_file, trees_dir or count." << endl;
                exit(0);
            }
            return operation_mode::solve;
        } else if (results.count("gen")) {
            return operation_mode::gen;
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
        } else if (results.count("rec")) {
            return operation_mode::rec;
        } else if (results.count("evaluate-reconstructions")) {
            return operation_mode::evaluate_reconstruction;
        } else {
            cout << "no mode of operation, see help" << endl;
            exit(1);
        }
    } catch (const cxxopts::OptionException& e) {
        cout << "error parsing options: " << e.what() << endl;
        exit(1);
    }
}
