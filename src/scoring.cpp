#include <limits>
#include "scoring.h"

using namespace std;

double get_score(History* history, string atoms_filename, string align_dir) {
  if (scoring == scoring_enum::num_events) {
    return history->get_history_score_num_events();
  } else if (scoring == scoring_enum::likelihood) {
    return history->get_history_score_likelihood(atoms_filename, align_dir);
  } else {
    return 0;
  }
}

double get_initial_score() {
  if (scoring == scoring_enum::num_events) {
    return numeric_limits<int>::max();
  } else if (scoring == scoring_enum::likelihood) {
    return -numeric_limits<double>::max();
  } else {
    return 0;
  }
}

double get_negative_difference(double previous_score, double current_score) {
  if (scoring == scoring_enum::num_events) {
    return previous_score - current_score;
  } else if (scoring == scoring_enum::likelihood) {
    return current_score - previous_score;
  } else {
    return 0;
  }
}

bool is_better_score(double previous_score, double current_score) {
  if (scoring == scoring_enum::num_events) {
    return previous_score >= current_score;
  } else if (scoring == scoring_enum::likelihood) {
    return previous_score <= current_score;
  } else {
    return false;
  }
}
