#ifndef SCORING_H
#define SCORING_H

#include <string>

#include "history.h"

double get_score(History* history, std::string atoms_filename, std::string align_dir);
double get_initial_score();
double get_negative_difference(double previous_score, double current_score);
bool is_better_score(double previous_score, double current_score);

#endif
