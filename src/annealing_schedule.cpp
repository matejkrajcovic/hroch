#include <limits>

#include "annealing_schedule.h"
#include "constants.h"
#include "random.h"
#include "scoring.h"

using namespace std;

bool AnnealingSchedule::accept_new(const double score) {
  bool will_accept;
  if (score > previous_score) {
    will_accept = true;
  } else {
    will_accept = (temperature / (previous_score - score)) > random_double();
  }

  if (will_accept) {
    previous_score = score;
  }

  return will_accept;
}

SimpleAnnealingSchedule::SimpleAnnealingSchedule() {
  temperature = starting_temperature;
  total_steps = annealing_steps;
  current_step = 1;
  change = starting_temperature / total_steps;
  previous_score = get_initial_score();
}

bool SimpleAnnealingSchedule::finished() {
  return current_step == total_steps;
}

void SimpleAnnealingSchedule::update_temperature() {
  temperature -= change;
}

bool SimpleAnnealingSchedule::accept_new(const double score) {
  bool will_accept = AnnealingSchedule::accept_new(score);
  update_temperature();
  current_step++;
  return will_accept;
}

AdvancedAnnealingSchedule::AdvancedAnnealingSchedule() {
  temperature = starting_temperature;
  previous_score = get_initial_score();
}

bool AdvancedAnnealingSchedule::finished() {
  return temperature < 0.01;
}

void AdvancedAnnealingSchedule::update_temperature() {
  temperature *= 0.95;
}

bool AdvancedAnnealingSchedule::accept_new(const double score) {
  bool will_accept = AnnealingSchedule::accept_new(score / 10);
  update_temperature();
  return will_accept;
}

BaselineAdvancedAnnealingSchedule::BaselineAdvancedAnnealingSchedule() {
  temperature = starting_temperature;
  previous_score = get_initial_score();
}

bool BaselineAdvancedAnnealingSchedule::finished() {
  return temperature < 0.01;
}

void BaselineAdvancedAnnealingSchedule::update_temperature() {
  temperature *= 0.95;
}

bool BaselineAdvancedAnnealingSchedule::accept_new(const double score) {
  bool will_accept = score > previous_score;
  if (will_accept) {
    previous_score = score;
  }
  update_temperature();
  return will_accept;
}

NewAdvancedAnnealingSchedule::NewAdvancedAnnealingSchedule(double initial_score) {
  temperature = 0.4;
  current_step = 0;
  total_steps = annealing_steps;
  previous_score = initial_score;
}

bool NewAdvancedAnnealingSchedule::finished() {
  return current_step == total_steps;
}

void NewAdvancedAnnealingSchedule::update_temperature() {
  temperature /= 3;
}

bool NewAdvancedAnnealingSchedule::accept_new(const double score) {
  bool will_accept;
  if (is_better_score(previous_score, score)) {
    will_accept = true;
  } else {
    will_accept = random_double() < exp(get_negative_difference(previous_score, score) / temperature);
    if (will_accept) {
      cout << "accepting worse" << endl;
    }
  }
  if (will_accept) {
    previous_score = score;
  }
  current_step++;
  return will_accept;
}

int NewAdvancedAnnealingSchedule::get_progress() {
  int old_progress = floor((double) (current_step - 1) / ceil((double) total_steps / 3));
  int new_progress = floor((double) current_step / ceil((double) total_steps / 3));
  if (old_progress != new_progress) {
    update_temperature();
  }
  return new_progress;
}
