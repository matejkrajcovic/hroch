#ifndef ANNEALING_SCHEDULE_H
#define ANNEALING_SCHEDULE_H

class AnnealingSchedule {
protected:
  double temperature;
  double previous_score;
  virtual void update_temperature() = 0;
public:
  virtual ~AnnealingSchedule() {};
  virtual bool accept_new(const double score) = 0;
  virtual bool finished() = 0;
};

class SimpleAnnealingSchedule : public AnnealingSchedule {
  int current_step;
  int total_steps;
  double change;
  void update_temperature();
public:
  SimpleAnnealingSchedule();
  bool finished();
  bool accept_new(const double score);
};

class AdvancedAnnealingSchedule : public AnnealingSchedule {
  void update_temperature();
public:
  AdvancedAnnealingSchedule();
  bool accept_new(const double score);
  bool finished();
};

class BaselineAdvancedAnnealingSchedule : public AnnealingSchedule {
  void update_temperature();
public:
  BaselineAdvancedAnnealingSchedule();
  bool accept_new(const double score);
  bool finished();
};

#endif
