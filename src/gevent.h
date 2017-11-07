// Representation of events used when generating data

#ifndef GEVENT_H
#define GEVENT_H

#include <string>
#include <ostream>

#include "sequence.h"

class GEvent {
protected:
    double time_interval;
    int from, to;
public:
    virtual Sequence* perform(Sequence* sequence) = 0;
    virtual std::string name() = 0;
    virtual std::string get_string() { return this->name(); }
    int get_length() { return to-from; }
    double get_time(double time_start = 0.0);
    GEvent(double time_interval);
    virtual ~GEvent() {};
};

class GEventDup: public GEvent {
protected:
    int cpos;
    Sequence* iperform(Sequence* sequence, bool invert);
public:
    virtual std::string name() { return "dup"; }
    virtual Sequence* perform(Sequence* sequence);
    virtual std::string get_string() {
        return this->name() + " " + std::to_string(from) + " " + std::to_string(to) +
            " " + std::to_string(cpos);
    }
    GEventDup(int from, int to, int cpos, double time_interval);
};

class GEventDupi: public GEventDup {
protected:
public:
    virtual std::string name() { return "dupi"; }
    virtual Sequence* perform(Sequence* sequence);
    GEventDupi(int from, int to, int cpos, double time_interval);
};

class GEventDel: public GEvent {
public:
    virtual std::string name() { return "del"; }
    virtual Sequence* perform(Sequence* sequence);
    virtual std::string get_string() { return this->name() + std::to_string(from) + " " + std::to_string(to); }
    GEventDel(int from, int to, double time_interval);
};

class GEventLeaf: public GEvent {
public:
    virtual std::string name() { return "leaf"; }
    virtual Sequence* perform(Sequence* sequence);
    GEventLeaf(double time_interval = 0.0);
};

class GEventRoot: public GEvent {
public:
    virtual std::string name() { return "root"; }
    virtual Sequence* perform(Sequence* sequence);
    GEventRoot();
};

std::ostream& operator<<(std::ostream& os, GEvent& event);

#endif
