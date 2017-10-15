// Representation of events used when generating data

#ifndef GEVENT_H
#define GEVENT_H

#include"constants.h"

class GEvent {
protected:
    double time_interval;
    int from, to;
public:    
    virtual Sequence* perform(Sequence* sequence) = 0;
    virtual string name() = 0;
    virtual string get_string() { return this->name(); }
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
    virtual string name() { return "dup"; }
    virtual Sequence* perform(Sequence* sequence);
    virtual string get_string() {
        return this->name() + " " + to_string(from) + " " + to_string(to) +
            " " + to_string(cpos);
    }
    GEventDup(int from, int to, int cpos, double time_interval);
};

class GEventDupi: public GEventDup {
protected:
public:
    virtual string name() { return "dupi"; }
    virtual Sequence* perform(Sequence* sequence);
    GEventDupi(int from, int to, int cpos, double time_interval);
};

class GEventDel: public GEvent {
public:
    virtual string name() { return "del"; }
    virtual Sequence* perform(Sequence* sequence);
    virtual string get_string() { return this->name() + to_string(from) + " " + to_string(to); }
    GEventDel(int from, int to, double time_interval);
};

class GEventLeaf: public GEvent {
public: 
    virtual string name() { return "leaf"; }
    virtual Sequence* perform(Sequence* sequence);
    GEventLeaf(double time_interval = 0.0);
};

class GEventRoot: public GEvent {
public:
    virtual string name() { return "root"; }
    virtual Sequence* perform(Sequence* sequence);
    GEventRoot(); 
};

ostream& operator<<(ostream& os, GEvent& event);

#endif
