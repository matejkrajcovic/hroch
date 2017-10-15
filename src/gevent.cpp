#include"gevent.h"

double GEvent::get_time(double time_start) {
    return time_start + time_interval;
}

GEvent::GEvent(double time_interval) {
    this->time_interval = time_interval;
}

GEventDup::GEventDup(int from, int to, int cpos, double time_interval) : 
    GEvent(time_interval) {
    this->from = from;
    this->to = to;
    this->cpos = cpos;
}

GEventDupi::GEventDupi(int from, int to, int cpos, double time_interval) : 
    GEventDup(from, to, cpos, time_interval) {
}

GEventDel::GEventDel(int from, int to, double time_interval) : 
    GEvent(time_interval) {
    this->from = from;
    this->to = to;
}
GEventLeaf::GEventLeaf(double time_interval) :
    GEvent(time_interval) {
}
GEventRoot::GEventRoot() :
    GEvent(0) {
}

Sequence* GEventLeaf::perform(Sequence* sequence) {
    Sequence* s = new Sequence(sequence);
    s->mutate(time_interval);
    return s;    
}

Sequence* GEventRoot::perform(Sequence* sequence) {
    return sequence;
}

Sequence* GEventDup::iperform(Sequence* sequence, bool invert) {
    sequence->split_breakpoints({from, to, cpos});
    Sequence* s = new Sequence(sequence);
    s->mutate(time_interval);
    int pos = 0;
    GAtom *cpos_atom = nullptr, *cfront_atom = nullptr, *cback_atom = nullptr;
    GAtom *n_atom = s->first, *p_atom = nullptr;
    while(n_atom != nullptr) {
        if (pos == cpos) cpos_atom = p_atom;
        if (pos == from) {
            cfront_atom = new GAtom(n_atom->parent, n_atom->parent->get_dna());
            cback_atom = cfront_atom;
        }
        if (pos > from && pos < to) {
            cback_atom->next = new GAtom(n_atom->parent, n_atom->parent->get_dna());
            cback_atom = cback_atom->next;
        }
        pos += n_atom->length();
        p_atom = n_atom;
        n_atom = n_atom->next;
    }
    if (invert) {
        GAtom *a = nullptr, *b = cfront_atom;
        while(b != cback_atom) {
            b->invert();
            GAtom *c = b;
            b = b->next;
            c->next = a;
            a = c; 
        }
        b->invert();
        b->next = a;
        swap(cfront_atom, cback_atom);
    }
    if (cpos_atom == nullptr) {
       cback_atom->next = s->first; 
       s->first = cfront_atom;
    } else {
        cback_atom->next = cpos_atom->next;
        cpos_atom->next = cfront_atom;
    }
    return s;
}

Sequence* GEventDup::perform(Sequence* sequence) {
    return iperform(sequence, false);
}
Sequence* GEventDupi::perform(Sequence* sequence) {
    return iperform(sequence, true);
}

Sequence* GEventDel::perform(Sequence* sequence) {
    sequence->split_breakpoints({from, to});
    Sequence* s = new Sequence(sequence);
    s->mutate(time_interval);
    int pos = 0;
    GAtom* n_atom = s->first;
    GAtom* p_atom = nullptr;
    while(n_atom != nullptr) {
        if (pos >= from && pos < to) {
            if (p_atom) p_atom->next = n_atom->next;
            else s->first = n_atom->next;
            pos += n_atom->length();
            n_atom->unlink_parent();
            delete n_atom;
            n_atom = (p_atom)?p_atom->next:s->first;
        } else {
            pos += n_atom->length();
            p_atom = n_atom;
            n_atom = n_atom->next;
        }
    }
    return s;
}

ostream& operator<<(ostream& os, GEvent& event) {
    return os << event.get_string();
}
