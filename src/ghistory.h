// Representation of histories used for generated data

#ifndef GHISTORY_H
#define GHISTORY_H

#include <vector>
#include <set>
#include <string>
#include <ostream>

#include "sequence.h"
#include "hevent.h"
#include "random.h"

class GHistory {
    std::vector<Sequence*> sequences;
    std::set<GAtomType*> types;
    std::vector<HEvent*> events;

    void clean();
    void write_stats(std::ostream& os);
    void write_final_sequence(std::ostream& os, const std::string& sep = "\n");
    void write_atoms(std::ostream& os);
    void write_atoms_align(std::string basepath);
    void write_atoms_align_nexus(std::string basepath);
    void write_events(std::ostream& os);
public:
    ~GHistory();
    void generate_random(double time, int sequence_length);
    void save_to_files(std::string basename, std::string id = "");
};

#endif
