#include <fstream>
#include <cassert>

#include "ghistory.h"
#include "files.h"

using namespace std;

GHistory::~GHistory() {
    clean();
}

void GHistory::clean() {
    if (sequences.size()) {
        this->types = sequences[0]->retype_atoms(0);
        for(auto t : types) delete t;
        types.clear();
    }
    for(auto s : sequences) delete s;
    sequences.clear();
    for(auto e : events) delete e;
    events.clear();
    GAtomType::reset_ids();
}

void GHistory::generate_random(double time, int sequence_length) {
    clean();
    sequences.push_back(new Sequence(0.0, sequence_length));
    int event_count = 0;
    GEventRoot root;
    events.push_back(new HEvent("e"+to_string(++event_count), &root, sequences.back()));

    double current_time = 0.0;
    while(current_time < time - epsilon) {
        Sequence* seq = sequences.back();
        GEvent* event = Model::instance()->get_random_event(seq->length());
        while (event->get_length() == 0) {
            delete event;
            event = Model::instance()->get_random_event(seq->length());
        }
        if (event->get_time(current_time) > time) {
            delete event;
            event = new GEventLeaf(time-current_time);
        }
        current_time = event->get_time(current_time);
        Sequence* nextseq = event->perform(seq);
        sequences.push_back(nextseq);
        events.push_back(new HEvent("e"+to_string(++event_count), event, nextseq));
        delete event;
        if (nextseq->length() == 0) break;
    }
    this->types = sequences.back()->retype_atoms(Model::instance()->length_threshold);
    assert(events.size() == sequences.size());
    For(i, events.size())
        events[i]->compute_atoms(i?events[i-1]:nullptr, i?sequences[i-1]:nullptr, sequences[i]);
    for(int i = events.size()-1; i>=0; --i)
        events[i]->compute_atom_ids(sequences[i]);

    if (debugging) {
        for(auto s : sequences) cout << *s;
        cout << endl;
        for(auto s : sequences) s->write_atoms_short();
        cout << endl;
        for(auto e : events) cout << *e;
        cout << endl;
        sequences.back()->write_atoms_short();
        cout << sequences.size() << " " << sequences.back()->atom_count() << " "
            << sequences.back()->length() << endl;
    }
}

void GHistory::save_to_files(string basename, string id) {
    if (id.size()) basename += "-" + id;
    ofstream f;
    cout << "Saving " << basename << endl;

    for(string species : {"unicorn"}) {
        f.open(basename+"-"+species+".dna", fstream::out);
        f << ">" << species << endl;
        write_final_sequence(f);
        f.close();
    }
    f.open(basename+".atoms", fstream::out);
    write_atoms(f);
    f.close();

    string basename_dir = basename + "/";
    remove_directory(basename_dir, "aln");
    create_directory(basename_dir);
    write_atoms_align(basename_dir);
    write_atoms_align_nexus(basename_dir);

    f.open(basename+".nhistory", fstream::out);
    write_events(f);
    f.close();
    f.open(basename+".stats", fstream::out);
    write_stats(f);
    f.close();

    cout << "       " << basename << " saved" << endl;
}

void GHistory::write_stats(ostream& os) {
    os << "Number of events: " << events.size() << endl;
    int atc = 0;
    set<int> att;
    ForGAtom(atom, sequences.back()) {
        if (atom->get_id() != 0) {
            atc++;
            att.insert(abs(atom->get_id()));
        }
    }
    os << "Number of atoms: " << atc << endl;
    os << "Number of types: " << att.size() << endl;
    int dashes = 0;
    ForGAtom(atom, sequences.back()) {
        if (atom->get_id() == 0) os << "(";
        os << atom->length();
        if (atom->get_id() == 0) os << ")";
        os << " ";
        dashes += atom->length() - atom->com_length();
    }
    os << endl;
    ForGAtom(atom, sequences.back()) {
        if (atom->get_id() == 0) os << "(";
        os << atom->com_length();
        if (atom->get_id() == 0) os << ")";
        os << " ";
    }
    os << endl;
    os << "length=" << sequences.back()->com_length() << " dels=" << dashes << endl;
}

void GHistory::write_final_sequence(ostream& os, const string& sep) {
    assert(sequences.size());
    sequences.back()->write_dna(os, sep);
}

void GHistory::write_atoms(ostream& os) {
    assert(sequences.size());
    sequences.back()->write_atoms(os);
}

void GHistory::write_atoms_align(string basepath) {
    assert(sequences.size());
    map<int, fstream> files;
    for(const auto& type : types)
        files[type->id].open(basepath+to_string(type->id)+".aln", fstream::out);
    ForGAtom(atom, sequences.back()) {
        int id = atom->get_type()->id;
        if (id) {
            if (!files.count(id)) {
                cerr << "Error: Not found type for atom\n";
                exit(1);
            }
            files[id] << ">" << atom->get_name() << endl
                      << atom->get_dna() << endl;
        }
    }
    for(auto& f : files) f.second.close();
}

// Nexus format is used as input to MrBayes and its syntax is described here:
// http://mrbayes.sourceforge.net/wiki/index.php/Tutorial#Getting_Data_into_MrBayes
void GHistory::write_atoms_align_nexus(string basepath) {
    assert(sequences.size());
    map<int, vector<GAtom*>> atoms_to_write;

    ForGAtom(atom, sequences.back()) {
        int typeId = atom->get_type()->id;
        if (typeId) {
            atoms_to_write[typeId].push_back(atom);
        }
    }

    for (const auto& [typeId, atoms] : atoms_to_write) {
        fstream file;
        file.open(basepath + to_string(typeId) + ".nex", fstream::out);

        file << "#NEXUS" << endl << "begin data;" << endl;
        file << "  dimensions ntax=" << atoms.size() << " nchar="
            << atoms[0]->length() << ";" << endl;
        file << "  format datatype=dna interleave=no gap=-;" << endl;
        file << "  matrix" << endl;

        for (auto atom : atoms) {
            file << "  " << atom->get_name() << " " << atom->get_dna() << endl;
        }

        file << "  ;" << endl << "end;" << endl;
        file.close();
    }
}

void GHistory::write_events(ostream& os) {
    for(auto e : events) os << *e;
}
