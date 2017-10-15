#include"ghistory.h"
#include"files.h"

void GHistory::clean() {
    if (SIZE(sequences)) {
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
    assert(SIZE(events) == SIZE(sequences));
    For(i, SIZE(events)) 
        events[i]->compute_atoms(i?events[i-1]:nullptr, i?sequences[i-1]:nullptr, sequences[i]);
    for(int i = SIZE(events)-1; i>=0; --i)
        events[i]->compute_atom_ids(sequences[i]);

    if (debuging) {
        for(auto s : sequences) cout << *s;
        cout << endl;
        for(auto s : sequences) s->write_atoms_short();
        cout << endl;
        for(auto e : events) cout << *e;
        cout << endl;
        sequences.back()->write_atoms_short();
        cout << SIZE(sequences) << " " << sequences.back()->atom_count() << " " 
            << sequences.back()->length() << endl;
    }
}

void GHistory::save_to_files(string basename, string id) {
    if (SIZE(id)) basename += "-" + id;
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
    write_atoms_align(basename+"/"); 
    f.open(basename+".nhistory", fstream::out);
    write_events(f);
    f.close();
    f.open(basename+".stats", fstream::out);
    write_stats(f);
    f.close();

    cout << "       " << basename << " saved" << endl; 
}

void GHistory::write_stats(ostream& os) {
    os << "Number of events: " << SIZE(events) << endl;
    int atc = 0;
    set<int> att;
    ForGAtom(atom, sequences.back()) {
        if (atom->get_id() != 0) {
            atc++;
            att.insert(abs(atom->get_id()));
        }
    }
    os << "Number of atoms: " << atc << endl;
    os << "Number of types: " << SIZE(att) << endl;
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
    assert(SIZE(sequences));
    sequences.back()->write_dna(os, sep);
}

void GHistory::write_atoms(ostream& os) {
    assert(SIZE(sequences));
    sequences.back()->write_atoms(os); 
}

void GHistory::write_atoms_align(string basepath) {
    assert(SIZE(sequences));
    map<int, fstream> files;
    remove_directory(basepath, "aln");
    create_directory(basepath);
    for(const auto& type : types) 
        files[type->id].open(basepath+to_string(type->id)+".aln", fstream::out);
    ForGAtom(atom, sequences.back()) {
        int id = atom->get_type()->id;
        if (id) {
            mustbe(files.count(id), "Not found type for atom");
            files[id] << ">" << atom->get_name() << endl
                      << atom->get_dna() << endl;
        }
    }
    for(auto& f : files) f.second.close();
}

void GHistory::write_events(ostream& os) {
    for(auto e : events) os << *e;
}
