#include"constants.h"
#include"history.h"
#include<fstream>
using namespace std;

namespace likelihood {
double calculate_reconstruction_likelihood(string atoms_filename, string align_dir, stringstream& reconstruction_stream) {
    random_init();

    HistoryLikelihood* h = new HistoryLikelihood();
    string f_atoms = atoms_filename;
    string p_align = align_dir;
    h->scoring_length = 1;
    h->scoring_similar= 1;
    h->scoring_cherryness = 0;
    h->scoring_boundary = 1;
    h->scoring_oldhistory = 1;
    h->scoring_divisor = 1;
    h->felsenstein_on = 1;

    // make sure this values are the same as in random.h
    h->set_parameters(200, 306718, 14307, 0.3902, 0.05);

    ifstream in;
    in.open(f_atoms, ifstream::in);
    mustbe(in.is_open(), "can't open atoms file: "+f_atoms);
    h->read_atom_lists(in);
    in.close();

    for(auto atom : h->existing_atoms){
        string f_align = p_align+"/"+to_string(atom)+".aln.fa";
        in.open(f_align, ifstream::in);
        mustbe(in.is_open(), "can't open alignment file: " + f_align);
        h->read_atom_alignments(in, atom);
        in.close();
    }

    h->read_events(reconstruction_stream);

    h->build_trees();
    //h->write_events(cout);
    //Number lik = h->likelihood_all();
    int old_felsenstein_on = h->felsenstein_on;
    h->felsenstein_on = 1;
    //Number fulllik = h->likelihood_all();
    h->felsenstein_on = old_felsenstein_on;

    //cout << "orig lik: " << lik.data << endl;
    //cout << "full: " << fulllik.data << endl;

    For(i,100) {
        h->edit_edge_lengths();
        //Number lik2 = h->likelihood_all();
        //cout << i << " lik po zmene: " << lik2.data << endl;
    }
    Number lik2 = h->likelihood_all();
    //cout << "povodna likelihood bola: " << lik.data << endl;
    //cout << "lik po vsetkych zmenach: " << lik2.data << endl;
    h->clear_events();

    return lik2.data;
}
}
