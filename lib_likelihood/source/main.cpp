#include"constants.h"
#include"tests.h"
#include"history.h"
#include<fstream>
using namespace std;

char bases[] = {'A', 'C', 'G', 'T'};
int debuging = 0;

string f_params, f_atoms, f_hist, p_align, f_output, f_result, f_orig;

string next_line(istream& is){
    string line;
    do {
        getline(is, line);
    } while (line[0] == '#');
    return line;
}

void cleanup(map<vi, double>& M, double threshold, double multiply = 1.0){
    vector<pair<vi, double>> V;
    for(auto x : M) {
        x.second*=multiply;
        if (x.second > threshold) V.push_back(x);
    }
    M.clear();
    M.insert(V.begin(), V.end());
}

int main(int argc, char **argv){
    random_init();
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " config-number test-number iterations " <<  endl;
        return 1;
    }
    string configname = argv[1];
    string testnumber = argv[2];
    int iterations = stoi(string(argv[3]));

    History* h = new History();
    ifstream fc;
    fc.open("configs/"+configname + ".config");
    mustbe(fc.is_open(), "can't open config file: configs/"+configname+".config");
    f_params = next_line(fc);
    f_atoms = next_line(fc)+testnumber+".atoms";
    p_align = next_line(fc)+testnumber;
    f_hist = next_line(fc)+testnumber+".nhistory";
    f_output = next_line(fc)+configname+"-"+testnumber;
    f_result = f_output+".history";
    f_orig = f_output+".orig";
    f_output += ".samples";
    h->scoring_length = stoi(next_line(fc));
    h->scoring_similar= stoi(next_line(fc));
    h->scoring_cherryness = stoi(next_line(fc));
    h->scoring_boundary = stoi(next_line(fc));
    h->scoring_oldhistory = stoi(next_line(fc));
    h->scoring_divisor = stoi(next_line(fc));
    h->felsenstein_on = stoi(next_line(fc));
    
    ifstream in;
    in.open(f_params, ifstream::in);
    mustbe(in.is_open(), "can't open parameter file: "+f_params);
    h->read_parameters(in);
    in.close();
    in.open(f_atoms, ifstream::in);
    mustbe(in.is_open(), "can't open atoms file: "+f_atoms);
    h->read_atom_lists(in);
    in.close();
        
    for(auto atom : h->existing_atoms){
        string f_align = p_align+"/"+to_string(atom)+".aln";
        in.open(f_align, ifstream::in);
        mustbe(in.is_open(), "can't open alignment file: " + f_align);
        h->read_atom_alignments(in, atom);
        in.close();
    }

    in.open(f_hist, ifstream::in);
    if(in.is_open()){
        h->read_events(in);
        in.close();
        h->build_trees();
        h->write_events(cout);
        Number lik = h->likelihood_all();
        int old_felsenstein_on = h->felsenstein_on;
        h->felsenstein_on = 1;
        Number fulllik = h->likelihood_all();
        h->felsenstein_on = old_felsenstein_on;

        cout << "orig lik: " << lik.data << endl;
        cout << "full: " << fulllik.data << endl;

        ofstream orig (f_orig, ofstream::out);
        orig << "Original likelihood: " << lik.data << " " << fulllik.data << endl;
        h->write_events(orig);
        orig.close();
        //zaba
        For(i,100) {
            h->edit_edge_lengths();
            Number lik2 = h->likelihood_all();
            cout << i << " lik po zmene: " << lik2.data << endl;
        }
        h->write_events(cout);
        Number lik2 = h->likelihood_all();
        cout << "povodna likelihood bola: " << lik.data << endl;
        cout << "lik po vsetkych zmenach: " << lik2.data << endl;
        //end zaba
        h->clear_events();
    }else{
        cout << "Original history is unknown" << endl;
    }
    ofstream output (f_output, ofstream::out);
    mustbe(output.is_open(), "can't open output file: "+f_output);
   
    map<string, Event*> bestev;
    Number bestlike = logNumber(-1e15), lastlike = logNumber(-1e15);
    double sum = 0; int cnt = 0;
    For(ite, iterations){
        h->compute_events();
        cerr << "computed" << endl;
        h->write_events(cerr);
        h->build_trees();
        cerr << "builded" << endl;
        Number lik = h->likelihood_all();
        lik.print(stderr);
        output << "New sample with likelihood: " << lik.data << endl;
        h->write_events(output);
        output << endl;        

        cout << "lik: " << lik.data;
        sum+=lik.data; cnt++;
        if (lastlike < lik || random_double()>=0.5+0.5/iterations*ite){
            cout << "    !! novy !!";
            for(auto x : h->old_events) {
                for(auto y : x.first) fprintf(stderr,"%d ",y);
                fprintf(stderr,"with value %lf\n", x.second);
            }
            if (h->scoring_oldhistory){
                cleanup(h->old_events, 0.1, 0.6);
                cleanup(h->old_seqences, 0.1, 0.6);
                for(auto ev : h->events) h->add_to_old(ev.second,1+
                        max(lik.data,bestlike.data)/lik.data);
                for(auto ev : bestev) h->add_to_old(ev.second,0.1);
            }
        }
        cleanup(h->old_events, 0.1, 0.9);
        cleanup(h->old_seqences, 0.1, 0.9);
        lastlike = lik;
        cout << "            " << SIZE(h->old_events) << " " << SIZE(h->old_seqences) << endl;


        if (bestlike<lik) {
            bestlike = lik;
            for(auto ev : bestev) delete ev.second;
            bestev.clear();
            bestev = h->events;
            h->events.clear(); // toto je nutne tu mat
            h->clear_events();
        }else{
            h->clear_events();
        }
    }
    output.close();
    
    cout << "Best found history: " << endl;
    bestlike.print();
    h->events = bestev;
    h->write_events(cout);
    h->trust_me();
    h->build_trees();
    Number lik = h->likelihood_all();
    cout << "avge lik: " << sum/cnt << endl;
    cout << "best lik: "; lik.print();
    int old_felsenstein_on = h->felsenstein_on;
    h->felsenstein_on = 1;
    Number fulllik = h->likelihood_all();
    cout << "total: "; fulllik.print();
    h->felsenstein_on = old_felsenstein_on;

    ofstream result (f_result);
    mustbe(result.is_open(), "can't open output file: "+f_result);
    result << "Best likelihood: " << lik.data << " " << fulllik.data << endl;
    h->write_events(result);
    result.close();
}
