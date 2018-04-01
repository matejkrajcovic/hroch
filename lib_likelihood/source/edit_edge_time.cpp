#include "constants.h"
#include "function.h"
#include "event.h"
#include "atomtree.h"

double get_random_double(double upper) {
    double x = (double)rand()/(double)RAND_MAX;
    return upper*x;
}

void History::edit_edge_lengths() {
    this->build_trees();
    int number_of_events = events.size() - 2;
    int choosen_one = rand()%number_of_events;
    Event* edit = NULL;
    for(auto e : events) {
        if(e.second->type == "root" || e.second->type == "leaf") continue;
        if(choosen_one == 0) edit = e.second;
        choosen_one--;
    }
    double max_time = 0;
    for(auto e : events) {
        if(e.second == edit || e.second->parent == edit) max_time += e.second->edge_time;
    }
    cerr << "picked event: " << edit->name << " " << max_time << endl;
    //upravujem tuto udalost
    find_better_time(edit, max_time);
}

double newton_rapson(double x, Function* f, Function* fd, Function* fdd, int h, int nh, double lambda) {
    For(i,50) {
        Number derf = fd->evaluate(x), dderf = fdd->evaluate(x);
        Number shift2 = derf/dderf;
        double shift = shift2.true_value();
        x = x - shift;
    }
    return x;
}

double find_max(Function* f, Function* fd, Function* fdd, int h, int nh, double max_time, double lambda) {
    vector<double> O;
    For(i,100) {
        double x_0 = get_random_double(max_time);
        double best = newton_rapson(x_0, f, fd, fdd, h, nh, lambda);
        if(best<0) best = 0;
        if(best > max_time) best = max_time;
        O.push_back(best);
    }
    int p=0;
    For(i,O.size()) {
        if(f->evaluate(O[p]) < f->evaluate(O[i])) p=i;
    }
    return O[p];
}

void History::find_better_time(Event* edit, double max_time) {
    Number old_likelihood = likelihood_all();
    //poratam funkcnu likelihood pre kazdy genovy strom
    cerr << "ratam funkcnu likelihood" << endl;
    Function f = Function(1.0);
    f.print();
    for(auto tree : group_trees) {
        Function pom = tree.second->function_likelihood(edit, max_time);
        f = f * pom;
    }
    Function fd = f.derivate();
    Function fdd = fd.derivate();
    //zisti ako vyzera likelihood eventov
    cerr << "likelihood eventov" << endl;
    int happened=1,not_happend=0;
    for(auto e : events)
        if(e.second->parent == edit && e.second->parent->type == "leaf") not_happend++;
        else if(e.second->parent == edit) happened++;
    //najdi vhodny cas - maximum funkcie f
    cerr << "Newton-Rapson" << endl;
    double better_time = find_max(&f, &fd, &fdd, happened, not_happend, max_time, this->parameters["likelihood_event_rate"]);
    cerr << "better time: " << better_time << endl;
    //nastav cas eventu na toto
    double new_event_time;
    for(auto e : events) {
        if(e.second->parent == edit) new_event_time = e.second->event_time;
    }
    double cas_zaloha = edit->event_time;
    fprintf(stderr,"stary cas: %lf\n",edit->event_time);
    new_event_time -= better_time;
    fprintf(stderr,"novy cas: %lf\n",new_event_time);
    edit->event_time = new_event_time;
    build_trees();
    Number new_likelihood = likelihood_all();
    if(new_likelihood < old_likelihood) {
        edit->event_time = cas_zaloha;
        build_trees();
    }
}

//funkciova pravdepodobnost
Function function_change_probability(int a, int b, double time, Node* vertex, Event* edit, double max_time) {
    if(vertex->event == edit) {
        if(a==b) return Function(1.0/4.0,3.0*exp(-4.0*FEL_ALPHA*max_time)/4.0,0);
        return Function(1.0/4.0,-1.0*exp(-4.0*FEL_ALPHA*max_time)/4.0,0);
    }
    else if(vertex->parent != NULL && vertex->parent->event == edit) {
        if(a==b) return Function(3.0/4.0,1.0/4.0,-1);
        return Function(-1.0/4.0,1.0/4.0,-1);
    }
    else {
        if(a==b) return Function((1+3* exp(-4.0*FEL_ALPHA*time))/4.0);
        return Function((1-exp(-4.0*FEL_ALPHA*time))/4.0);
    }
}

//ratam funkcioveho Felsenstaina
Function AtomTree::residue_function_likelihood(Event* edit, double max_time, Node* vertex, int position) {
    For(i,4) vertex->function_probability[i] = Function(0.0,true);
    if(vertex->is_leaf()) {
        char current_base = vertex->is_del()?'N':vertex->get_residue(position);
        For(i,4) {
            if(current_base == 'N' || current_base == bases[i]) vertex->function_probability[i] = Function(1.0);
            else vertex->function_probability[i] = Function(0.0,true);
        }
    }
    else if(vertex->left == NULL)
        fprintf(stderr, "I should have left son\n");
    else if(vertex->right == NULL) {
        residue_function_likelihood(edit, max_time, vertex->left, position);
        For(a, 4) For(b, 4) {
            vertex->function_probability[a] = vertex->function_probability[a] +
                vertex->left->function_probability[b] * 
                function_change_probability(b, a, vertex->left->edge_time, vertex, edit, max_time);
        }
    }
    else {
        residue_function_likelihood(edit, max_time, vertex->left, position);
        residue_function_likelihood(edit, max_time, vertex->right, position);
        For(a, 4) For(b, 4) For(c, 4)
            vertex->function_probability[a] = vertex->function_probability[a] +
            vertex->left->function_probability[b] *
            vertex->right->function_probability[c] *
            function_change_probability(b, a, vertex->left->edge_time, vertex, edit, max_time) *
            function_change_probability(c, a, vertex->right->edge_time, vertex, edit, max_time);
    }

    Function res = Function(0.0,true);
    For(i,4) res = res + vertex->function_probability[i];
    return res;
}

Function AtomTree::function_likelihood(Event* edit, double max_time) {
    Function res = Function(1.0);
    For(i, atom_length) {
        res = res * residue_function_likelihood(edit, max_time, root, i);
    }
    return res;
}
