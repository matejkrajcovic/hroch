#ifndef ATOMTREE_H
#define ATOMTREE_H

#include<string>
#include"number.h"
#include "function.h"
using namespace std;

namespace likelihood {
class Node{
public:
    Node *left, *right, *parent;
    string data;                        // v liste je to sekvencia atomu, inde ""
    double edge_time;                   // dlzka hrany do otca
    double probability[4];
    Node();
    bool is_leaf();
    bool is_del();
    char get_residue(int position);

  //begin{Optimalizacia}
    int leaf_from, leaf_to;             // ktory interval listov zodpoveda mojmu podstromu

    bool trained[4];                    // uz poznam hodnotu pre jednotvarny podstrom so znakov v indexe
    double train_value[4][4];           // ake su naucena hodnoty

  //end{Optimalizacia}

    //uprava dlzok hran
    Event* event;
    Function function_probability[4];

    vdo allprobs[4];                    // ma podobnu funkciu ako probability, pouzivame vtedy,
                                        // ked si potrebujeme pamatat naraz cely riadok

    //vypis
    int print(int num, int next_num);
};

double change_probability(int a, int b, double t);

class AtomTree{
    int dg_docnt, dg_tracnt, dg_skipcnt;

  //begin{Optimalizacie}
    vector<Node*> leaf_order;           // listy v DFS poradi zlava doprava
    vector<int> last_same;              // posledne rovnake pismenko v stlpci
    void dfs_lf(Node* v);

    int last_ist_pos;                   // posledna hodnota pos, s ktorou bola zavolana is_single_type
    int get_single_type(int pos, int from, int to);

    void optimize(Node* vertex = NULL); // skomprimuje dlhe cesty
    void compute_leaf_order();          // porata leaf_order a opt. data v Node
    void clear_leaf_order();            // zmaze leaf_order
  //end{Optimalizacie}

    Number residue_likelihood(Node* vertex, int position);
public:
    Node* root;
    int atom_length;

    Number likelihood();
    ~AtomTree();
    void dfs_delete(Node* v);

    //zmena dlzok hran
    Function function_likelihood(Event* edit, double max_time);
    Function residue_function_likelihood(Event* edit, double max_time, Node* vertex, int position);
    //ternarna zmena dlzok
    Number ternary_likelihood(Event* edit, double max_time, double time);
    Number residue_ternary_likelihood(Event* edit, double max_time, double time, Node* vertex, int position);

    //vypis
    void print();
};
}
#endif
