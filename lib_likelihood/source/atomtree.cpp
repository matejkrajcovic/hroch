#include"atomtree.h"

Node::Node(){
    left = right = parent = NULL;
    edge_time = 0.0;
    data = "";
    For(i, 4) trained[i] = false;
    event = NULL;
}

bool Node::is_leaf() {
    if(left == NULL && right == NULL) return true;
    return false;
}
bool Node::is_del() {
    return data.size()==0 && is_leaf();
}

char Node::get_residue(int position) {
    //return data.size()?data[position]:'N';
    return data[position];
}

void AtomTree::optimize(Node* vertex){
    if (vertex == NULL) vertex = root;
    if (vertex->is_leaf()) 
        return;
    if (vertex->left == NULL) swap(vertex->left, vertex->right);
    if (vertex->left == NULL) {
        cerr << "I don't have sons!" << endl;
        return;
    }
    optimize(vertex->left);

    if (vertex->right != NULL){
        optimize(vertex->right); 
    }else if(vertex!=root) {
        vertex->left->parent = vertex->parent;
        if (vertex->parent->left == vertex)
            vertex->parent->left = vertex->left;
        if (vertex->parent->right == vertex)
            vertex->parent->right = vertex->left;
        vertex->left->edge_time += vertex->edge_time;
        delete vertex;
    }
}


// volanie s rovnakym pos ako minule trva O(1)
// volanie s inym pos trva O(SIZE(leaf_order))
// vrati -1 alebo 0..3 toho, co je rovnake
int AtomTree::get_single_type(int pos, int from, int to){
    if (last_ist_pos != pos){
        last_same = vector<int>(SIZE(leaf_order));
        if (SIZE(last_same)){
            last_same[0] = 0;
            For(i, SIZE(last_same)-1){
                last_same[i+1] = (leaf_order[i]->get_residue(pos)==leaf_order[i+1]->get_residue(pos))?
                        last_same[i]:i+1;
            }
        }
        last_ist_pos = pos;
    }
    if (from<to && from>=last_same[to-1]){
        For(i, 4) if (leaf_order[from]->get_residue(pos)==bases[i])
            return i;
    }
    return -1;
}

void AtomTree::dfs_lf(Node* v){
    v->leaf_from = SIZE(leaf_order);
    if (v->is_leaf() && SIZE(v->data)){
        leaf_order.push_back(v);
    }

    For(i, 4) v->trained[i] = false;

    if(v->left != NULL) dfs_lf(v->left);
    if(v->right != NULL) dfs_lf(v->right);

    v->leaf_to = SIZE(leaf_order);
}

void AtomTree::compute_leaf_order(){
    leaf_order.clear();
    last_ist_pos = -1;
    dfs_lf(root);
}

void AtomTree::clear_leaf_order(){
    leaf_order.clear();
}

double change_probability(int a, int b, double t) {
    if(a==b) return (1.0 + 3.0 * exp(-4.0 * FEL_ALPHA * t))/4.0;
    return (1.0 - exp(-4.0 * FEL_ALPHA * t))/4.0;
}

Number AtomTree::residue_likelihood(Node *vertex, int position) {
    For(i,4) vertex->probability[i] = 0.0;
    if(vertex->is_leaf()) {
        char current_base = vertex->is_del()?'N':vertex->get_residue(position);
        For(i,4) {
            if(current_base == 'N' || current_base == '-') 
                vertex->probability[i] = 1.0;
            else if(current_base == bases[i]) vertex->probability[i] = 1.0;
            else vertex->probability[i] = 0.0;
        }
    } else if(vertex->left == NULL) {
        fprintf(stderr, "Something is wrong, I should have left son\n");
    } else if(vertex->right == NULL) {
        assert(vertex->left->edge_time > 1e-9);
        residue_likelihood(vertex->left, position);
        For(a, 4) For(b, 4)
            vertex->probability[a] += 
                vertex->left->probability[b] * 
                change_probability(b, a, vertex->left->edge_time);
        
    } else {
        assert(vertex->left->edge_time > 1e-9);
        assert(vertex->right->edge_time > 1e-9);
        int training_on = SIZE(leaf_order)?get_single_type(position, 
                vertex->leaf_from, vertex->leaf_to):-1;

        if (training_on<0 || (!vertex->trained[training_on])){
            residue_likelihood(vertex->left, position);
            residue_likelihood(vertex->right, position);
            For(a,4) For(b,4) For(c,4)
                vertex->probability[a] += 
                    vertex->left->probability[b] *
                    vertex->right->probability[c] *
                    change_probability(b, a, vertex->left->edge_time) *
                    change_probability(c, a, vertex->right->edge_time);
            dg_docnt++;
        }
        if (training_on>=0){
            if (vertex->trained[training_on]){
                For(a,4) vertex->probability[a] = vertex->train_value[training_on][a];
                dg_skipcnt++;
            }else{
                vertex->trained[training_on] = true;
                For(a,4) vertex->train_value[training_on][a] = vertex->probability[a];
                dg_tracnt++;
            }
        }
    }
    double result_likelihood = 0.0;
    For(i,4) result_likelihood += vertex->probability[i];// / 4.0;
    
    /*cerr << int(vertex)%1000 << " " << int(vertex->parent)%1000 
         << " " << position << " " << result_likelihood << ":";
    For(i, 4) cerr << " " << vertex->probability[i];
    cerr << "    " << int(vertex->left)%1000 << " " << int(vertex->right)%1000 
         << " '" << vertex->data << endl;
    */
    
    return Number(result_likelihood);
}

Number AtomTree::likelihood(){
    optimize();
    compute_leaf_order();
    dg_docnt = dg_tracnt = dg_skipcnt = 0;

    Number result_likelihood = Number(1.0);
    For(i,atom_length) {
        Number x;
        result_likelihood *= (x = residue_likelihood(root, i));
        //x.print();
    }

    //printf("%6d%6d%6d x%6d %.3lf   ", dg_docnt, dg_tracnt, dg_skipcnt,
    //        atom_length, dg_docnt/double(atom_length*SIZE(leaf_order)));
    clear_leaf_order();
    //result_likelihood.print();
    return result_likelihood;
}

AtomTree::~AtomTree(){
    dfs_delete(root);
}

void AtomTree::dfs_delete(Node* node){
    if (node==NULL) return;
    dfs_delete(node->left);
    dfs_delete(node->right);
}

//vypis

int Node::print(int num, int next_num) {
    int sons=0;
    if(left!=NULL) sons++;
    if(right!=NULL) sons++;
    fprintf(stderr,"%d:",num);
    For(i,sons) fprintf(stderr," %d",next_num+i); fprintf(stderr,"\n");
    int p=next_num+sons;
    if(left!=NULL) p=left->print(next_num,p);
    if(right!=NULL) p=right->print(next_num+1,p);
    return p;
}

void AtomTree::print() {
    root->print(0,1);
}
