#include<random>
#include"random.h"
#include<unistd.h>

default_random_engine generator;
uniform_real_distribution<double> ur_distribution;

void random_init(){
    //srand(time(NULL)+getpid());
    // 466414 - guaranted to be random, generated by 6 rolls of fair dice :)
    srand(466414);
    generator = default_random_engine(rand());
    ur_distribution = uniform_real_distribution<double>();
}

double random_double() {
    return ur_distribution(generator);
}
double random_double(double from, double to) {
    return double(generator()-generator.min())/double(generator.max()-generator.min()+1)*(to-from)+from;
}
int random_int(int from, int to) {
    return from + rand() % (to-from);
}

Model* Model::_instance = nullptr;

Model::Model() {
    length_distribution = geometric_distribution<int>(1./mean_len);
    distance_distribution = geometric_distribution<int>(1./mean_dist);
    time_distribution = exponential_distribution<double>(event_rate);
}
Model* Model::instance() {
    if (Model::_instance == nullptr) Model::_instance = new Model();
    return Model::_instance;
}
bool Model::is_random_del() {
    return random_double() < prob_del;
}
bool Model::is_random_inv() {
    return random_double() < prob_inv;
}
int Model::get_random_len() {
    return length_distribution(generator);
}
int Model::get_random_dist() {
    return distance_distribution(generator);
}
double Model::get_random_time() {
    return time_distribution(generator);
}
char Model::get_mutated_base(char base, double time) {
    return (random_double() < exp(-4.*mut_alpha*time))?base:bases[rand()%BASES];
}
bool Model::get_indel_happened(double time) {
    return (random_double() > exp(-0.02*time));
}

GEvent* Model::get_random_event(int sequence_length) {
    double time_interval = this->get_random_time();
    int length = this->get_random_len();
    while(length > sequence_length) length = this->get_random_len();
    int from = random_int(0, sequence_length - length + 1);

    if (this->is_random_del()) {
        return new GEventDel(from, from + length, time_interval);
    } else {
        int to = -1;
        while (to < 0 || to > sequence_length) {
            int dist = this->get_random_dist();
            to = (rand()%2)?from-dist:from+length+dist;
        }
        if (this->is_random_inv()) {
            return new GEventDupi(from, from+length, to, time_interval);
        } else {
            return new GEventDup(from, from+length, to, time_interval);
        }
    }
}
