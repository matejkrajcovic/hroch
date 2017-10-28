#ifndef CONSTANTS_H
#define CONSTANTS_H

#include<cstdio>
#include<iostream>
#include<sstream>
#include<fstream>
#include<string>
#include<algorithm>
#include<vector>
#include<map>
#include<set>
#include<cmath>
#include<cassert>
#include<queue>
#include<tuple>
using namespace std;

class GAtom;
class GAtomType;
class GEvent;
class Sequence;
class Candidate;
class Dynamics;
class HAtom;
class HEvent;
class History;
class GHistory;
class CherryForest;
class Machine;
class ScoringData;

typedef pair<int, int> pii;

typedef vector<int> vi;
typedef vector<double> vdo;
typedef vector<vdo> vvdo;

#include"utils.h"
#include"random.h"
#include"gatom.h"
#include"gevent.h"
#include"hatom.h"
#include"hevent.h"
#include"cherry.h"
#include"sequence.h"
#include"history.h"
#include"ghistory.h"
#include"dynamics.h"
#include"score.h"
#include"machine.h"

extern string datapath;

extern int LOWER_RANGE;
extern int UPPER_RANGE;
extern string TEST_CASE;

extern int error_happened;
extern int debuging;
extern int stats;
extern int do_cheeryness;
#define SPECIAL_TRAINING 47
extern int strict_compare;
extern int num_bases;
extern char bases[];
extern int base_id[256];
extern char base_inv[256];
extern double epsilon;

void setup_constants();
set<string> parse_arguments(int argc, char **argv);

#endif
