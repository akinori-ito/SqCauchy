#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include "SQC_train.h"
#include "Adam.h"

using namespace std;

class SQCUpdate: public SQCauchy {
public:
  int count;
  SQCauchy *sqc;
  Adam *adam_weight;
  Adam **adam_mean;
  Adam **adam_var;
  SQCUpdate(SQCauchy *s) : SQCauchy(s->dim,s->nmix),sqc(s) {
    adam_weight = new Adam(nmix);
    adam_mean = new Adam*[nmix];
    adam_var = new Adam*[nmix];
    for (int i = 0; i < nmix; i++) {
      adam_mean[i] = new Adam(dim);
      adam_var[i] = new Adam(dim);
    }
  }
  void clear() {
    count = 0;
    for (int i = 0; i < nmix; i++) {
      weight[i] = 0;
      for (int k = 0; k < dim; k++) {
        mean[i][k] = 0;
        var[i][k] = 0;
      }
    }
  }
  void feed1(float *x) {
    double px = sqc->prob(x);
    for (int i = 0; i < nmix; i++) {
      double pxi = sqc->prob_i(x,i);
      weight[i] += pxi/px;
      for (int k = 0; k < dim; k++) {
        double m = sqc->mean[i][k];
        double v = sqc->var[i][k];
        double d = x[k]-m;
        mean[i][k] += sqc->weight[i]*pxi*4*d/(d*d+v)/px;
        var[i][k] += sqc->weight[i]*pxi*(3*d*d-v)/(2*v*(d*d+v))/px;
      }
    }
    count++;
  }
  double update_weight(int i) {
    double grad = weight[i]/count;
    adam_weight->update(i,grad);
    double u = adam_weight->optGradient(i,grad);
    //cout << "Weight " << i << ": grad="<<u<<endl;
    return u;
  }
  double update_mean(int i, int k) {
    double grad = mean[i][k]/count;
    adam_mean[i]->update(k,grad);
    double u = adam_mean[i]->optGradient(k,grad);
    //cout << "Mean " << i << "," << k << ": grad="<< u << endl;
    return u;
  }
  double update_var(int i, int k) {
    double grad = var[i][k]/count;
    adam_var[i]->update(k,grad);
    double u = adam_var[i]->optGradient(k,grad);
    //cout << "Var  " << i << "," << k << ": grad="<< u << endl;
    return u;
  }
};


void read_feature(float *x, int dim, istream &f) {
  f.read((char*)x,dim*sizeof(float));
}

int main(int argc, char *argv[])
{
  int i = 1;
  int dim = 13;
  int mix = 4;
  double delta = 1;
  double var_limit = 1e-2;
  int iter = 10;
  string initfile;
  string datafile;
  string outfile;
  while (i < argc && argv[i][0] == '-') {
    string arg(argv[i]);
    if (arg == "-l") {
      dim = atoi(argv[++i]);
    } else if (arg == "-m") {
      mix = atoi(argv[++i]);
    } else if (arg == "-g") {
      initfile = argv[++i];
    } else if (arg == "-i") {
      datafile = argv[++i];
    } else if (arg == "-o") {
      outfile = argv[++i];
    } else if (arg == "-D") {
      delta = atof(argv[++i]);
    } else if (arg == "-I") {
      iter = atoi(argv[++i]);
    } else if (arg == "-L") {
      var_limit = atof(argv[++i]);
    } else {
      cerr << "Unknown option " << arg << endl;
      exit(1);
    }
    i++;
  }
  if (initfile == "" || datafile == "" || outfile == "") {
    cerr << "initfile=" << initfile << endl;
    cerr << "datafile=" << datafile << endl;
    cerr << "outfile =" << outfile  << endl;
    cerr << "usage: SQCupdate -l dim -g initfile -o outfile -i datafile [-I iter] [-L var_limit]" << endl;
    exit(2);
  }
  SQCauchy sqc(dim,mix);
  sqc.readBinary(initfile);
  sqc.set_var_limit(var_limit);
  ifstream ifile;
  SQCUpdate update(&sqc);
  double lprob,plprob = -1e100;
  for (int t = 0; t < iter; t++) {
    ifile.open(datafile.c_str());
    lprob = 0;
    update.clear();
    while (!ifile.eof()) {
      float x[dim];
      read_feature(x,dim,ifile);
      double lp = sqc.logprob(x);
      //cout << lp << endl;
      lprob += lp;
      update.feed1(x);
    }
    cout << t << ": " << lprob << endl;
    ifile.close();
    if (lprob < plprob || isnan(lprob))  // problem calculating probability
      break;
    else
      sqc.writeBinary(outfile);
    plprob = lprob;
    for (int i = 0; i < mix; i++) {
      sqc.update_weight(i,update.update_weight(i),delta);
      for (int k = 0; k < dim; k++) {
	sqc.update_mean(i,k,update.update_mean(i,k),delta);
	sqc.update_var(i,k,update.update_var(i,k),delta);
      }
    }
  }
  if (!isnan(lprob) && lprob > plprob)
    sqc.writeBinary(outfile);
  return 0;
}
