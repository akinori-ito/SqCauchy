#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cmath>

using namespace std;

#include "SQC_train.h"


SQCauchy::SQCauchy(int Dim, int NMix) {
  dim = Dim;
  nmix = NMix;
  mean = new float*[nmix];
  var = new float*[nmix];
  weight = new float[nmix];
  for (int i = 0; i < nmix; i++) {
     mean[i] = new float[dim];
     var[i] = new float[dim];
  }
  var_limit = 1e-2;
}
SQCauchy::~SQCauchy() {
  for (int i = 0; i < nmix; i++) {
    delete[] mean[i];
    delete[] var[i];
  }
  delete[] mean;
  delete[] var;
  delete[] weight;
}
void SQCauchy::readBinary(string filename) {
  ifstream f;
  f.open(filename.c_str(),ifstream::in);
  if (!f) {
    cerr << "Can't open " << filename << endl;
    exit(1);
  }
  f.read((char*)weight,nmix*sizeof(float));
  for (int i = 0; i < nmix; i++) {
    f.read((char*)mean[i],dim*sizeof(float));
    f.read((char*)var[i],dim*sizeof(float));
  }
  f.close();
}
void SQCauchy::writeBinary(string filename) {
  ofstream f(filename.c_str());
  if (!f) {
    cerr << "Can't open " << filename << endl;
    exit(1);
  }
  f.write((char*)weight,nmix*sizeof(float));
  for (int i = 0; i < nmix; i++) {
    f.write((char*)mean[i],dim*sizeof(float));
    f.write((char*)var[i],dim*sizeof(float));
  }
  f.close();
}

double SQCauchy::prob_i(float *x,int i) {
  double y = 1.0;
  double z;
  for (int k = 0; k < dim; k++) {
    double denom = x[k]-mean[i][k];
    denom = var[i][k]+denom*denom;
    //    cout << "x[" << k << "]=" << x[k] << " ";
    //    cout << "mean["<<i<<"]["<<k<<"]="<<mean[i][k]<<" ";
    //    cout << "var["<<i<<"]["<<k<<"]="<<var[i][k]<<" ";
    z = 2*var[i][k]*sqrt(var[i][k])/(3.14159265*denom*denom);
    //    cout << "prob=" << z << endl;
    y *= z;
  }
  return y;
}

double SQCauchy::prob(float *x) {
  double y = 0;
  for (int i = 0; i < nmix; i++) {
    y += weight[i]*prob_i(x,i);
  }
  return y;
}

double SQCauchy::logprob(float *x) {
  return log(prob(x));
}

void SQCauchy::update_weight(int i, double diff, double delta) {
  //std::cout << "Weight " << i << " delta=" << delta*diff <<" " << weight[i];
  weight[i] += delta*diff;
  //std::cout << " => " << weight[i] << std::endl;
}
void SQCauchy::update_mean(int i, int k, double diff, double delta) {
  //std::cout << "Mean " << i << "," << k << " delta=" << delta*diff <<" " << mean[i][k];
  mean[i][k] += delta*diff;
  //std::cout << " => " << mean[i][k] << std::endl;
}
void SQCauchy::update_var(int i, int k, double diff, double delta) {
  //std::cout << "Var  " << i << "," << k << " delta=" << delta*diff <<" " << var[i][k];
  double prev = var[i][k];
  var[i][k] += delta*diff;
  if (var[i][k] <= var_limit) {
    //cout << "var(" << i << "," << k << ") become small or non-positive: unchanged" << endl;
    var[i][k] = prev;
  }
  //std::cout << " => " << var[i][k] << std::endl;
}
