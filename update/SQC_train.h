#ifndef SQC_TRAIN_H
#define SQC_TRAIN_H
#include <string>
class SQCauchy {
 public:
  int dim;
  int nmix;
  float **mean;
  float **var;
  float *weight;
  double var_limit;
  SQCauchy(int Dim, int NMix);
  ~SQCauchy();
  void readBinary(std::string filename);
  void writeBinary(std::string filename);
  double prob_i(float *x,int i);
  double prob(float *x);
  double logprob(float *x);
  void update_weight(int i, double diff, double delta = 1.0);
  void update_mean(int i, int k, double diff, double delta = 1.0);
  void update_var(int i, int k, double diff, double delta = 1.0);
  void set_var_limit(double x) {
    var_limit = x;
  }
};
#endif
