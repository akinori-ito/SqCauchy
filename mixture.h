#ifndef MIXTURE_H
#define MIXTURE_H

#include <cstring>
#include "distribution.h"
#include "addlog.h"


class MixtureDensity {
 protected:
  int nmixture;  // number of mixture
  double *weight;
  DiagonalMultivariateDistribution **dist;
 public:
  MixtureDensity(int n) {
    nmixture = n;
    weight = new double[n];
    dist = new DiagonalMultivariateDistribution*[n];
    for (int i = 0; i < nmixture; i++)
      dist[i] = 0;
  }
  virtual ~MixtureDensity() {
    for (int i = 0; i < nmixture; i++)
      if (dist[i] != 0)
	delete dist[i];
    delete[] dist;
  }
  virtual void ready() {
    for (int i = 0; i < nmixture; i++)
      dist[i]->ready();
  }
  void set(int i, double w, DiagonalMultivariateDistribution *d) {
    if (dist[i] != 0)
      delete dist[i];
    weight[i] = w;
    dist[i] = d;
  }
  virtual void setDistribution(int i, int dim, double w, double *mean, double *var) {}
  virtual void setDistribution(int i, int dim, float w, float *mean, float *var) {}
  virtual void setParam(const char *param, double val) {};
  virtual double prob(double *x) {
    double p = 0.0;
    for (int i = 0; i < nmixture; i++) {
//      if (dist[i] != 0)
	p += weight[i]*dist[i]->prob(x);
    }
    return p;
  }
  virtual double lprob(double *x) {
    //return log(prob(x));
    return AddLog::Qlog2(prob(x));
  }
};

 
class GaussianMixture: public MixtureDensity {
 public:
  GaussianMixture(int n): MixtureDensity(n) {}
  void setDistribution(int i, int dim, double w, double *mean, double *var) {
    set(i, w, (DiagonalMultivariateDistribution*)new DiagonalGaussian(dim,mean,var));
  }
  void setDistribution(int i, int dim, float w, float *mean, float *var) {
    set(i, w, (DiagonalMultivariateDistribution*)new DiagonalGaussian(dim,mean,var));
  }
};

// Class of Gaussian mixture which used addlog operation for probability calculation.
class ALGaussianMixture: public GaussianMixture {
  double *lweight;
  void set_lweight();
 public:
  ALGaussianMixture(int n): GaussianMixture(n) {
    lweight = 0;
  }
  void ready () {
    GaussianMixture::ready();
    set_lweight();
  }
  double lprob(double *x) {
    double p = -1.0e30;
    for (int i = 0; i < nmixture; i++) {
      p = AddLog::addlog(p,dist[i]->lprob(x)+lweight[i]);
    }
    return p;
  }
  double prob(double *x) {
    return exp(lprob(x));
  }
};

class SqCauchyMixture: public MixtureDensity {
  double varianceMagnitude;
 public:
  SqCauchyMixture(int n): MixtureDensity(n) {
    varianceMagnitude = 1.0;
  }
  void setParam(const char *param, double x) {
    if (strcmp(param, "varMagnitude") == 0)
      varianceMagnitude = x;
  }
  void setDistribution(int i, int dim, double w, double *mean, double *var) {
    set(i, w, new DiagonalSqCauchy(dim,mean,var,varianceMagnitude));
  }
  void setDistribution(int i, int dim, float w, float *mean, float *var) {
    set(i, w, new DiagonalSqCauchy(dim,mean,var,varianceMagnitude));
  }
  void set(int i, double w, DiagonalSqCauchy *d) {
    if (dist[i] != 0)
      delete dist[i];
    weight[i] = w; // not needed
    dist[i] = (DiagonalMultivariateDistribution*)d;
    d->multiplyWeight(w);
  } 
  double prob(double *x) {
    double p = 0.0;
    for (int i = 0; i < nmixture; i++) {
      p += dist[i]->prob(x);
    }
    return p;
  }
};

class CauchyMixture: public MixtureDensity {
  double varianceMagnitude;
 public:
  CauchyMixture(int n): MixtureDensity(n) {
    varianceMagnitude = 1.0;
  }
  void setParam(const char *param, double x) {
    if (strcmp(param, "varMagnitude") == 0)
      varianceMagnitude = x;
  }
  void setDistribution(int i, int dim, double w, double *mean, double *var) {
    set(i, w, (DiagonalMultivariateDistribution*)new DiagonalCauchy(dim,mean,var,varianceMagnitude));
  }
  void setDistribution(int i, int dim, float w, float *mean, float *var) {
    set(i, w, (DiagonalMultivariateDistribution*)new DiagonalCauchy(dim,mean,var,varianceMagnitude));
  }
};

#endif    
