#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H
#include <cmath>
#include <vector>
#include "addlog.h"

#define RAPIDCAUCHY

class Distribution {
  public:
  virtual double prob(double x) = 0;
  virtual double lprob(double x) = 0;
};

class GaussianDistribution: public Distribution {
friend class DiagonalGaussian;
  private:
    double mean;
    double var2;
    double rvar2;
    double coef;
    double lcoef;
  public:
    GaussianDistribution(double m, double v)  {
        mean = m;
        var2 = 2*v;
        rvar2 = 1.0/var2;
        coef = 1.0/sqrt(2*3.14159265*v);
	lcoef = log(coef);
    }
    double prob(double x) {
      double r = x-mean;
      return coef*exp(-r*r*rvar2);
    }
    double lprob(double x) {
      double r = x-mean;
      return lcoef-r*r*var2;
    }
};

class CauchyDistribution: public Distribution {
  friend class DiagonalCauchy;
  private:
    double mean;
    double s2;
    double s1;
    static const double coef = 1.0/3.14159265;
  public:
    CauchyDistribution(double m, double v) {
        mean = m;
        s2 = v;
        s1 = sqrt(v)*coef;
    }
    double prob(double x) {
        double r = x-mean;
        r = s2+r*r;
        return s1/(r*r);
    }
    double lprob(double x) {
      return log(prob(x));
    }
};

class SquareCauchyDistribution: public Distribution {
  friend class DiagonalSqCauchy;
  private:
    double mean;
    double s2;
    double s3;
    static const double coef = 2.0/3.14159265;
  public:
    SquareCauchyDistribution(double m, double v) {
        mean = m;
        s2 = v;
        s3 = s2*sqrt(v)*coef;
    }
    double prob(double x) {
        double r = x-mean;
        r = s2+r*r;
        return s3/(r*r);
    }
    double lprob(double x) {
      return log(prob(x));
    }
};

class DiagonalMultivariateDistribution {
 protected:
  int dimension;
  Distribution **dist;
 public:
  DiagonalMultivariateDistribution(int n) {
    dimension = n;
    dist = new Distribution*[n];
  }
  virtual ~DiagonalMultivariateDistribution() {}
  void setDistribution(int i, Distribution *d) {
    dist[i] = d;
  }
  virtual void ready() {}
  virtual double prob(double *x) {
    double p = 1.0;
    for (int i = 0; i < dimension; i++) {
      p *=  dist[i]->prob(x[i]);
    }
    return p;
  }
  virtual double prob(std::vector<double> x) {
    double p = 1.0;
    for (int i = 0; i < dimension; i++) {
      p *=  dist[i]->prob(x[i]);
    }
    return p;
  }
  virtual double lprob(double *x) {
    double p = 0.0;
    for (int i = 0; i < dimension; i++) {
      p += dist[i]->lprob(x[i]);
    }
    return p;
  }
  virtual double lprob(std::vector<double> x) {
    double p = 0.0;
    for (int i = 0; i < dimension; i++) {
      p += dist[i]->lprob(x[i]);
    }
    return p;
  }
};

class DiagonalGaussian: public DiagonalMultivariateDistribution {
  double coef;
  double lcoef;
 
  void setCoef() {
    coef = 1.0;
    for (int i = 0; i < dimension; i++) {
      coef *= ((GaussianDistribution*)dist[i])->coef;
    }
    lcoef = log(coef);
  }
 public:
  DiagonalGaussian(int n) : DiagonalMultivariateDistribution(n) {}
  DiagonalGaussian(int n, double *means, double *vars) : DiagonalMultivariateDistribution(n) {
    for (int i = 0; i < dimension; i++) {
      setDistribution(i,means[i],vars[i]);
    }
    setCoef();
  }
  DiagonalGaussian(int n, float *means, float *vars) : DiagonalMultivariateDistribution(n) {
    for (int i = 0; i < dimension; i++) {
      setDistribution(i,(double)means[i],(double)vars[i]);
    }
    setCoef();
  }
  DiagonalGaussian(int n, std::vector<double> means, std::vector<double> vars) : DiagonalMultivariateDistribution(n) {
    for (int i = 0; i < dimension; i++) {
      setDistribution(i,means[i],vars[i]);
    }
    setCoef();
  }
  void setDistribution(int i, double mean, double var) {
    dist[i] = (Distribution*)new GaussianDistribution(mean,var);
  }
  double prob(double *x) {
    double p = 1.0;
    for (int i = 0; i < dimension; i++) {
      double r = x[i]-((GaussianDistribution*)dist[i])->mean;
      p *= exp(-r*r*((GaussianDistribution*)dist[i])->rvar2);
    }
    return p*coef;
  }
  double lprob(double *x) {
    double lp = 0.0;
    for (int i = 0; i < dimension; i++) {
      double r = x[i]-((GaussianDistribution*)dist[i])->mean;
      lp += -r*r*((GaussianDistribution*)dist[i])->rvar2;
    }
    return lp+lcoef;
  }
      
};

class DiagonalSqCauchy: public DiagonalMultivariateDistribution {
#ifdef RAPIDCAUCHY
  double numer;
#endif
 public:
  DiagonalSqCauchy(int n) : DiagonalMultivariateDistribution(n) {}
  DiagonalSqCauchy(int n, float *means, float *vars, double mag = 1.0) : DiagonalMultivariateDistribution(n) {
    for (int i = 0; i < dimension; i++) {
      setDistribution(i,(double)means[i],(double)vars[i]*mag);
    }
  }
 DiagonalSqCauchy(int n, double *means, double *vars, double mag = 1.0) : DiagonalMultivariateDistribution(n) {
    for (int i = 0; i < dimension; i++) {
      setDistribution(i,means[i],vars[i]*mag);
    }
  }
  DiagonalSqCauchy(int n, std::vector<double> means, std::vector<double> vars, double mag = 1.0) : DiagonalMultivariateDistribution(n) {
    for (int i = 0; i < dimension; i++) {
      setDistribution(i,means[i],vars[i]*mag);
    }
  }
  void setDistribution(int i, double mean, double var, double mag=1.0) {
    dist[i] = (Distribution*)new SquareCauchyDistribution(mean,var*mag);
  }
#ifdef RAPIDCAUCHY
  void ready() {
    numer = 1.0;
    for (int i = 0; i < dimension; i++) {
     SquareCauchyDistribution *sq = (SquareCauchyDistribution*)dist[i];
      numer *= sq->s3;
    }
  }
  void multiplyWeight(double weight) {
    numer *= weight;
  }
  double prob(double *x) {
    double denom = 1.0;
    SquareCauchyDistribution **sq = (SquareCauchyDistribution**)dist;
    for (int i = 0; i < dimension; i++) {
      double r = x[i]-(*sq)->mean;
      denom *= (*sq)->s2+r*r;
      sq++;
    }
    return numer/(denom*denom);
  }
/*
  double lprob(double *x) {
    return AddLog::Qlog2(prob(x));
  }
*/
#endif

};

class DiagonalCauchy: public DiagonalMultivariateDistribution {
  double numer;
 public:
  DiagonalCauchy(int n) : DiagonalMultivariateDistribution(n) {}
  DiagonalCauchy(int n, float *means, float *vars, double mag = 1.0) : DiagonalMultivariateDistribution(n) {
    for (int i = 0; i < dimension; i++) {
      setDistribution(i,(double)means[i],(double)vars[i]*mag);
    }
  }
 DiagonalCauchy(int n, double *means, double *vars, double mag = 1.0) : DiagonalMultivariateDistribution(n) {
    for (int i = 0; i < dimension; i++) {
      setDistribution(i,means[i],vars[i]*mag);
    }
  }
  DiagonalCauchy(int n, std::vector<double> means, std::vector<double> vars, double mag = 1.0) : DiagonalMultivariateDistribution(n) {
    for (int i = 0; i < dimension; i++) {
      setDistribution(i,means[i],vars[i]*mag);
    }
  }
  void setDistribution(int i, double mean, double var, double mag=1.0) {
    dist[i] = (Distribution*)new SquareCauchyDistribution(mean,var*mag);
  }
  void ready() {
    numer = 1.0;
    for (int i = 0; i < dimension; i++) {
      CauchyDistribution *cd = (CauchyDistribution*)dist[i];
      numer *= cd->s1;
    }
  }
  double prob(double *x) {
    double denom = 1.0;
    for (int i = 0; i < dimension; i++) {
      CauchyDistribution *sq = (CauchyDistribution*)dist[i];
      double r = x[i]-sq->mean;
      denom *= sq->s2+r*r;
    }
    return numer/denom;
  }
};

#endif
