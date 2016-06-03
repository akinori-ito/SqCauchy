#include "mixture.h"
#include <iostream>
#include <cstring>
#include <cstdio>

using namespace std;
  
void ALGaussianMixture::set_lweight() 
 {
    lweight = new double[nmixture];
    for (int i = 0; i < nmixture; i++) {
      lweight[i] = log(weight[i]);
    }
  }

#if 0
#include <cmath>
#include <random>

int main(int argc, char *argv[])
{
  default_random_engine gen;
  uniform_real_distribution<double> dis(0.0,100.0);

  for (int i = 0; i < 100; i++) {
    double x = exp(-dis(gen));
    double y1 = log(x)/log(2.0);
    double y2 = AddLog::Qlog2(x);
    cout << x << "\t" << y1 << "\t" << y2 << "\t" << y1-y2 << endl;
  }
return 0;
}
#endif
