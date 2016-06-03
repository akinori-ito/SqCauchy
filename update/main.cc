#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>
#include "distribution.h"
#include "mixture.h"

//#define DEBUG

extern "C" {
#include <unistd.h>
#include <time.h>
#include <sys/times.h>
#include <sys/resource.h>
}

using namespace std;

int DIM = 13;
const int MAXFRAME = 10000;

enum DistributionName {Gaussian, AddlogGaussian, Cauchy, SquareCauchy};

/*
GaussianMixture *readGMM_ascii(string filename) {
  ifstream f;
  f.open(filename.c_str(),ios::in | ios::binary);
  if (!f) {
    cerr << "Can't open " << filename << endl;
    exit(2);
  }
  int nmix;
  f >> nmix;
  GaussianMixture *gmm = new GaussianMixture(nmix);
  double weight[nmix];
  double mean[nmix][DIM];
  double var[nmix][DIM];
  int i,j;
  for (i = 0; i < nmix; i++)
    f >> weight[i];
  for (i = 0; i < nmix; i++) {
    for (j = 0; j < DIM; j++) {
      f >> mean[i][j];
    }    
  }
  for (i = 0; i < nmix; i++) {
    for (j = 0; j < DIM; j++) {
      f >> var[i][j];
    }    
  }
  for (i = 0; i < nmix; i++) {
    for (j = 0; j < DIM; j++) {
      gmm->setDistribution(i,DIM,weight[i],mean[i],var[i]);
    }    
  }
  return gmm;

}
*/

MixtureDensity *readGMM_binary(int nmix, string filename, 
      enum DistributionName distname = Gaussian, double mag = 0.0) {
  ifstream f;
  f.open(filename.c_str());
  if (!f) {
    cerr << "Can't open " << filename << endl;
    exit(2);
  }
  MixtureDensity *gmm = 0;
  switch (distname) {
	case Gaussian:
		gmm = (MixtureDensity*)new GaussianMixture(nmix);
		break;
	case AddlogGaussian:
		gmm = (MixtureDensity*)new ALGaussianMixture(nmix);
		break;
	case Cauchy:
		gmm = (MixtureDensity*)new CauchyMixture(nmix);
		gmm->setParam("varMagnitude",mag);
		break;
	case SquareCauchy:
		gmm = (MixtureDensity*)new SqCauchyMixture(nmix);
		gmm->setParam("varMagnitude",mag);
		break;
  }

  float weight[nmix];
  float mean[nmix][DIM];
  float var[nmix][DIM];
  int i,j;
  f.read((char*)weight,nmix*sizeof(float));
  for (i = 0; i < nmix; i++) {
    f.read((char*)mean[i],DIM*sizeof(float));
    f.read((char*)var[i],DIM*sizeof(float));
  }
  for (i = 0; i < nmix; i++) {
    for (j = 0; j < DIM; j++) {
      gmm->setDistribution(i,DIM,weight[i],mean[i],var[i]);
    }    
  }
  gmm->ready();
  return gmm;
}

void readFeature_ascii(ifstream& f, double *vec)
{
  for (int i = 0; i < DIM; i++) {
    f >> vec[i];
  }
}

bool readFeature(ifstream& f, double *vec)
{
  float x;
  for (int i = 0; i < DIM; i++) {
    f.read((char*)&x,sizeof(float));
    if (f.eof())
	return true;
    vec[i] = x;
  }
  return false;
}

struct likelihood {
  int n;
  double lprob;
};

bool compare(const likelihood& x, const likelihood &y) {
  return y.lprob < x.lprob;
}

double timeval2sec(struct timeval &tv)
{
  return (double)tv.tv_sec+(double)tv.tv_usec/1000000;
}

int main(int argc, char *argv[])
{
  int i;
  int nmix = 64;
  string gmmlist = "files-gmm.txt";
  double mag = 0.0;
  enum DistributionName distname = Gaussian;

  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-')
	break;
    if (strcmp(argv[i],"-m") == 0)
	nmix = atoi(argv[++i]);
    else if (strcmp(argv[i],"-l") == 0)
	DIM = atoi(argv[++i]);
    else if (strcmp(argv[i],"-g") == 0)
	gmmlist = argv[++i];
    else if (strcmp(argv[i],"-cauchy") == 0) {
		distname = Cauchy;
        mag = atof(argv[++i]);
    }
    else if (strcmp(argv[i],"-sqcauchy") == 0) {
		distname = SquareCauchy;
        mag = atof(argv[++i]);
    }
    else if (strcmp(argv[i],"-gauss") == 0)
        distname = Gaussian;
    else if (strcmp(argv[i],"-addlog") == 0)
        distname = AddlogGaussian;
    else {
        cerr << "Unknown option " << argv[i] << endl;
        exit(1);
    }
  }
  if (i >= argc) {
    cerr << "Usage: a.out [-m mix] [-l dim] [-g gmmlist] [-cauchy mag] [-gauss] [-addlog] file.mfcc" << endl;
    exit(1);
  }
  char *testfile = argv[i];
  ifstream istr;
  istr.open(gmmlist.c_str(),ifstream::in);
  vector<string> fn;
  vector<MixtureDensity*> gmm;
  i = 0;
  do {
    string x;
    istr >> x;
    if (x.length() == 0)
      continue;
    fn.push_back(x);
    gmm.push_back(readGMM_binary(nmix,fn[i],distname,mag));
    //cout << i << ": " << fn[i] << endl;
    i++;
  } while (!istr.eof());
  // cout << "Number of file = " << fn.size() << endl;
  istr.close();

  int nmodel = fn.size();

  vector<likelihood> lprobs(nmodel);
  for (int i = 0; i < nmodel; i++) {
    lprobs[i].n = i;
    lprobs[i].lprob = 0;
  }
  istr.open(testfile,ios::in|ios::binary);
  int nframe = 0;
  double total_tick = 0;
  double *frame[MAXFRAME];
  struct rusage tbuf1,tbuf2;
  while (nframe < MAXFRAME) {
    frame[nframe] = (double*)malloc(sizeof(double)*DIM);
    if (readFeature(istr,frame[nframe]))
	break;
    nframe++;
  }
  if (nframe == MAXFRAME) {
    cerr << "Error: Signal too long" << endl;
    exit(1);
  }
  getrusage(RUSAGE_SELF,&tbuf1);
  for (int f = 0; f < nframe; f++) {
    for (int i = 0; i < nmodel; i++) {
      lprobs[i].lprob += gmm[i]->lprob(frame[f]);
#ifdef DEBUG
       cout << f << ": lprobs[" << fn[i] << "]=" << lprobs[i].lprob << endl;
#endif
    }
  }
  getrusage(RUSAGE_SELF,&tbuf2);
  total_tick = timeval2sec(tbuf2.ru_utime)-timeval2sec(tbuf1.ru_utime);
  sort(lprobs.begin(), lprobs.end(), compare);

#ifdef DEBUG
  for (int i = 0; i < nmodel; i++) {
    cout << lprobs[i].n << " " << lprobs[i].lprob << endl;
  }
#else
  cout << lprobs[0].n << " " << lprobs[0].lprob << endl;
#endif
  cout << "Time: " << total_tick<< endl;


  return 0;
}
