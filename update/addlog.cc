using namespace std;

#include <cstring>
#include "addlog.h"

double AddLog::addlogtbl[AddLog::tblsize] = {
#include "addlogtbl.d"
};
double AddLog::log2tbl[AddLog::tblsize] = {
#include "addlogtbl2.d"
};

// IEEE double
//
// seeeeeee eeeemmmm mmmmmmmm mmmmmmmm mmmmmmmm mmmmmmmm mmmmmmmm mmmmmmmm
double AddLog::Qlog2(double x) {
  unsigned char p[8];
  memcpy(p, (const void*)&x, 8);
  //  printf("%02x%02x%02x%02x\n",p[0],p[1],p[2],p[3]);
  double ex = ((int)(p[7]&0x7f)<<4)+((int)(p[6]&0xf0)>>4)-1023;
  p[7] = 0x3f;
  p[6] |= 0xf0;
  double mant = *(double*)p-1.0;
  mant *= tick;
  //cout << "mantissa "<<mant<<endl;
  if (mant >= AddLog::tblsize)
    return ex;
  return ex+log2tbl[(int)mant];
}

#ifdef TEST
#include <iostream>
#include <cmath>
extern "C" {
#include <unistd.h>
#include <time.h>
#include <sys/times.h>
}

int main()
{
  const int NMAX = 100000000;
  struct tms tbuf1,tbuf2;
  double x;

  x = 0;
  times(&tbuf1);
  for (int i = 1; i <= NMAX; i++) {
    x += log((double)i);
  }
  times(&tbuf2);
  cout << x << endl;
  cout << "Time: " << (double)(tbuf2.tms_utime - tbuf1.tms_utime)/sysconf(_SC_CLK_TCK) << endl;

  x = 0;
  times(&tbuf1);
  for (int i = 1; i <= NMAX; i++) {
    x += AddLog::Qlog2((double)i);
  }
  times(&tbuf2);
  cout << x*log(2.0) << endl;
  cout << "Time: " << (double)(tbuf2.tms_utime - tbuf1.tms_utime)/sysconf(_SC_CLK_TCK) << endl;

  return 0;
}
#endif
