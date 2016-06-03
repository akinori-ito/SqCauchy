#ifndef ADDLOG_H
#define ADDLOG_H

class AddLog {
  static const int tblsize = 1000;
  static const int tick = 100;
  static double addlogtbl[tblsize];
  static double log2tbl[tblsize];
 public:
  // calculate log(x+exp(-x)) using table
  static double t_log(double x) {
    //if (x < 0) return 0;  x should be non-negative
    x = x*tick;
    if (x >= tblsize) return 0;
    return addlogtbl[(int)x];
  }
  // calculate log(exp(x)+exp(y)) using table
  static double addlog(double x, double y) {
    double a,b;
    if (x > y) {
      a = x; b = y;
    } else {
      a = y; b = x;
    }
    return a+t_log(a-b);
  }
  static double Qlog2(double x);
};

#endif