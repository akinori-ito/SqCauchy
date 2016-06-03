class Adam {
  int dim;
  double *m;
  double *v;
  double alpha;
  double beta1;
  double beta2;
  double beta1t;
  double beta2t;
  double epsilon;
public:
  Adam(int Dim) {
    dim = Dim;
    m = new double[dim];
    v = new double[dim];
    alpha = 0.001;
    beta1 = 0.9;
    beta2 = 0.999;
    beta1t = beta1;
    beta2t = beta2;
    epsilon = 1e-8;
    for (int i = 0; i < dim; i++) {
      m[i] = 0;
      v[i] = 0;
    }
  }
  ~Adam() {
    delete[] m;
    delete[] v;
  }
  void update(int d, double grad) {
    m[d] = beta1*m[d]+(1-beta1)*grad;
    v[d] = beta2*v[d]+(1-beta2)*grad*grad;
  }
  double optGradient(int d, double grad) {
    double hm = m[d]/(1-beta1t);
    double hv = v[d]/(1-beta2t);
    beta1t = beta1t*beta1;
    beta2t = beta2t*beta2;
    return alpha*hm/(sqrt(hv)+epsilon);
  }
};
