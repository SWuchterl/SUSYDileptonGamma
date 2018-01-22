#ifndef MT2FUNCTOR_H_
#define MT2FUNCTOR_H_

#include <iostream>
#include <math.h>
#include "TObject.h"

// using namespace std;

class MT2Functor{
  // class MT2Functor : public TObject {
 public:
  static const float RELATIVE_PRECISION; 
  static const float ABSOLUTE_PRECISION;
  static const float MIN_MASS;
  static const float ZERO_MASS;
  static const float SCANSTEP;
  
  MT2Functor();
  virtual ~MT2Functor();
  void   mt2_bisect();
  void   mt2_massless();
  void   set_momenta(double *pa0, double *pb0, double* pmiss0, unsigned long evtNr);
  void   set_mn(double mn);
  //~ inline void set_verbose(int vlevel){verbose = vlevel;};
  double get_mt2();
  void   print();
  unsigned long    nevt;

 private:

  int verbose;
  bool   solved;
  bool   momenta_set;
  double mt2_b;

  int    nsols(double Dsq);
  int    nsols_massless(double Dsq);
  //inline
  int    signchange_n( long double t1, long double t2, long double t3, long double t4, long double t5);
  //inline
  int    signchange_p( long double t1, long double t2, long double t3, long double t4, long double t5);
  //~ int scan_high(double &Deltasq_high);
  int find_high(double &Deltasq_high);
  //data members
  double pax, pay, ma, Ea;
  double pmissx, pmissy;
  double pbx, pby, mb, Eb;
  double mn, mn_unscale;

  //auxiliary definitions
  double masq, Easq;
  double mbsq, Ebsq;
  double pmissxsq, pmissysq;
  double mnsq;

  //auxiliary coefficients
  double a1, b1, c1, a2, b2, c2, d1, e1, f1, d2, e2, f2;
  double d11, e11, f12, f10, d21, d20, e21, e20, f22, f21, f20;

  double scale;
  double precision;
  // ClassDef(MT2Functor,1)
};

#endif
