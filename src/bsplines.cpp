#include "bsplines.h"

int bsplines(int n, int k, 
             std::vector<double> &gl_x,
             std::vector<double> &knots,
             std::vector<double> &splines) {
  int i1 ;
  double dl, sl, x;
  int len = (k+1)*(k+2)/2;
  std::vector<double> Db(k);
  std::vector<double> work(len);

  for(auto i=k-1; i<n; ++i){
    dl = knots[i+1] - knots[i];
    sl = knots[i+1] + knots[i];

    for(int p=0; p<k; ++p){
      x = dl*0.5 * gl_x[p] + sl*0.5;    //x-transformation
      i1 = i + 1 ;
      dbspvd_(&knots[0], k, 1, x, i1, k, &Db[0], &work[0]);

      splines.insert(std::end(splines), std::begin(Db), std::end(Db));
    }
  }
  return 0;
}

int bsplinesp(int n, int k, 
              std::vector<double> &gl_x,
              std::vector<double> &knots,
              std::vector<double> &splinesp) {
  int i1 ;
  double dl, sl, x;
  int len = (k+1)*(k+2)/2;
  std::vector<double> Db(k*2);
  std::vector<double> work(len);

  for(auto i=k-1; i<n; ++i){
    dl = knots[i+1] - knots[i];
    sl = knots[i+1] + knots[i];

    for(int p=0; p<k; ++p){
      x = dl*0.5 * gl_x[p] + sl*0.5;    //x-transformation
      i1 = i + 1 ;
      dbspvd_(&knots[0], k, 2, x, i1, k, &Db[0], &work[0]);

      splinesp.insert(std::end(splinesp), std::begin(Db)+k, std::end(Db));
    }
  }
  return 0;
}

int bspline_int(int n, int k,
                std::vector<double> &gl_w,
                std::vector<double> &gl_x,
                std::vector<double> &ov,
                int di, int dj,
                std::vector<double> &bsp,
                ModelV &V) {
  if (di==0 && dj==0) {
    for() {

      for() {

        for(int p=0; p<k; ++p) {

        }
      } 
    }
  } else if (di==1 && dj==1) {
    for() {

      for() {

        for(int p=0; p<k; ++p) {
        
        }
      } 
    }
  } else {  

  }  
  return 0;
}