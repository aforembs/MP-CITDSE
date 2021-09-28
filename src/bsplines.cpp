#include "bsplines.h"

int bsp::GenKnots(int n, int k, double r_max, 
                  double fkn, char type,
                  std::vector<double> &kkn) {
  double par=0.0;
  static constexpr double PI = 3.141592653589793238463;
  kkn.reserve(n+k);
  for(int i=0; i<k; ++i) kkn.emplace_back(0.0);

  if(type=='l') {
    par = r_max/(double)(n-k+1);

    for(int i=k; i<n; ++i)
      kkn.emplace_back{par};

  } else if(type=='e') {
    par = log(r_max/fkn)/(double)(n-k+1);

    for(int i=k; i<n; ++i)
      kkn.emplace_back(fkn*exp((double)(i-k_)*par));

  } else if(type=='s') {
    par = -log((2./PI)*asin(fkn/r_max))/log((double)(n-k+1));

    for(int i=k; i<n; ++i) 
      kkn.emplace_back(r_max*sin((PI/2.)
      *pow((double)(i-k_+1)/(double)(n_-k_+1),par)));

  } else {
    std::cout << "Invalid knot type !\n";
    return 1;
  }
  for(int i=n; i<n+k; ++i) {
    kkn.emplace_back(r_max);
  }
  return 0;
}

int bsp::GenKnots(int n, int k, 
                  std::string file, char type, 
                  std::vector<double> &kkn) {
  double pt=0.0;
  kkn.reserve(n+k);
  for(int i=0; i<k; ++i) kkn.emplace_back(0.0);

  if(type=='t') {
    std::ifstream fl(file);

    while(fl) {
      fl >> pt;
      kkn.emplace_back(pt);
    }
  } else if(type=='b') {
    std::ifstream fl(file, sdt::ios::in | std::ios::binary);
    fl.unsetf(std::ios::skipws);

    kkn.insert(std::end(kkn),
              std::istream_iterator<double>(file),
              std::istream_iterator<double>());
  } else {
    std::cout << "Invalid knot file type !\n";
    return 1;
  }

  for(int i=n; i<n+k; ++i) {
    kkn.emplace_back(r_max);
  }
  return 0;
}

int WrKnotsH5(int n, int k, double r_max, 
              double fkn, char type,
              std::string file,
              std::vector<double> &kkn) {
  std::string tp_string;

  if(type=='l') {
    tp_string="linear";
  } else if(type=='e') {
    tp_string="exponential";
  } else if(type=='s') {
    tp_string="sine";
  } else if(type=='c') {
    tp_string="custom";
  }

  H5::H5File *fl = new H5::H5File(file, H5F_ACC_TRUNC);
  // Check if the file was opened
  if (!file) {
    cerr << "# H5::H5File:: file couldn't opened: " << outFile.c_str() << "\n";
    exit(-1);
  }

  hsize_t nKnots_d[1] = {kkn.size()};
  hsize_t att_space[1] = {1};

  H5::Attribute N = file->createAttribute("N", H5::PredType::NATIVE_INT32, H5::DataSpace(1, att_space));
  H5::Attribute K = file->createAttribute("K", H5::PredType::NATIVE_INT32, H5::DataSpace(1, att_space));
  H5::Attribute R = file->createAttribute("R", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, att_space));

  return 0;
}

int bsp::Splines(int n, int k, 
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

int bsp::SplinesP(int n, int k, 
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

int bsp::SplineInt(int n, int k,
                  std::vector<double> &gl_w,
                  std::vector<double> &gl_x,
                  std::vector<double> &ov,
                  std::vector<double> &spl,
                  std::vector<double> &kkn,
                  ModelV &Vptr) {
  int j_max, iv_min, iv_max, ivm1;  
  double ovlp, bsum, gsum, dl, sl, x;

  for(int i=0; i<n; ++i) {
    j_max = min(i+k,n);
    for(int j=i; j<j_max; ++j) {
      iv_min = max(k, j+2);
      iv_max = min(i+k+1, n);
      sum_ij=0.0;
      ovlp=0.0;

      for(int iv=iv_min; iv<=iv_max; ++iv) {
        ivm1 = iv-1;
        dl = (kkn[iv] - kkn[ivm1])*0.5;
        sl = (kkn[iv] + kkn[ivm1])*0.5;
        bsum = 0.0;
        for(int p=0; p<k; ++p) {
          x=dl*gl_x[p]+sl;
          bsum += gl_w[p]*spl[i+k*(p+ivm1*k)]*spl[j+k*(p+ivm1*k)]*Vptr.V(x);
        }
        ovlp += dl*bsum;
      }
      // index into ov here needs Fortran (col major)
      ov[j*k+k-1+i-j] = ovlp;
    } 
  }
  return 0;
}

// int bsp::SplineInt(int n, int k,
//                   std::vector<double> &gl_w,
//                   std::vector<double> &gl_x,
//                   std::vector<double> &ov,
//                   int di, int dj,
//                   std::vector<double> &spl,
//                   std::vector<double> &splp,
//                   std::vector<double> &kkn,
//                   ModelV &V) {
//   if (di==1 && dj==0) {
//     for(int i=0; i<n; ++i) {

//       for() {

//         for(int p=0; p<k; ++p) {
        
//         }
//       } 
//     }
//   }  else if (di==0 && dj==1) {
//     for(int i=0; i<n; ++i) {

//       for() {

//         for(int p=0; p<k; ++p) {
        
//         }
//       } 
//     }
//   } else {

//   }
//   return 0;
// }