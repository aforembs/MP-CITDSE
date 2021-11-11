#include "bsplines.h"

int bsp::GenKnots(int n, int k, double r_max, 
                  double fkn, char type,
                  std::vector<double> &kkn) {
  double par=0.0, lin=0.0;
  static constexpr double PI = 3.141592653589793238463;
  kkn.reserve(n+k);
  for(int i=0; i<k; ++i) {kkn.emplace_back(0.0);}

  if(type=='l') {
    par = r_max/(double)(n-k+1);

    for(int i=k; i<n; ++i) {
      lin+=par;
      kkn.emplace_back(lin);
    }

  } else if(type=='e') {
    par = log(r_max/fkn)/(double)(n-k+1);

    for(int i=k; i<n; ++i) {kkn.emplace_back(fkn*exp((double)(i-k)*par));}

  } else if(type=='s') {
    par = -log((2./PI)*asin(fkn/r_max))/log((double)(n-k+1));

    for(int i=k; i<n; ++i) {
      kkn.emplace_back(r_max*sin((PI/2.)
      *pow((double)(i-k+1)/(double)(n-k+1),par)));}

  } else {
    std::cout << "Invalid knot type !\n";
    return 1;
  }
  for(int i=n; i<n+k; ++i) {kkn.emplace_back(r_max);}
  return 0;
}

int bsp::GenKnots(int n, int k, double r_max,
                  std::string file, char type, 
                  std::vector<double> &kkn) {
  std::string pt;
  int line_num=0;
  kkn.reserve(n+k);
  for(int i=0; i<k; ++i) {kkn.emplace_back(0.0);}

  if(type=='t') {
    std::ifstream fl(file);

    while(std::getline(fl, pt)) {
      ++line_num;
      kkn.emplace_back(std::stod(pt));
    }
    //assert(line_num==n);
  } else if(type=='b') {
    std::ifstream fl(file, std::ios::in | std::ios::binary);
    fl.unsetf(std::ios::skipws);

    kkn.insert(std::end(kkn),
              std::istream_iterator<double>(fl),
              std::istream_iterator<double>());
  } else {
    std::cout << "Invalid knot file type !\n";
    return -1;
  }

  for(int i=n; i<n+k; ++i) {kkn.emplace_back(r_max);}
  return 0;
}

int bsp::WrKnotsH5(int n, int k, double r_max, 
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

  auto fl = std::unique_ptr<H5::H5File>(new H5::H5File(file, H5F_ACC_TRUNC));
  // Check if the file was opened
  if (!fl) {
    std::cerr << "# H5::H5File:: file couldn't opened: " << file << "\n";
    exit(-1);
  }

  hsize_t nKnots_d[1]  = {kkn.size()};
  hsize_t att_space[1] = {1};
  hsize_t str_space[1] = {sizeof(tp_string)/sizeof(char)};
  H5::StrType str_ty(H5::PredType::C_S1, H5T_VARIABLE);

  H5::Attribute N = fl->createAttribute("N", H5::PredType::NATIVE_INT32, H5::DataSpace(1, att_space));
  H5::Attribute K = fl->createAttribute("K", H5::PredType::NATIVE_INT32, H5::DataSpace(1, att_space));
  H5::Attribute R = fl->createAttribute("R", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, att_space));
  H5::Attribute TP= fl->createAttribute("Type", str_ty, H5::DataSpace(1, str_space));
  H5::Attribute F = fl->createAttribute("FKN", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, att_space));

  H5::DataSet Knots = fl->createDataSet("Knots", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, nKnots_d));

  N.write(H5::PredType::NATIVE_INT32, &n);
  K.write(H5::PredType::NATIVE_INT32, &k);
  R.write(H5::PredType::NATIVE_DOUBLE, &r_max);
  TP.write(str_ty, &tp_string);
  F.write(H5::PredType::NATIVE_DOUBLE, &fkn);
  Knots.write(&kkn[0], H5::PredType::NATIVE_DOUBLE);

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

  splines.reserve(n*k*k);

  for(int i=0; i<(k-1); ++i) {
    for(int j=0; j<k*k; ++j) { 
      {splines.emplace_back(0.0);} 
    }
  }

  for(auto i=k-1; i<n; ++i){
    dl = knots[i+1] - knots[i];
    sl = knots[i+1] + knots[i];
    i1 = i + 1 ;

    for(int p=0; p<k; ++p){
      x = dl*0.5 * gl_x[p] + sl*0.5;    //x-transformation
      
      dbspvd_(&knots[0], k, 1, x, i1, k, &Db[0], &work[0]);

      splines.insert(std::end(splines), std::begin(Db), std::end(Db));
    }
  }
  return 0;
}

int bsp::SimpSplines(int n, int k, 
                std::vector<double> &gl_x,
                std::vector<double> &knots,
                std::vector<double> &sisplines) {
  int i1 ;
  double dl, sl, x;
  int len = (k+1)*(k+2)/2;
  std::vector<double> Db(k);
  std::vector<double> work(len);

  sisplines.reserve(n*2*k*k);

  for(int i=0; i<(k-1); ++i) {
    for(int j=0; j<2*k*k; ++j) { 
      {sisplines.emplace_back(0.0);} 
    }
  }

  double xm1 = 0.0;
  double xs38a, xs38b;
  auto ia=0;
  auto ib=0;

  for(auto i=k-1; i<n; ++i){
    i1 = i + 1 ;
    dl = knots[i1] - knots[i];
    sl = knots[i1] + knots[i];

    for(int p=0; p<k; ++p){
      x = dl*0.5 * gl_x[p] + sl*0.5;    //x-transformation
      xs38a = (2*xm1+x)/3;
      xs38b = (xm1+2*x)/3;

      ia=i1-(knots[i]>xs38a);
      ib=i1-(knots[i]>xs38b);
      dbspvd_(&knots[0], k, 1, xs38a, ia, k, &Db[0], &work[0]);
      sisplines.insert(std::end(sisplines), std::begin(Db), std::end(Db));

      dbspvd_(&knots[0], k, 1, xs38b, ib, k, &Db[0], &work[0]);
      sisplines.insert(std::end(sisplines), std::begin(Db), std::end(Db));
      xm1=x;
    }
  }
  return 0;
}

int bsp::LobSplines(int n, int k, 
                std::vector<double> &gl_x,
                std::vector<double> &knots,
                std::vector<double> &lobsplines) {
  int i1 ;
  double dl, sl, x;
  int len = (k+1)*(k+2)/2;
  std::vector<double> Db(k);
  std::vector<double> work(len);

  lobsplines.reserve(n*2*k*k);

  for(int i=0; i<(k-1); ++i) {
    for(int j=0; j<2*k*k; ++j) { 
      {lobsplines.emplace_back(0.0);} 
    }
  }

  double xm1 = 0.0;
  double Lob4p = 0.4472135954999579392818347e0;
  double dlob=0.0;
  double slob=0.0;
  double xloba, xlobb;
  auto ia=0;
  auto ib=0;

  for(auto i=k-1; i<n; ++i){
    i1 = i + 1 ;
    dl = (knots[i1] - knots[i])*0.5;
    sl = (knots[i1] + knots[i])*0.5;

    for(int p=0; p<k; ++p){
      x = dl*gl_x[p] + sl;    //x-transformation
      dlob = (x-xm1)*0.5*Lob4p;
      slob = (x+xm1)*0.5;
      xloba = -dlob + slob;
      xlobb = dlob  + slob;

      ia=i1-(knots[i]>xloba);
      ib=i1-(knots[i]>xlobb);
      dbspvd_(&knots[0], k, 1, xloba, ia, k, &Db[0], &work[0]);
      lobsplines.insert(std::end(lobsplines), std::begin(Db), std::end(Db));

      dbspvd_(&knots[0], k, 1, xlobb, ib, k, &Db[0], &work[0]);
      lobsplines.insert(std::end(lobsplines), std::begin(Db), std::end(Db));
      xm1=x;
    }
  }
  return 0;
}

int bsp::GL2Splines(int n, int k, int k2,
                std::vector<double> &gl_outer,
                std::vector<double> &gl_inner,
                std::vector<double> &knots,
                std::vector<double> &glsplines) {
  int i1 ;
  double dl, sl, x;
  int len = (k+1)*(k+2)/2;
  std::vector<double> Db(k);
  std::vector<double> work(len);

  glsplines.reserve(n*k2*k*k);

  for(int i=0; i<(k-1); ++i) {
    for(int j=0; j<k2*k*k; ++j) { 
      {glsplines.emplace_back(0.0);} 
    }
  }

  double xm1 = 0.0;
  double gl2x;
  auto ia=0;

  for(auto i=k-1; i<n; ++i){
    i1 = i + 1 ;
    dl = knots[i1] - knots[i];
    sl = knots[i1] + knots[i];

    for(int p=0; p<k; ++p){
      x = dl*0.5 * gl_outer[p] + sl*0.5;    //x-transformation

      for (int j=0; j<k2; ++j) {
        gl2x=(x-xm1)*0.5*gl_inner[j] + (x+xm1)*0.5;

        ia=i1-(knots[i]>gl2x);
        dbspvd_(&knots[0], k, 1, gl2x, ia, k, &Db[0], &work[0]);
        glsplines.insert(std::end(glsplines), std::begin(Db), std::end(Db));
      }
      xm1=x;
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

  splinesp.reserve(n*k*k);

  for(int i=0; i<(k-1)*k*k; ++i) {splinesp.emplace_back(0.0);}

  for(auto i=k-1; i<n; ++i){
    dl = knots[i+1] - knots[i];
    sl = knots[i+1] + knots[i];
    i1 = i + 1 ;

    for(int p=0; p<k; ++p){
      x = dl*0.5 * gl_x[p] + sl*0.5;    //x-transformation
      
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
                  std::unique_ptr<ModelV> &Vptr) {
  int j_max, t_min, t_max, tm1;  
  double ovlp, bsum, dl, sl, x;

  for(int i=0; i<n; ++i) {
    j_max = std::min(i+k,n);
    for(int j=i; j<j_max; ++j) {
      t_min = std::max(k, j+2);
      t_max = std::min(i+k+1, n+2);
      ovlp=0.0;

      for(int t=t_min; t<=t_max; ++t) {
        tm1 = t-1;
        dl = (kkn[t] - kkn[tm1])*0.5;
        sl = (kkn[t] + kkn[tm1])*0.5;
        bsum = 0.0;
        for(int p=0; p<k; ++p) {
          x=dl*gl_x[p]+sl;
          bsum += gl_w[p]*spl[i+1-t+k+k*(p+(tm1)*k)]*spl[j+1-t+k+k*(p+(tm1)*k)]*Vptr->V(x);
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