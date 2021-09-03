#include "V_12.h"

int V12(uint L_max, std::vector<uint> &N_sz) {
  uint min_dir=0, min_exc=0;
  double Y_norm=0.0;
  double F_dir=1.0; // dummy slater integrals
  double F_exc=F_dir;
  double v_mat; // standin for V_12 matrix
  double sum_k=0.0;
  idx4 e12, e12p;

  int n=0;   // no. of points
  int k=0;   // max B-spline order
  int nkn=0; // no. of knots
  int nCf=0; // no. of coefficients
  int nSt=0; // no. of states
  std::vector<double> kkn;
  std::vector<double> Cf;
  std::vector<double*> C(L_max+1);

  uint L_sz=0, v_sz=0;
  std::string filename;
  std::string outfile_name;
  H5::H5File *file=nullptr;
  H5::H5File *outfile=nullptr;
  std::vector<idx4> idx_data;
  H5::DataSet *L_set=nullptr;
  H5::DataSet *V_set=nullptr;
  H5::DataSpace cspace;
  hsize_t offset[2], count[2], stride[2], block[2];
  hsize_t dimms[2];
  offset[0]=0; offset[1]=0;
  count[0] =0;
  stride[0]=1; stride[1]=1;
  block[0] =1; block[1]=1;
  dimms[0] =0;
  H5::DataSpace memspace;

  int tot_states = 0;
  for(auto &n : N_sz) {
    tot_states += n;
  }

  filename = cpot + std::to_string(0) + ".h5";
  file = new H5::H5File(filename, H5F_ACC_RDONLY);
  file->openAttribute("N").read(H5::PredType::NATIVE_INT32, &n);
  file->openAttribute("K").read(H5::PredType::NATIVE_INT32, &k);
  rset = new H5::DataSet(file->openDataSet("Knots"));
  nkn = rset->getSpace().getSimpleExtentNpoints();

  kkn.reserve(nkn);
  rset->read(&kkn[0], H5::PredType::NATIVE_DOUBLE); //read knots
  delete rset;

  rset = new H5::DataSet(file->openDataSet("En"));
  nSt = rset->getSpace().getSimpleExtentNpoints();
  delete rset;
  delete file;

  Cf.reserve(tot_states*n);
  C[0]=&Cf[1];
  Cf[0]=0.0; Cf[nSt+1] =0.0;
  int nst_prev = 0;
  for(int i=1; i<=L_max; ++i) {
    nst_prev += N_max[i-1];
    C[i] = &Cf[nst_prev*n+1];
    Cf[nst_prev*n]=0.0;
    Cf[i*n-1]=0.0;
  }

  // read coefficients for all l
  for(int l=0; l<=L_max; ++l) {
    count[0] = N_max[l]; count[1] = nSt;
    dimms[0] = count[0]; dimms[1] = count[1];
    memspace.setExtentSimple(2, dimms, NULL);

    filename = cpot + std::to_string(l) + ".h5";
    file = new H5::H5File(filename, H5F_ACC_RDONLY);
    rset = new H5::DataSet(file->openDataSet("Coeff"));
    cspace = rset->getSpace();
    cspace.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    rset->read(C[l], H5::PredType::NATIVE_DOUBLE, memspace, cspace);
    delete rset;
    delete file;
  }

  std::vector<double> gl_x(k);
  for(int i=1; i<=k; ++i) {
    gl_x[k-i] = fastgl::GLPair(k, i).x(); // generate GL nodes over B-splines support
  }

  // generate B-splines
  int len = (k+1)*(k+2)/2;
  std::vector<double> Db(k);
  std::vector<double> work(len);
  std::vector<double> Bsplines;

  Bsplines.reserve(nSt*k*k);

  int i1 ;
  double dl, sl, x;
             
  for(auto i=k-1; i<n; ++i){
    dl = kkn[i+1] - kkn[i];
    sl = kkn[i+1] + kkn[i];

    for(int p=0; p<k; ++p){
      x = dl*0.5 * gl_x[p] + sl*0.5;    //x-transformation
      i1 = i + 1 ;
      dbspvd_(&kkn[0], k, 1, x, i1, k, &Db[0], &work[0]);

      Bsplines.insert(std::end(Bsplines), std::begin(Db), std::end(Db));
    }
  }

  filename = pot + std::to_string(Lf_i) + "idx.h5";
  file = new H5::H5File(filename, H5F_ACC_RDONLY);
  L_set = new H5::DataSet(file->openDataSet("idx"));
  L_sz = L_set->getSpace().getSimpleExtentNpoints()/4;
  std::vector<idx4> L_idx(Lidx_sz);
  delete L_set;
  delete file;

  std::vector<double> v_mat;

  for(uint L=0; L<=L_max; ++L) {
    // Read indices n1l1;n2l2
    filename = pot + std::to_string(L) + "idx.h5";
    file = new H5::H5File(filename, H5F_ACC_RDONLY);
    L_set = new H5::DataSet(file->openDataSet("idx"));
    L_sz = L_set->getSpace().getSimpleExtentNpoints()/4;
    L_set->read(&L_idx[0], H5::PredType::NATIVE_UINT32);
    delete L_set;
    delete file;

    v_sz = L_sz*(L_sz+1)/2;
    v_mat.reserve(v_sz);

    for(uint NL2=0; NL2<L_sz; ++NL2) {
      //set n1'l1';n2'l2'
      e12p = L_idx[NL2];

      for(uint NL1=NL2; NL1<L_sz; ++NL1){
        //set n1l1;n2l2
        e12 = L_idx[NL1];

        Y_norm = sqrt((2*e12.l1+1)*(2*e12p.l1+1)*(2*e12.l2+1)*(2*e12p.l2+1));
        sum_k=0.0;
        for(uint k=0; k<=l1e_max; ++k) {
          min_dir = ((L+e12.l2+e12p.l1) >> 0) & 1;
          if(min_dir ==(((L+e12.l1+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l1)<=k) && (k<=e12.l1+e12p.l1)) &&
             ((abs(e12.l2-e12p.l2)<=k) && (k<=e12.l2+e12p.l2))) {
            sum_k += pow(-1,min_dir)*F_dir*wigner_3j0(e12.l1,k,e12p.l1)
                  *wigner_3j0(e12.l2,k,e12p.l2)*wigner_6j(e12p.l1,e12.l2,L,e12.l2,e12.l1,k);
          }
          min_exc = ((L+e12.l1+e12p.l1) >> 0) & 1;
          if(min_exc ==(((L+e12.l2+e12p.l2) >> 0) & 1) &&
             ((abs(e12.l1-e12p.l2)<=k) && (k<=e12.l1+e12p.l2)) &&
             ((abs(e12.l2-e12p.l1)<=k) && (k<=e12.l2+e12p.l1))) {
            sum_k += pow(-1,min_exc)*F_exc*wigner_3j0(e12.l1,k,e12p.l2)
                  *wigner_3j0(e12.l2,k,e12p.l1)*wigner_6j(e12p.l1,e12p.l2,L,e12.l1,e12.l2,k);
          }
        }
        // write symmetric V_12 as upper triangular
        v_mat[(2*L_sz-NL2-1)*NL2/2 + NL1] = pow(-1,(l1+l2))*Y_norm*sum_k;
      }
    }
    // save upper triangular V_12
    outfile_name = pot + "V12" + std::to_string(L) + ".h5";
    outfile = new H5::H5File(outfile_name, H5F_ACC_TRUNC);
    V_set = new H5::DataSet(outfile->createDataSet("V_12", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, v_sz)));
    V_set->write(&v_mat[0], H5::PredType::NATIVE_DOUBLE);
    delete V_set;
    delete outfile;

    v_mat.clear();
  }

  return 0;
}