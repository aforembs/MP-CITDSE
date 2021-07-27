#include "V_12.h"

int V12(uint L_max) {
  uint min_dir=0, min_exc=0;
  double Y_norm=0.0;
  double F_dir=1.0; // dummy slater integrals
  double F_exc=F_dir;
  double v_mat; // standin for V_12 matrix
  double sum_k=0.0;
  idx4 e12, e12p;

  uint L_sz=0, v_sz=0;
  std::string filename;
  std::string outfile_name;
  H5::H5File *file=nullptr;
  H5::H5File *outfile=nullptr;
  std::vector<idx4> idx_data;
  H5::DataSet *L_set=nullptr;
  H5::DataSet *V_set=nullptr;

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