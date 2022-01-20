#include "CIdiag.h"

int CalcCI(std::string pot, char gauge, int L_max) {
  std::string filename;
  std::unique_ptr<H5::H5File> efile=nullptr, file=nullptr;
  std::unique_ptr<H5::DataSet> edata=nullptr, v12data=nullptr, diag_set=nullptr;
  std::unique_ptr<H5::DataSet> L_set=nullptr;
  std::vector<double> ens, v12, eig, vecs;
  std::vector<int> ifail;
  std::vector<idx4> L_idx;
  int m, L_sz, v_sz, w_sz;
  std::ofstream dat_fl;

  for(auto L=0; L<L_max; ++L) {

    // read sum energies
    filename = pot + "2_" + std::to_string(L) + std::to_string(L+1) + gauge + ".h5";
    efile = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    edata = std::unique_ptr<H5::DataSet>(new H5::DataSet(efile->openDataSet("e_i")));
    L_sz = edata->getSpace().getSimpleExtentNpoints();
    ens.reserve(L_sz);
    edata->read(&ens[0], H5::PredType::NATIVE_DOUBLE);
    efile->close();

    v_sz = L_sz*(L_sz+1)/2;
    v12.reserve(v_sz);
    vecs.resize(L_sz*L_sz);

    // read V_12
    filename = pot + "V12_" + std::to_string(L) + ".h5";
    file = std::unique_ptr<H5::H5File>(
            new H5::H5File(filename, H5F_ACC_RDONLY));
    v12data=std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("V_12")));
    v12data->read(&v12[0], H5::PredType::NATIVE_DOUBLE);
    file->close();

    // get idx 
    filename = pot + std::to_string(L) + "idx.h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
    L_idx.reserve(L_sz);
    L_set->read(&L_idx[0], H5::PredType::NATIVE_INT32);
    file->close();

    dat_fl.open("dat/r12-"+std::to_string(L)+".dat");
    for(auto i=0; i<L_sz; ++i) {
      v12[(2*L_sz-i-1)*i/2+i] += ens[i];
      dat_fl
            // <<L_idx[i].n1<<" "<<L_idx[i].l1<<" "<<L_idx[i].n2<<" "
            // <<L_idx[i].l2<<" "<<ens[i]<<" "
            << std::abs(v12[(2*L_sz-i-1)*i/2+i]) << "\n";
      vecs[i*L_sz+i] = v12[(2*L_sz-i-1)*i/2+i];
      for(auto j=i+1; j<L_sz; ++j) {
        vecs[i*L_sz+j] = v12[(2*L_sz-i-1)*i/2+j];
        dat_fl
              // <<L_idx[i].n1<<" "<<L_idx[i].l1<<" "<<L_idx[i].n2<<" "
              // <<L_idx[i].l2<<" "
              // <<L_idx[j].n1<<" "<<L_idx[j].l1<<" "<<L_idx[j].n2<<" "
              // <<L_idx[j].l2<<" "
              << std::abs(v12[(2*L_sz-i-1)*i/2+j]) << "\n";
        // vecs[j*L_sz+i] = vecs[i*L_sz+j];
      }
    }
    dat_fl.close();

    eig.resize(L_sz);
    // ifail.reserve(L_sz);
    // w_sz = 10*L_sz;

    // std::cout << LAPACKE_dsyevx(LAPACK_ROW_MAJOR, 'V', 'I', 'U', L_sz, &v12[0], L_sz,
    //   0.0, 0.0, 1, L_sz, 0.0, &m, &eig[0], &vecs[0], L_sz, &ifail[0]) << "\n";

    std::cout << LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', L_sz, &vecs[0], L_sz, &eig[0]) << "\n";

    hsize_t d1[1] = {L_sz};
    hsize_t d2[2] = {L_sz,L_sz};

    filename = pot + "_diag" + std::to_string(L) + ".h5";
    file = std::unique_ptr<H5::H5File>(
            new H5::H5File(filename, H5F_ACC_TRUNC));
    diag_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->createDataSet(
                    "En", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, d1))));
    diag_set->write(&eig[0], H5::PredType::NATIVE_DOUBLE);
    // diag_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->createDataSet(
    //                 "EnVec", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, d2))));
    // diag_set->write(&vecs[0], H5::PredType::NATIVE_DOUBLE);
    file->close();
    vecs.clear();
  }

  filename = pot + "2_" + std::to_string(L_max-1) + std::to_string(L_max) + gauge + ".h5";
  efile = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  edata = std::unique_ptr<H5::DataSet>(new H5::DataSet(efile->openDataSet("e_f")));
  L_sz = edata->getSpace().getSimpleExtentNpoints();
  ens.reserve(L_sz);
  edata->read(&ens[0], H5::PredType::NATIVE_DOUBLE);
  efile->close();

  v_sz = L_sz*(L_sz+1)/2;
  v12.reserve(L_sz*L_sz);
  vecs.resize(L_sz*L_sz);

  // read V_12
  filename = pot + "V12_" + std::to_string(L_max) + ".h5";
  file = std::unique_ptr<H5::H5File>(
          new H5::H5File(filename, H5F_ACC_RDONLY));
  v12data=std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("V_12")));
  v12data->read(&v12[0], H5::PredType::NATIVE_DOUBLE);
  file->close();

  // get idx 
  filename = pot + std::to_string(L_max) + "idx.h5";
  file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
  L_idx.reserve(L_sz);
  L_set->read(&L_idx[0], H5::PredType::NATIVE_INT32);
  file->close();

  dat_fl.open("dat/r12-"+std::to_string(L_max)+".dat");
  for(auto i=0; i<L_sz; ++i) {
    v12[(2*L_sz-i-1)*i/2+i] += ens[i];
    dat_fl
          <<L_idx[i].n1<<" "<<L_idx[i].l1<<" "<<L_idx[i].n2<<" "
          <<L_idx[i].l2<<" "
          << std::abs(v12[(2*L_sz-i-1)*i/2+i]) << "\n";
    vecs[i*L_sz+i] = v12[(2*L_sz-i-1)*i/2+i];
    for(auto j=i+1; j<L_sz; ++j) {
      vecs[i*L_sz+j] = v12[(2*L_sz-i-1)*i/2+j];
      dat_fl
            <<L_idx[i].n1<<" "<<L_idx[i].l1<<" "<<L_idx[i].n2<<" "
            <<L_idx[i].l2<<" "
            <<L_idx[j].n1<<" "<<L_idx[j].l1<<" "<<L_idx[j].n2<<" "
            <<L_idx[j].l2<<" "
            << std::abs(v12[(2*L_sz-i-1)*i/2+j]) << "\n";
      // vecs[j*L_sz+i] = vecs[i*L_sz+j];
    }
  }
  dat_fl.close();
  
  eig.resize(L_sz);

  std::cout << LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', L_sz, &vecs[0], L_sz, &eig[0]) << "\n";

  hsize_t d1[1] = {L_sz};
  hsize_t d2[2] = {L_sz,L_sz};

  filename = pot + "_diag" + std::to_string(L_max) + ".h5";
  file = std::unique_ptr<H5::H5File>(
          new H5::H5File(filename, H5F_ACC_TRUNC));
  diag_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->createDataSet(
                  "En", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, d1))));
  diag_set->write(&eig[0], H5::PredType::NATIVE_DOUBLE);
  // diag_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->createDataSet(
  //                 "EnVec", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, d2))));
  // diag_set->write(&vecs[0], H5::PredType::NATIVE_DOUBLE);
  file->close();

  return 0;
}