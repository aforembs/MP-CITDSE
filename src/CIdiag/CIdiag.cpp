int CalcCI(std::string pot, int L_max) {
  std::string filename;
  std::unique_ptr<H5::H5File> efile=nullptr, file=nullptr;
  std::unique_ptr<H5::DataSet> edata=nullptr;
  std::vector<double> ens, v12, eig, vecs, work, iwork;
  auto m=0, w_sz=0;

  for(auto L=0; L<L_max; ++L) {

    // read sum energies
    filename = pot + "2_" + std::to_string(L) + std::to_string(L+1) + gauge + ".h5";
    efile = std::unique_ptr<H5::H5File>(
            new H5::H5File(filename, H5F_ACC_TRUNC));
    edata = std::unique_ptr<H5::DataSet>(new H5::DataSet(efile->openDataSet("e_i")));
    L_sz = edata->getSpace().getSimpleExtentNpoints();
    ens.reserve(L_sz);
    edata->read(&ens[0], H5::PredType::NATIVE_DOUBLE);

    v_sz = L_sz*(L_sz+1)/2;
    v12.reserve(v_sz);

    // read V_12
    filename = pot + "V12_" + std::to_string(L) + ".h5";
    file = std::unique_ptr<H5::H5File>(
            new H5::H5File(filename, H5F_ACC_TRUNC));
    v12data=std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("V_12")));
    v12data->read(&v12[0], H5::PredType::NATIVE_DOUBLE);

    for(auto i=0; i<L_sz; ++i) {
      v12[(2*L_sz-i-1)*i/2+i] += ens[i];
    }

    eig.reserve(L_sz);
    vecs.reserve(L_sz*L_sz);
    w_sz = 8*L_sz;
    work.reserve(w_sz);
    iwork.reserve(5*L_sz);

    LAPACKE_dsyevx(LAPACK_ROW_MAJOR, 'V', 'I', 'U', L_sz, &v12[0], L_sz,
      0.0, 0.0, 1, L_sz, 0.0, &m, &eig[0], &vecs[0], L_sz, &work[0], 
      w_sz, &iwork[0], &ifail);

    filename = pot + "_diag" + std::to_string(L) + ".h5";
    file = std::unique_ptr<H5::H5File>(
            new H5::H5File(filename, H5F_ACC_TRUNC));
    diag_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->createDataSet(
                    "En", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, L_sz))));
    diag_set->write(&eig[0], H5::PredType::NATIVE_DOUBLE);
    diag_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->createDataSet(
                    "EnVec", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, {L_sz,L_sz}))));
    diag_set->write(&vecs[0], H5::PredType::NATIVE_DOUBLE);
  }

  edata = std::unique_ptr<H5::DataSet>(new H5::DataSet(efile->openDataSet("e_f")));
  L_sz = edata->getSpace().getSimpleExtentNpoints();
  ens.reserve(L_sz);
  edata->read(&ens[0], H5::PredType::NATIVE_DOUBLE);

  v_sz = L_sz*(L_sz+1)/2;
  v12.reserve(v_sz);

  // read V_12
  filename = pot + "V12_" + std::to_string(L_max) + ".h5";
  file = std::unique_ptr<H5::H5File>(
          new H5::H5File(filename, H5F_ACC_TRUNC));
  v12data=std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("V_12")));
  v12data->read(&v12[0], H5::PredType::NATIVE_DOUBLE);

  for(auto i=0; i<L_sz; ++i) {
    v12[(2*L_sz-i-1)*i/2+i] += ens[i];
  }

  eig.reserve(L_sz);
  vecs.reserve(L_sz*L_sz);
  w_sz = 8*L_sz;
  work.reserve(w_sz);
  iwork.reserve(5*L_sz);

  LAPACKE_dsyevx(LAPACK_ROW_MAJOR, 'V', 'I', 'U', L_sz, &v12[0], L_sz,
    0.0, 0.0, 1, L_sz, 0.0, &m, &eig[0], &vecs[0], L_sz, &work[0], 
    w_sz, &iwork[0], &ifail);

  filename = pot + "_diag" + std::to_string(L_max) + ".h5";
  file = std::unique_ptr<H5::H5File>(
          new H5::H5File(filename, H5F_ACC_TRUNC));
  diag_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->createDataSet(
                  "En", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, L_sz))));
  diag_set->write(&eig[0], H5::PredType::NATIVE_DOUBLE);
  diag_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->createDataSet(
                  "EnVec", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, {L_sz,L_sz}))));
  diag_set->write(&vecs[0], H5::PredType::NATIVE_DOUBLE);

  return 0;
}