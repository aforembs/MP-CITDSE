#include "cfg_in.h"

int cfg::ReadCfg(std::string dir, int L, int &sym, int &ncf, std::vector<cfg::line> &cfgs) {
  std::string filename = dir + "/cfg-" + std::to_string(L) + ".inp"; 
  std::ifstream cfgfile(filename);
  std::string line;
  std::vector<int> vals;
  int val;

  if(!std::filesystem::exists(filename)) {
    std::cout << "Input file: " << filename << " does not exist!\n";
    return -1;
  }

  do {
    std::getline(cfgfile, line);
    std::istringstream iss(line);
    vals.clear();
    std::copy(std::istream_iterator<int>(iss),
              std::istream_iterator<int>(),
              std::back_inserter(vals));
  } while(vals.size()==1);

  sym = 2*vals[1]+1;

  std::getline(cfgfile, line);
  std::istringstream iss(line);
  iss >> ncf;

  cfg::line cf_l;
  cfgs.reserve(ncf);

  for(auto i=0; i<ncf; ++i) {
    std::getline(cfgfile, line);
    std::istringstream iss(line);
    if(iss >> val) {
      cf_l.n1 = val-1; // -1 for C indexing
      iss >> val;
      cf_l.l1 = val;
      iss >> val;
      cf_l.l2 = val;
      iss >> val;
      cf_l.n2min = val-1;
      iss >> val;
      cf_l.n2max = val;
      if(abs(cf_l.l1-cf_l.l2)<=L && L<=(cf_l.l1+cf_l.l2) && ((L>>0)&1)==(((cf_l.l1+cf_l.l2)>>0)&1)) 
        cfgs.emplace_back(cf_l);
      else
        std::cout << "cfg-"<< L <<".inp, line: (" << i+1 << ") has an invalid configuration, skipping\n";
    } else {
      --i; // ignore empty lines
    }
  }
  return 0;
}

int cfg::GenL_idx(std::string pot, char gauge, int L_max, std::string dir) {
  int l1=0;
  int l2=0;
  std::string filename;
  std::string outfile_name;
  std::unique_ptr<H5::DataSet> e1=nullptr, ei=nullptr;
  std::unique_ptr<H5::H5File> file=nullptr;
  std::unique_ptr<H5::H5File> outfile=nullptr;
  std::vector<double> en12, en_srtd;
  std::vector<idx4> idx_data;
  std::vector<en_data> Li_dat, Lf_dat;

  int ncf, sym;
  std::vector<cfg::line> cfgs;
  cfg::ReadCfg(dir, 0, sym, ncf, cfgs);

  en_data en_d;
  idx4 idx_elm;

  H5::DataSpace e_space;
  hsize_t offset[1], count[1], stride[1], block[1], dimms[1];
  offset[0]=0;
  stride[0]=1;
  block[0] =1;
  H5::DataSpace memspace_l1, memspace_l;

  hsize_t t_sz=0;
  for(auto i=0; i<ncf; ++i) {
    t_sz += cfgs[i].n2max-cfgs[i].n2min;
  }

  auto max_n2l = *std::max_element(cfgs.begin(), cfgs.end(),
        [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
  auto max_Nsz = max_n2l.n2max;
  en12.reserve(max_Nsz*2);
  count[0] =max_Nsz;
  dimms[0] =max_Nsz;
  memspace_l.setExtentSimple(1, dimms, NULL);

  // auto max_line = *std::max_element(cfgs.begin(), cfgs.end(),
  //       [](cfg::line const &a, cfg::line const &b) { return a.n1 < b.n1; });

  // auto l1_sz = max_line.n1;
  // count1[0]=l1_sz;
  // dimms1[0]=l1_sz;
  // memspace_l1.setExtentSimple(1, dimms1, NULL);

  // Sort the energies in parallel using C++17 built in parallel sort
  // std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(), 
  //       [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });

  // std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(), 
  //       [](cfg::line const &a, cfg::line const &b) {return (a.l2 < b.l2) && (a.l1==b.l1);});

  // L initial sort & save
  Li_dat.reserve(t_sz);
  en_srtd.reserve(t_sz);
  idx_data.reserve(t_sz);
  hsize_t write_sz[] = {t_sz};
  hsize_t idx_sz[] = {t_sz*4};

  filename = pot +"0.h5";
  file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("En")));
  e_space = e1->getSpace();
  e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
  e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);

  e_space = e1->getSpace();
  e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
  e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);

  int last_l1=0, last_l2=0;
  for(const auto &line : cfgs) {
    l1=line.l1;
    l2=line.l2;

    if(l1!=last_l1) {
      filename = pot + std::to_string(l1)+".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("En")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);
    }

    if(l2!=last_l2) {
      filename = pot + std::to_string(l2)+".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("En")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);
    }

    double ent, en1=en12[line.n1];
    for(int n2=line.n2min; n2<line.n2max; ++n2) {
      ent = en1 + en12[max_Nsz+n2];
      en_d = {ent,line.n1,l1,n2,l2};
      Li_dat.emplace_back(en_d);
    }
    last_l1=l1;
    last_l2=l2;
  }

  // Sort the energies in parallel using C++17 built in parallel sort
  // std::sort(std::execution::par_unseq, Li_dat.begin(), Li_dat.end(), 
  //       [](en_data const &a, en_data const &b) { return a.en < b.en; });

  // Separate the sorted energies from the sorted indices
  for(auto i : Li_dat) {
    en_srtd.emplace_back(i.en);
    idx_elm = {i.n1,i.l1,i.n2,i.l2};
    idx_data.emplace_back(idx_elm);
  }

  // write the energies to a file
  outfile_name = pot + "2_01" + gauge + ".h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("e_i", H5::PredType::NATIVE_DOUBLE, 
        H5::DataSpace(1, write_sz))));
  ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);
  outfile->close();

  en_srtd.clear();

  // write indices to file
  outfile_name = pot + "0idx.h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("idx", H5::PredType::NATIVE_INT32, 
        H5::DataSpace(1, idx_sz))));
  ei->write(&idx_data[0], H5::PredType::NATIVE_INT32);
  outfile->close();

  idx_data.clear();
  Li_dat.clear();
  cfgs.clear();

  // L final sort & save
  cfg::ReadCfg(dir, 1, sym, ncf, cfgs);

  t_sz=0;
  for(auto i=0; i<ncf; ++i) {
    t_sz += cfgs[i].n2max-cfgs[i].n2min;
  }

  max_n2l = *std::max_element(cfgs.begin(), cfgs.end(),
        [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
  max_Nsz = max_n2l.n2max;
  en12.reserve(max_Nsz*2);
  count[0] =max_Nsz;
  dimms[0] =max_Nsz;
  memspace_l.setExtentSimple(1, dimms, NULL);

  // max_line = *std::max_element(cfgs.begin(), cfgs.end(),
  //       [](cfg::line const &a, cfg::line const &b) { return a.n1 < b.n1; });

  // l1_sz = max_line.n1;
  // count1[0]=l1_sz;
  // dimms1[0]=l1_sz;
  // memspace_l1.setExtentSimple(1, dimms1, NULL);

  // Sort the energies in parallel using C++17 built in parallel sort
  // std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(), 
  //       [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });

  // std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(), 
  //       [](cfg::line const &a, cfg::line const &b) {return (a.l2 < b.l2) && (a.l1==b.l1);});

  Lf_dat.reserve(t_sz);
  en_srtd.reserve(t_sz);
  idx_data.reserve(t_sz);
  write_sz[0]=t_sz;
  idx_sz[0]=t_sz*4;

  filename = pot +"0.h5";
  file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("En")));
  e_space = e1->getSpace();
  e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
  e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);

  e_space = e1->getSpace();
  e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
  e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);
  file->close();

  last_l1=0; last_l2=0;
  for(const auto &line : cfgs) {
    l1=line.l1;
    l2=line.l2;

    if(l1!=last_l1) {
      filename = pot + std::to_string(l1)+ ".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("En")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);
      file->close();
    }

    if(l2!=last_l2) {
      filename = pot + std::to_string(l2)+".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("En")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);
      file->close();
    }

    double ent, en1=en12[line.n1];
    for(int n2=line.n2min; n2<line.n2max; ++n2) {
      ent = en1 + en12[max_Nsz+n2];
      en_d = {ent,line.n1,l1,n2,l2};
      Lf_dat.emplace_back(en_d);
    }
    last_l1=l1;
    last_l2=l2;
  }

  // std::sort(std::execution::par_unseq, Lf_dat.begin(), Lf_dat.end(), 
  //       [](en_data const &a, en_data const &b) { return a.en < b.en; });

  for(auto i : Lf_dat) {
    en_srtd.emplace_back(i.en);
    idx_elm = {i.n1,i.l1,i.n2,i.l2};
    idx_data.emplace_back(idx_elm);
  }

  outfile_name = pot + "2_01" + gauge + ".h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_RDWR));
  ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("e_f", H5::PredType::NATIVE_DOUBLE, 
        H5::DataSpace(1, write_sz))));
  ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);
  outfile->close();

  // write indices to file
  outfile_name = pot + "1idx.h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("idx", H5::PredType::NATIVE_INT32, 
        H5::DataSpace(1, idx_sz))));
  ei->write(&idx_data[0], H5::PredType::NATIVE_INT32);
  outfile->close();

  Lf_dat.clear();

  // Loop
  for(int L_itr=2; L_itr<=L_max; ++L_itr) {
    outfile_name = pot + "2_" + std::to_string(L_itr-1) + std::to_string(L_itr) + gauge + ".h5";
    outfile = std::unique_ptr<H5::H5File>(
              new H5::H5File(outfile_name, H5F_ACC_TRUNC));
    ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
          outfile->createDataSet("e_i", H5::PredType::NATIVE_DOUBLE, 
          H5::DataSpace(1, write_sz))));
    ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);

    en_srtd.clear();
    idx_data.clear();
    Lf_dat.clear();
    cfgs.clear();

    cfg::ReadCfg(dir, L_itr, sym, ncf, cfgs);

    t_sz=0;
    for(auto i=0; i<ncf; ++i) {
      t_sz += cfgs[i].n2max-cfgs[i].n2min;
    }

    max_n2l = *std::max_element(cfgs.begin(), cfgs.end(),
          [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
    max_Nsz = max_n2l.n2max;
    en12.reserve(max_Nsz*2);
    count[0] =max_Nsz;
    dimms[0] =max_Nsz;
    memspace_l.setExtentSimple(1, dimms, NULL);

    // max_line = *std::max_element(cfgs.begin(), cfgs.end(),
    //       [](cfg::line const &a, cfg::line const &b) { return a.n1 < b.n1; });

    // l1_sz = max_line.n1;
    // count1[0]=l1_sz;
    // dimms1[0]=l1_sz;
    // memspace_l1.setExtentSimple(1, dimms1, NULL);

    // Sort the energies in parallel using C++17 built in parallel sort
    // std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(), 
    //       [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });

    // std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(), 
    //       [](cfg::line const &a, cfg::line const &b) {return (a.l2 < b.l2) && (a.l1==b.l1);});

    Lf_dat.reserve(t_sz);
    en_srtd.reserve(t_sz);
    idx_data.reserve(t_sz);
    write_sz[0]=t_sz;
    idx_sz[0]=t_sz*4;

    filename = pot +"0.h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("En")));
    e_space = e1->getSpace();
    e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);

    e_space = e1->getSpace();
    e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);
    file->close();

    last_l1=0; last_l2=0;
    for(const auto &line : cfgs) {
      l1=line.l1;
      l2=line.l2;

      if(l1!=last_l1) {
        filename = pot + std::to_string(l1)+".h5";
        file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
        e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("En")));
        e_space = e1->getSpace();
        e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);
        file->close();
      }

      if(l2!=last_l2) {
        filename = pot + std::to_string(l2)+".h5";
        file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
        e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("En")));
        e_space = e1->getSpace();
        e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l, e_space);
        file->close();
      }

      double ent, en1=en12[line.n1];
      for(int n2=line.n2min; n2<line.n2max; ++n2) {
        ent = en1 + en12[max_Nsz+n2];
        en_d = {ent,line.n1,l1,n2,l2};
        Lf_dat.emplace_back(en_d);
      }
      last_l1=l1;
      last_l2=l2;
    }

    // std::sort(std::execution::par_unseq, Lf_dat.begin(), Lf_dat.end(), 
    //     [](en_data const &a, en_data const &b) { return a.en < b.en; });

    for(auto i : Lf_dat) {
      en_srtd.emplace_back(i.en);
      idx_elm = {i.n1,i.l1,i.n2,i.l2};
      idx_data.emplace_back(idx_elm);
    }

    ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
          outfile->createDataSet("e_f", H5::PredType::NATIVE_DOUBLE, 
          H5::DataSpace(1, write_sz))));
    ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);
    outfile->close();

    // write indices to file
    outfile_name = pot + std::to_string(L_itr) + "idx.h5";
    outfile = std::unique_ptr<H5::H5File>(
              new H5::H5File(outfile_name, H5F_ACC_TRUNC));
    ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
          outfile->createDataSet("idx", H5::PredType::NATIVE_INT32, 
          H5::DataSpace(1, idx_sz))));
    ei->write(&idx_data[0], H5::PredType::NATIVE_INT32);
    outfile->close();
  }

  return 0;
}