#include "dmx2e.h"

en_L DMX2e::make_enL(int L, int l1e_max) {
  int sz=0;
  int st=0; 
  int l1=0;
  int l2=0;
  std::vector<l_ab> lp;
  l_ab l_pair;
  
  if(L==0) { // For L=0, l1=l2
    sz = l1e_max+1;
    lp.reserve(sz);
    for(int i=0; i<sz; ++i) {
      l_pair = {i,i};
      lp.emplace_back(l_pair);
    }
  } else if(L==1) { // For L=1, l2=l1+1
    sz = l1e_max;
    lp.reserve(sz);
    for(int i=0; i<sz; ++i) {
      l_pair = {i, i+1};
      lp.emplace_back(l_pair);
    }
  } else if((L^1)==(L+1) && L!=0) { // L>=2, L even
    sz = l1e_max*2;
    lp.reserve(sz);
    l_pair = {0,L};
    lp.emplace_back(l_pair);
    for(l1=1; l1<=l1e_max; ++l1) {
      if (l1*2>=L) { st=l1; }
      else { st=L-l1; }
      for(l2=st; l2<=l1e_max; l2+=2) {
        l_pair = {l1,l2};
        lp.emplace_back(l_pair);
      }
    }
  } else if(L & 1 && L!=1) { // L>=3, L odd 
    sz = l1e_max*2;
    lp.reserve(sz);
    l_pair = {0,L};
    lp.emplace_back(l_pair);
    for(l1=1; l1<=l1e_max; ++l1) {
      if (l1*2+1>=L) { st=l1+1; }
      else { st=L-l1; }
      for(l2=st; l2<=l1e_max; l2+=2) {
        l_pair = {l1,l2};
        lp.emplace_back(l_pair);
      }
    }
  }

  auto en_vec = std::vector<en_data>();
  en_L Len = {L,lp,en_vec};
  return Len;
}

int DMX2e::sort_L(int L_max, std::vector<int> &N_sz) {
  int l1=0;
  int l2=0;
  std::string filename;
  std::string outfile_name;
  std::unique_ptr<H5::DataSet> e1, ei;
  std::unique_ptr<H5::H5File> file;
  std::unique_ptr<H5::H5File> outfile;
  std::vector<double> en12, en_srtd;
  std::vector<idx4> idx_data;
  int max_Nsz = *std::max_element(std::begin(N_sz), std::end(N_sz));
  en12.reserve(max_Nsz*2);
  en_data en_d;
  idx4 idx_elm;

  en_L Li = make_enL(0, L_max);
  en_L Lf = make_enL(1, L_max);
  int num_l_pairs = Li.l_pair.size();

  H5::DataSpace e_space;
  hsize_t offset[1], count[1], stride[1], block[1];
  hsize_t dimms[1];
  offset[0]=0;
  count[0] =N_sz[0];
  stride[0]=1;
  block[0] =1;
  dimms[0] =N_sz[0];
  H5::DataSpace memspace;

  // L initial sort & save
  hsize_t t_sz=0;
  for(auto i=0; i<num_l_pairs; ++i) {
    l1=Li.l_pair[i].l1;
    l2=Li.l_pair[i].l2;
    if (l1!=l2) {
      t_sz+=N_sz[l1]*N_sz[l2];
    } else {
      t_sz+=N_sz[l1]*(N_sz[l1]+1)/2;
    }
  }
  Li.en_dat.reserve(t_sz);
  en_srtd.reserve(t_sz);
  idx_data.reserve(t_sz);
  hsize_t write_sz[] = {t_sz};
  hsize_t idx_sz[] = {t_sz*4};

  for(auto i=0; i<num_l_pairs; ++i) {
    l1=Li.l_pair[i].l1;
    l2=Li.l_pair[i].l2;

    if(l1==l2) { // Calculate energies if l1=l2
      count[0]=N_sz[l1];
      dimms[0]=N_sz[l1];
      memspace.setExtentSimple(1, dimms, NULL);
      filename = pot + std::to_string(l1)+std::to_string(l1+1) + gauge + ".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace, e_space);

      double ent;
      for(int n1=0; n1<N_sz[l1]; ++n1) {
        for(int n2=n1; n2<N_sz[l2]; ++n2) {
          ent = en12[n1] + en12[n2];
          en_d = {ent,n1,l1,n2,l2};
          Li.en_dat.emplace_back(en_d);
        }
      }

    } else { // Calculate energies if l1!=l2
      count[0]=N_sz[l1];
      dimms[0]=N_sz[l1];
      memspace.setExtentSimple(1, dimms, NULL);
      filename = pot + std::to_string(l1)+std::to_string(l1+1) + gauge + ".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace, e_space);

      count[0]=N_sz[l2];
      dimms[0]=N_sz[l2];
      memspace.setExtentSimple(1, dimms, NULL);
      filename = pot + std::to_string(l2)+std::to_string(l2+1) + gauge + ".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace, e_space);

      double ent;
      for(int n1=0; n1<N_sz[l1]; ++n1) {
        for(int n2=0; n2<N_sz[l2]; ++n2) {
          ent = en12[n1] + en12[max_Nsz+n2];
          en_d = {ent,n1,l1,n2,l2};
          Li.en_dat.emplace_back(en_d);
        }
      }
    }
  }

  // Sort the energies in parallel using C++17 built in parallel sort
  std::sort(std::execution::par_unseq, Li.en_dat.begin(), Li.en_dat.end(), 
        [](en_data const &a, en_data const &b) { return a.en < b.en; });

  // Separate the sorted energies from the sorted indices
  for(auto i : Li.en_dat) {
    en_srtd.emplace_back(i.en);
    idx_elm = {i.n1,i.l1,i.n2,i.l2};
    idx_data.emplace_back(idx_elm);
  }

  // write the energies to a file
  outfile_name = pot + "2_" + std::to_string(Li.L) + std::to_string(Li.L+1) + gauge + ".h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("e_i", H5::PredType::NATIVE_DOUBLE, 
        H5::DataSpace(1, write_sz))));
  ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);

  en_srtd.clear();

  // write indices to file
  outfile_name = pot + std::to_string(Li.L) + "idx.h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("idx", H5::PredType::NATIVE_INT32, 
        H5::DataSpace(1, idx_sz))));
  ei->write(&idx_data[0], H5::PredType::NATIVE_INT32);

  idx_data.clear();

  // L final sort & save
  t_sz=0;
  num_l_pairs = Lf.l_pair.size();
  for(auto i=0; i<num_l_pairs; ++i) {
    l1=Lf.l_pair[i].l1;
    l2=Lf.l_pair[i].l2;
    if (l1!=l2) {
      t_sz+=N_sz[l1]*N_sz[l2];
    } else {
      t_sz+=N_sz[l1]*(N_sz[l1]+1)/2;
    }
  }
  Lf.en_dat.reserve(t_sz);
  en_srtd.reserve(t_sz);
  idx_data.reserve(t_sz);
  write_sz[0]=t_sz;
  idx_sz[0]=t_sz*4;

  for(auto i=0; i<num_l_pairs; ++i) {
    l1=Lf.l_pair[i].l1;
    l2=Lf.l_pair[i].l2;

    if(l1==l2) { // Calculate energies if l1=l2
      count[0]=N_sz[l1];
      dimms[0]=N_sz[l1];
      memspace.setExtentSimple(1, dimms, NULL);
      filename = pot + std::to_string(l1)+std::to_string(l1+1) + gauge + ".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace, e_space);

      double ent;
      for(int n1=0; n1<N_sz[l1]; ++n1) {
        for(int n2=n1; n2<N_sz[l2]; ++n2) {
          ent = en12[n1] + en12[n2];
          en_d = {ent,n1,l1,n2,l2};
          Lf.en_dat.emplace_back(en_d);
        }
      }

    } else { // Calculate energies if l1!=l2
      count[0]=N_sz[l1];
      dimms[0]=N_sz[l1];
      memspace.setExtentSimple(1, dimms, NULL);
      filename = pot + std::to_string(l1)+std::to_string(l1+1) + gauge + ".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace, e_space);

      count[0]=N_sz[l2];
      dimms[0]=N_sz[l2];
      memspace.setExtentSimple(1, dimms, NULL);
      filename = pot + std::to_string(l2)+std::to_string(l2+1) + gauge + ".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
      e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace, e_space);

      double ent;
      for(int n1=0; n1<N_sz[l1]; ++n1) {
        for(int n2=0; n2<N_sz[l2]; ++n2) {
          ent = en12[n1] + en12[max_Nsz+n2];
          en_d = {ent,n1,l1,n2,l2};
          Lf.en_dat.emplace_back(en_d);
        }
      }
    }
  }

  std::sort(std::execution::par_unseq, Lf.en_dat.begin(), Lf.en_dat.end(), 
        [](en_data const &a, en_data const &b) { return a.en < b.en; });

  for(auto i : Lf.en_dat) {
    en_srtd.emplace_back(i.en);
    idx_elm = {i.n1,i.l1,i.n2,i.l2};
    idx_data.emplace_back(idx_elm);
  }

  outfile_name = pot + "2_" + std::to_string(Li.L) + std::to_string(Li.L+1) + gauge + ".h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("e_f", H5::PredType::NATIVE_DOUBLE, 
        H5::DataSpace(1, write_sz))));
  ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);

  // write indices to file
  outfile_name = pot + std::to_string(Lf.L) + "idx.h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("idx", H5::PredType::NATIVE_INT32, 
        H5::DataSpace(1, idx_sz))));
  ei->write(&idx_data[0], H5::PredType::NATIVE_INT32);

  // Loop
  Li = Lf;
  Lf.en_dat.clear();
  for(int L_itr=1; L_itr<L_max; ++L_itr) {
    outfile_name = pot + "2_" + std::to_string(Li.L) + std::to_string(Li.L+1) + gauge + ".h5";
    outfile = std::unique_ptr<H5::H5File>(
              new H5::H5File(outfile_name, H5F_ACC_TRUNC));
    ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
          outfile->createDataSet("e_i", H5::PredType::NATIVE_DOUBLE, 
          H5::DataSpace(1, write_sz))));
    ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);

    en_srtd.clear();
    idx_data.clear();

    Lf = make_enL(L_itr+1, L_max);

    t_sz=0;
    num_l_pairs = Lf.l_pair.size();
    for(auto i=0; i<num_l_pairs; ++i) {
      l1=Lf.l_pair[i].l1;
      l2=Lf.l_pair[i].l2;
      if (l1!=l2) {
        t_sz+=N_sz[l1]*N_sz[l2];
      } else {
        t_sz+=N_sz[l1]*(N_sz[l1]+1)/2;
      }
    }
    Lf.en_dat.reserve(t_sz);
    en_srtd.reserve(t_sz);
    idx_data.reserve(t_sz);
    write_sz[0]=t_sz;
    idx_sz[0]=t_sz*4;

    for(auto i=0; i<num_l_pairs; ++i) {
      l1=Lf.l_pair[i].l1;
      l2=Lf.l_pair[i].l2;

      if(l1==l2) {
        count[0]=N_sz[l1];
        dimms[0]=N_sz[l1];
        memspace.setExtentSimple(1, dimms, NULL);
        filename = pot + std::to_string(l1)+std::to_string(l1+1) + gauge + ".h5";
        file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
        e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
        e_space = e1->getSpace();
        e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace, e_space);

        double ent;
        for(int n1=0; n1<N_sz[l1]; ++n1) {
          for(int n2=n1; n2<N_sz[l2]; ++n2) {
            ent = en12[n1] + en12[n2];
            en_d = {ent,n1,l1,n2,l2};
            Lf.en_dat.emplace_back(en_d);
          }
        }

      } else {
        count[0]=N_sz[l1];
        dimms[0]=N_sz[l1];
        memspace.setExtentSimple(1, dimms, NULL);
        filename = pot + std::to_string(l1)+std::to_string(l1+1) + gauge + ".h5";
        file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
        e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
        e_space = e1->getSpace();
        e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace, e_space);

        count[0]=N_sz[l2];
        dimms[0]=N_sz[l2];
        memspace.setExtentSimple(1, dimms, NULL);
        filename = pot + std::to_string(l2)+std::to_string(l2+1) + gauge + ".h5";
        file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
        e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
        e_space = e1->getSpace();
        e_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
        e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace, e_space);

        double ent;
        for(int n1=0; n1<N_sz[l1]; ++n1) {
          for(int n2=0; n2<N_sz[l2]; ++n2) {
            ent = en12[n1] + en12[max_Nsz+n2];
            en_d = {ent,n1,l1,n2,l2};
            Lf.en_dat.emplace_back(en_d);
         }
        }
      }
    }

    std::sort(std::execution::par_unseq, Lf.en_dat.begin(), Lf.en_dat.end(), 
        [](en_data const &a, en_data const &b) { return a.en < b.en; });

    for(auto i : Lf.en_dat) {
      en_srtd.emplace_back(i.en);
      idx_elm = {i.n1,i.l1,i.n2,i.l2};
      idx_data.emplace_back(idx_elm);
    }

    ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
          outfile->createDataSet("e_f", H5::PredType::NATIVE_DOUBLE, 
          H5::DataSpace(1, write_sz))));
    ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);

    // write indices to file
    outfile_name = pot + std::to_string(Lf.L) + "idx.h5";
    outfile = std::unique_ptr<H5::H5File>(
              new H5::H5File(outfile_name, H5F_ACC_TRUNC));
    ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
          outfile->createDataSet("idx", H5::PredType::NATIVE_INT32, 
          H5::DataSpace(1, idx_sz))));
    ei->write(&idx_data[0], H5::PredType::NATIVE_INT32);

    Li=Lf;
  }

  return 0;
}

int DMX2e::sort_L(int L_max, std::string dir) {
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
  hsize_t offset[1], count1[1], count2[1], stride[1], block[1];
  hsize_t dimms2[1], dimms1[1];
  offset[0]=0;
  stride[0]=1;
  block[0] =1;
  H5::DataSpace memspace_l1, memspace_l2;

  hsize_t t_sz=0;
  for(auto i=0; i<ncf; ++i) {
    t_sz += cfgs[i].n2max-cfgs[i].n2min;
  }

  auto max_n2l = *std::max_element(cfgs.begin(), cfgs.end(),
        [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
  auto max_Nsz = max_n2l.n2max;
  en12.reserve(max_Nsz*2);
  count2[0] =max_Nsz;
  dimms2[0] =max_Nsz;
  memspace_l2.setExtentSimple(1, dimms2, NULL);

  auto max_line = *std::max_element(cfgs.begin(), cfgs.end(),
        [](cfg::line const &a, cfg::line const &b) { return a.n1 < b.n1; });

  auto l1_sz = max_line.n1;
  count1[0]=l1_sz;
  dimms1[0]=l1_sz;
  memspace_l1.setExtentSimple(1, dimms1, NULL);

  // Sort the energies in parallel using C++17 built in parallel sort
  std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(), 
        [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });

  std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(), 
        [](cfg::line const &a, cfg::line const &b) {return (a.l2 < b.l2) && (a.l1==b.l1);});

  // L initial sort & save
  Li_dat.reserve(t_sz);
  en_srtd.reserve(t_sz);
  idx_data.reserve(t_sz);
  hsize_t write_sz[] = {t_sz};
  hsize_t idx_sz[] = {t_sz*4};

  filename = pot + std::to_string(0)+std::to_string(1) + gauge + ".h5";
  file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
  e_space = e1->getSpace();
  e_space.selectHyperslab(H5S_SELECT_SET, count1, offset, stride, block);
  e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l1, e_space);

  e_space = e1->getSpace();
  e_space.selectHyperslab(H5S_SELECT_SET, count2, offset, stride, block);
  e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l2, e_space);

  int last_l1=0, last_l2=0;
  for(const auto &line : cfgs) {
    l1=line.l1;
    l2=line.l2;

    if(l1!=last_l1) {
      filename = pot + std::to_string(l1)+std::to_string(l1+1) + gauge + ".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count1, offset, stride, block);
      e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l1, e_space);
    }

    if(l2!=last_l2) {
      filename = pot + std::to_string(l2)+std::to_string(l2+1) + gauge + ".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count2, offset, stride, block);
      e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l2, e_space);
    }
    std::cout << line.n1 << "\n";
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
  std::sort(std::execution::par_unseq, Li_dat.begin(), Li_dat.end(), 
        [](en_data const &a, en_data const &b) { return a.en < b.en; });

  // Separate the sorted energies from the sorted indices
  for(auto i : Li_dat) {
    en_srtd.emplace_back(i.en);
    idx_elm = {i.n1,i.l1,i.n2,i.l2};
    idx_data.emplace_back(idx_elm);
  }

  // write the energies to a file
  outfile_name = pot + "2_" + std::to_string(0) + std::to_string(1) + gauge + ".h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("e_i", H5::PredType::NATIVE_DOUBLE, 
        H5::DataSpace(1, write_sz))));
  ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);

  en_srtd.clear();

  // write indices to file
  outfile_name = pot + std::to_string(0) + "idx.h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("idx", H5::PredType::NATIVE_INT32, 
        H5::DataSpace(1, idx_sz))));
  ei->write(&idx_data[0], H5::PredType::NATIVE_INT32);

  idx_data.clear();
  Li_dat.clear();

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
  count2[0] =max_Nsz;
  dimms2[0] =max_Nsz;
  memspace_l2.setExtentSimple(1, dimms2, NULL);

  max_line = *std::max_element(cfgs.begin(), cfgs.end(),
        [](cfg::line const &a, cfg::line const &b) { return a.n1 < b.n1; });

  l1_sz = max_line.n1;
  count1[0]=l1_sz;
  dimms1[0]=l1_sz;
  memspace_l1.setExtentSimple(1, dimms1, NULL);

  // Sort the energies in parallel using C++17 built in parallel sort
  std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(), 
        [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });

  // std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(), 
  //       [](cfg::line const &a, cfg::line const &b) {return (a.l2 < b.l2) && (a.l1==b.l1);});

  Lf_dat.reserve(t_sz);
  en_srtd.reserve(t_sz);
  idx_data.reserve(t_sz);
  write_sz[0]=t_sz;
  idx_sz[0]=t_sz*4;

  filename = pot + std::to_string(0)+std::to_string(1) + gauge + ".h5";
  file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
  e_space = e1->getSpace();
  e_space.selectHyperslab(H5S_SELECT_SET, count1, offset, stride, block);
  e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l1, e_space);

  e_space = e1->getSpace();
  e_space.selectHyperslab(H5S_SELECT_SET, count2, offset, stride, block);
  e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l2, e_space);

  last_l1=0; last_l2=0;
  for(const auto &line : cfgs) {
    l1=line.l1;
    l2=line.l2;

    if(l1!=last_l1) {
      filename = pot + std::to_string(l1)+std::to_string(l1+1) + gauge + ".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count1, offset, stride, block);
      e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l1, e_space);
    }

    if(l2!=last_l2) {
      filename = pot + std::to_string(l2)+std::to_string(l2+1) + gauge + ".h5";
      file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
      e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
      e_space = e1->getSpace();
      e_space.selectHyperslab(H5S_SELECT_SET, count2, offset, stride, block);
      e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l2, e_space);
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

  std::sort(std::execution::par_unseq, Lf_dat.begin(), Lf_dat.end(), 
        [](en_data const &a, en_data const &b) { return a.en < b.en; });

  for(auto i : Lf_dat) {
    en_srtd.emplace_back(i.en);
    idx_elm = {i.n1,i.l1,i.n2,i.l2};
    idx_data.emplace_back(idx_elm);
  }

  outfile_name = pot + "2_" + std::to_string(0) + std::to_string(1) + gauge + ".h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("e_f", H5::PredType::NATIVE_DOUBLE, 
        H5::DataSpace(1, write_sz))));
  ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);

  // write indices to file
  outfile_name = pot + std::to_string(1) + "idx.h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("idx", H5::PredType::NATIVE_INT32, 
        H5::DataSpace(1, idx_sz))));
  ei->write(&idx_data[0], H5::PredType::NATIVE_INT32);

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

    cfg::ReadCfg(dir, L_itr, sym, ncf, cfgs);

    t_sz=0;
    for(auto i=0; i<ncf; ++i) {
      t_sz += cfgs[i].n2max-cfgs[i].n2min;
    }

    max_n2l = *std::max_element(cfgs.begin(), cfgs.end(),
          [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
    max_Nsz = max_n2l.n2max;
    en12.reserve(max_Nsz*2);
    count2[0] =max_Nsz;
    dimms2[0] =max_Nsz;
    memspace_l2.setExtentSimple(1, dimms2, NULL);

    max_line = *std::max_element(cfgs.begin(), cfgs.end(),
          [](cfg::line const &a, cfg::line const &b) { return a.n1 < b.n1; });

    l1_sz = max_line.n1;
    count1[0]=l1_sz;
    dimms1[0]=l1_sz;
    memspace_l1.setExtentSimple(1, dimms1, NULL);

    // Sort the energies in parallel using C++17 built in parallel sort
    std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(), 
          [](cfg::line const &a, cfg::line const &b) { return a.l1 < b.l1; });

    // std::sort(std::execution::par_unseq, cfgs.begin(), cfgs.end(), 
    //       [](cfg::line const &a, cfg::line const &b) {return (a.l2 < b.l2) && (a.l1==b.l1);});

    Lf_dat.reserve(t_sz);
    en_srtd.reserve(t_sz);
    idx_data.reserve(t_sz);
    write_sz[0]=t_sz;
    idx_sz[0]=t_sz*4;

    filename = pot + std::to_string(0)+std::to_string(1) + gauge + ".h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
    e_space = e1->getSpace();
    e_space.selectHyperslab(H5S_SELECT_SET, count1, offset, stride, block);
    e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l1, e_space);

    e_space = e1->getSpace();
    e_space.selectHyperslab(H5S_SELECT_SET, count2, offset, stride, block);
    e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l2, e_space);

    last_l1=0; last_l2=0;
    for(const auto &line : cfgs) {
      l1=line.l1;
      l2=line.l2;

      if(l1!=last_l1) {
        filename = pot + std::to_string(l1)+std::to_string(l1+1) + gauge + ".h5";
        file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
        e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
        e_space = e1->getSpace();
        e_space.selectHyperslab(H5S_SELECT_SET, count1, offset, stride, block);
        e1->read(&en12[0], H5::PredType::NATIVE_DOUBLE, memspace_l1, e_space);
      }

      if(l2!=last_l2) {
        filename = pot + std::to_string(l2)+std::to_string(l2+1) + gauge + ".h5";
        file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
        e1 = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("e_i")));
        e_space = e1->getSpace();
        e_space.selectHyperslab(H5S_SELECT_SET, count2, offset, stride, block);
        e1->read(&en12[max_Nsz], H5::PredType::NATIVE_DOUBLE, memspace_l2, e_space);
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

    std::sort(std::execution::par_unseq, Lf_dat.begin(), Lf_dat.end(), 
        [](en_data const &a, en_data const &b) { return a.en < b.en; });

    for(auto i : Lf_dat) {
      en_srtd.emplace_back(i.en);
      idx_elm = {i.n1,i.l1,i.n2,i.l2};
      idx_data.emplace_back(idx_elm);
    }

    ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
          outfile->createDataSet("e_f", H5::PredType::NATIVE_DOUBLE, 
          H5::DataSpace(1, write_sz))));
    ei->write(&en_srtd[0], H5::PredType::NATIVE_DOUBLE);

    // write indices to file
    outfile_name = pot + std::to_string(L_itr) + "idx.h5";
    outfile = std::unique_ptr<H5::H5File>(
              new H5::H5File(outfile_name, H5F_ACC_TRUNC));
    ei = std::unique_ptr<H5::DataSet>(new H5::DataSet(
          outfile->createDataSet("idx", H5::PredType::NATIVE_INT32, 
          H5::DataSpace(1, idx_sz))));
    ei->write(&idx_data[0], H5::PredType::NATIVE_INT32);
  }

  return 0;
}

int DMX2e::calc_dmx(int L_max, std::vector<int> &N_max) {
  int l1i=0, l2i=0, l1f=0, l2f=0;
  int Li_sz=0;
  int Lf_sz=0;
  double Lsq=0.0;
  std::string filename;
  std::string outfile_name;
  std::unique_ptr<H5::H5File> file, outfile;
  std::unique_ptr<H5::DataSet> dmx, L_set;
  H5::DataSpace dmx_space;
  hsize_t offset[2], count[2], stride[2], block[2];
  hsize_t dimms[2], T_dimms[2];
  offset[0]=0;
  offset[1]=0;
  stride[0]=1;
  stride[1]=1;
  block[0] =1;
  block[1] =1;
  H5::DataSpace memspace;

  std::vector<int> dmx_sz(L_max);
  std::vector<dmx_dim> D_dim(L_max);
  std::vector<double> D_data;
  std::vector<double*> D(L_max);
  std::vector<double> T;

  int max_Ldim = *std::max_element(N_max.begin(), N_max.end());
  int max_Lsz = max_Ldim*max_Ldim;
  std::vector<idx4> Li_idx(max_Lsz*L_max);
  std::vector<idx4> Lf_idx(max_Lsz*L_max);
  std::vector<idx4*> buffs(2);

  int tot_sz=0;
  int sz_i=0;
  dmx_dim dim_i;
  // Caluclate the sizes of the 1e dmx subblocks
  for(int i=0; i<L_max; ++i) {
    dim_i = {N_max[i+1],N_max[i]};
    sz_i = N_max[i]*N_max[i+1];
    D_dim[i] = dim_i;
    dmx_sz[i] = sz_i;
    tot_sz += sz_i;
  }

  D_data.reserve(tot_sz);

  D[0] = &D_data[0];
  for(int i=1; i<L_max; ++i) {
    D[i] = &D_data[dmx_sz[i-1]];
  }

  // Read all 1e dipoles
  for(int i=0; i<L_max; ++i) {
    count[0] = D_dim[i].row;
    count[1] = D_dim[i].col;
    dimms[0] = D_dim[i].row;
    dimms[1] = D_dim[i].col;

    filename = pot + std::to_string(i) + std::to_string(i+1) + gauge + ".h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    dmx = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("d_if")));
    dmx_space = dmx->getSpace();
    dmx_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    memspace.setExtentSimple(2, dimms, NULL);
	  dmx->read(D[i], H5::PredType::NATIVE_DOUBLE, memspace, dmx_space);
  }

  // for loop over L's
  int L=0;
  int Lf_i=1;

  // read L indices
  filename = pot + std::to_string(L) + "idx.h5";
  file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
  Li_sz = L_set->getSpace().getSimpleExtentNpoints()/4; // get the size of L
  L_set->read(&Li_idx[0], H5::PredType::NATIVE_INT32);

  filename = pot + std::to_string(Lf_i) + "idx.h5";
  file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
  Lf_sz = L_set->getSpace().getSimpleExtentNpoints()/4; // get the size of L+1
  L_set->read(&Lf_idx[0], H5::PredType::NATIVE_INT32);

  T.reserve(Li_sz*Lf_sz);
  T_dimms[0]=Lf_sz;
  T_dimms[1]=Li_sz;

  // calculate 2e dipoles
  Lsq = 2*sqrt(Lf_i);
  for(int idLf=0; idLf<Lf_sz; ++idLf) { // up to sz Lf
    #pragma omp parallel for private(l1i,l2i,l1f,l2f)
    for(int idLi=0; idLi<Li_sz; ++idLi) { // up to sz Li
      l1i=Li_idx[idLi].l1;
      l2i=Li_idx[idLi].l2;
      l1f=Lf_idx[idLf].l1;
      l2f=Lf_idx[idLf].l2;
      if(l1i+1==l1f && l2i==l2f) {
        T[idLi+Li_sz*idLf] = Lsq*pow(-1,l2i+l1f)*sqrt(4*l1f*l1f-1)*wigner_6j_2e(L,l1f,l1i,l2i)*
              D[l1i][Li_idx[idLi].n1+D_dim[l1i].col*Lf_idx[idLf].n1]/sqrt(l1f);
      } else if(l2i+1==l2f && l1i==l1f) {
        T[idLi+Li_sz*idLf] = Lsq*pow(-1,l1i+l2f)*sqrt(4*l2f*l2f-1)*wigner_6j_2e(L,l2f,l2i,l1i)*
              D[l2i][Li_idx[idLi].n2+D_dim[l2i].col*Lf_idx[idLf].n2]/sqrt(l2f);
      } else {T[idLi+Li_sz*idLf]=0.0;}
    }
  }

  // Save TL1;L2
  outfile_name = pot + "2_" + std::to_string(L) + std::to_string(Lf_i) + gauge + ".h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  dmx = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("d_if", H5::PredType::NATIVE_DOUBLE, 
        H5::DataSpace(2, T_dimms))));
  dmx->write(&T[0], H5::PredType::NATIVE_DOUBLE);

  T.clear();

  // Set L=1, L+1=2 and change index buffer pointers
  int buf_Li=1;
  int buf_Lf=0;
  Li_sz = Lf_sz;
  buffs[0] = &Li_idx[0];
  buffs[1] = &Lf_idx[0];

  // Loop for L>0 
  for(L=1; L<L_max; ++L) {
    Lf_i=L+1;

    // read L indices
    filename = pot + std::to_string(Lf_i) + "idx.h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
    Lf_sz = L_set->getSpace().getSimpleExtentNpoints()/4;
    L_set->read(buffs.at(buf_Lf), H5::PredType::NATIVE_INT32);

    T.reserve(Li_sz*Lf_sz);
    T_dimms[0]=Lf_sz;
    T_dimms[1]=Li_sz;

    // calculate 2e dipoles
    Lsq = 2*sqrt(Lf_i);
    for(int idLf=0; idLf<Lf_sz; ++idLf) { // up to sz Lf
      #pragma omp parallel for private(l1i,l2i,l1f,l2f)
      for(int idLi=0; idLi<Li_sz; ++idLi) { // up to sz Li
        l1i=buffs.at(buf_Li)[idLi].l1;
        l2i=buffs.at(buf_Li)[idLi].l2;
        l1f=buffs.at(buf_Lf)[idLf].l1;
        l2f=buffs.at(buf_Lf)[idLf].l2;
        if(l1i+1==l1f && l2i==l2f) {
          T[idLi+Li_sz*idLf] = Lsq*pow(-1,l2i+l1f)*sqrt(4*l1f*l1f-1)*wigner_6j_2e(L,l1f,l1i,l2i)*
              D[l1i][buffs.at(buf_Li)[idLi].n1+D_dim[l1i].col*buffs.at(buf_Lf)[idLf].n1]/sqrt(l1f);
        } else if(l2i+1==l2f && l1i==l1f) {
          T[idLi+Li_sz*idLf] = Lsq*pow(-1,l1i+l2f)*sqrt(4*l2f*l2f-1)*wigner_6j_2e(L,l2f,l2i,l1i)*
              D[l2i][buffs.at(buf_Li)[idLi].n2+D_dim[l2i].col*buffs.at(buf_Lf)[idLf].n2]/sqrt(l2f);
        } else {T[idLi+Li_sz*idLf]=0.0;}
      }
    }

    outfile_name = pot + "2_" + std::to_string(L) + std::to_string(Lf_i) + gauge + ".h5";
    outfile = std::unique_ptr<H5::H5File>(
              new H5::H5File(outfile_name, H5F_ACC_TRUNC));
    dmx = std::unique_ptr<H5::DataSet>(new H5::DataSet(
          outfile->createDataSet("d_if", H5::PredType::NATIVE_DOUBLE, 
          H5::DataSpace(2, T_dimms))));
    dmx->write(&T[0], H5::PredType::NATIVE_DOUBLE);

    T.clear();

    buf_Li = buf_Lf;
    buf_Lf = 1 - buf_Lf;
    Li_sz = Lf_sz;
  }
  return 0;
}

int DMX2e::calc_dmx(int L_max, std::string dir) {
  int l1i=0, l2i=0, l1f=0, l2f=0;
  int Li_sz=0;
  int Lf_sz=0;
  double Lsq=0.0;
  std::string filename;
  std::string outfile_name;
  std::unique_ptr<H5::H5File> file=nullptr, outfile=nullptr;
  std::unique_ptr<H5::DataSet> dmx=nullptr, L_set=nullptr;
  H5::DataSpace dmx_space;
  hsize_t offset[2], count[2], stride[2], block[2];
  hsize_t dimms[2], T_dimms[2];
  offset[0]=0;
  offset[1]=0;
  stride[0]=1;
  stride[1]=1;
  block[0] =1;
  block[1] =1;
  H5::DataSpace memspace;

  std::vector<int> dmx_sz(L_max);
  std::vector<double> D_data;
  std::vector<double> T;

  int ncf, sym;
  std::vector<cfg::line> cfgs;
  cfg::ReadCfg(dir, 0, sym, ncf, cfgs);

  std::vector<idx4*> buffs(2);

  int t_sz=0;
  for(auto i=0; i<ncf; ++i) {
    t_sz += cfgs[i].n2max-cfgs[i].n2min;
  }

  auto max_n2l = *std::max_element(cfgs.begin(), cfgs.end(),
        [](cfg::line const &a, cfg::line const &b) { return a.n2max < b.n2max; });
  auto max_N = max_n2l.n2max;
  auto max_Nsz = ncf*max_N;
  std::vector<idx4> Li_idx(max_Nsz*L_max);
  std::vector<idx4> Lf_idx(max_Nsz*L_max);

  auto max_line = *std::max_element(cfgs.begin(), cfgs.end(),
        [](cfg::line const &a, cfg::line const &b) { return a.l2 < b.l2; });
  auto l_m = max_line.l2;

  count[0]=max_N;
  count[1]=max_N;
  dimms[0]=max_N;
  dimms[1]=max_N;
  D_data.reserve(max_N*max_N*l_m);
  auto max2 = max_N*max_N;

  // Read all 1e dipoles
  for(auto i=0; i<l_m; ++i) {
    filename = pot + std::to_string(i) + std::to_string(i+1) + gauge + ".h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    dmx = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("d_if")));
    dmx_space = dmx->getSpace();
    dmx_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
    memspace.setExtentSimple(2, dimms, NULL);
	  dmx->read(&D_data[i*max2], H5::PredType::NATIVE_DOUBLE, memspace, dmx_space);
  }

  // for loop over L's
  int L=0;
  int Lf_i=1;

  // read L indices
  filename = pot + std::to_string(L) + "idx.h5";
  file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
  Li_sz = L_set->getSpace().getSimpleExtentNpoints()/4; // get the size of L
  L_set->read(&Li_idx[0], H5::PredType::NATIVE_INT32);

  filename = pot + std::to_string(Lf_i) + "idx.h5";
  file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
  L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
  Lf_sz = L_set->getSpace().getSimpleExtentNpoints()/4; // get the size of L+1
  L_set->read(&Lf_idx[0], H5::PredType::NATIVE_INT32);

  T.reserve(Li_sz*Lf_sz);
  T_dimms[0]=Lf_sz;
  T_dimms[1]=Li_sz;

  // calculate 2e dipoles
  Lsq = 2*sqrt(Lf_i);
  #pragma omp parallel
  {
  for(int idLf=0; idLf<Lf_sz; ++idLf) { // up to sz Lf
    #pragma omp for private(l1i,l2i,l1f,l2f)
    for(int idLi=0; idLi<Li_sz; ++idLi) { // up to sz Li
      l1i=Li_idx[idLi].l1;
      l2i=Li_idx[idLi].l2;
      l1f=Lf_idx[idLf].l1;
      l2f=Lf_idx[idLf].l2;
      if(l1i+1==l1f && l2i==l2f) {
        T[idLi+Li_sz*idLf] = Lsq*pow(-1,l2i+l1f)*sqrt(4*l1f*l1f-1)*wigner_6j_2e(L,l1f,l1i,l2i)*
              D_data[l1i*max2 + Li_idx[idLi].n1+max_N*Lf_idx[idLf].n1]/sqrt(l1f);
      } else if(l2i+1==l2f && l1i==l1f) {
        T[idLi+Li_sz*idLf] = Lsq*pow(-1,l1i+l2f)*sqrt(4*l2f*l2f-1)*wigner_6j_2e(L,l2f,l2i,l1i)*
              D_data[l2i*max2 + Li_idx[idLi].n2+max_N*Lf_idx[idLf].n2]/sqrt(l2f);
      } else {T[idLi+Li_sz*idLf]=0.0;}
    }
  }
  }
  // Save TL1;L2
  outfile_name = pot + "2_" + std::to_string(L) + std::to_string(Lf_i) + gauge + ".h5";
  outfile = std::unique_ptr<H5::H5File>(
            new H5::H5File(outfile_name, H5F_ACC_TRUNC));
  dmx = std::unique_ptr<H5::DataSet>(new H5::DataSet(
        outfile->createDataSet("d_if", H5::PredType::NATIVE_DOUBLE, 
        H5::DataSpace(2, T_dimms))));
  dmx->write(&T[0], H5::PredType::NATIVE_DOUBLE);

  T.clear();

  // Set L=1, L+1=2 and change index buffer pointers
  int buf_Li=1;
  int buf_Lf=0;
  Li_sz = Lf_sz;
  buffs[0] = &Li_idx[0];
  buffs[1] = &Lf_idx[0];

  // Loop for L>0 
  for(Lf_i=2; Lf_i<=L_max; ++Lf_i) {
    L=Lf_i-1;

    // read L indices
    filename = pot + std::to_string(Lf_i) + "idx.h5";
    file = std::unique_ptr<H5::H5File>(new H5::H5File(filename, H5F_ACC_RDONLY));
    L_set = std::unique_ptr<H5::DataSet>(new H5::DataSet(file->openDataSet("idx")));
    Lf_sz = L_set->getSpace().getSimpleExtentNpoints()/4;
    L_set->read(buffs.at(buf_Lf), H5::PredType::NATIVE_INT32);

    T.reserve(Li_sz*Lf_sz);
    T_dimms[0]=Lf_sz;
    T_dimms[1]=Li_sz;

    // calculate 2e dipoles
    Lsq = 2*sqrt(Lf_i);
    #pragma omp parallel
    {
    for(int idLf=0; idLf<Lf_sz; ++idLf) { // up to sz Lf
      #pragma omp for private(l1i,l2i,l1f,l2f)
      for(int idLi=0; idLi<Li_sz; ++idLi) { // up to sz Li
        l1i=buffs.at(buf_Li)[idLi].l1;
        l2i=buffs.at(buf_Li)[idLi].l2;
        l1f=buffs.at(buf_Lf)[idLf].l1;
        l2f=buffs.at(buf_Lf)[idLf].l2;
        if(l1i+1==l1f && l2i==l2f) {
          T[idLi+Li_sz*idLf] = Lsq*pow(-1,l2i+l1f)*sqrt(4*l1f*l1f-1)*wigner_6j_2e(L,l1f,l1i,l2i)*
              D_data[l1i*max2 + buffs.at(buf_Li)[idLi].n1+max_N*buffs.at(buf_Lf)[idLf].n1]/sqrt(l1f);
        } else if(l2i+1==l2f && l1i==l1f) {
          T[idLi+Li_sz*idLf] = Lsq*pow(-1,l1i+l2f)*sqrt(4*l2f*l2f-1)*wigner_6j_2e(L,l2f,l2i,l1i)*
              D_data[l2i*max2 + buffs.at(buf_Li)[idLi].n2+max_N*buffs.at(buf_Lf)[idLf].n2]/sqrt(l2f);
        } else {T[idLi+Li_sz*idLf]=0.0;}
      }
    }
    }

    outfile_name = pot + "2_" + std::to_string(L) + std::to_string(Lf_i) + gauge + ".h5";
    outfile = std::unique_ptr<H5::H5File>(
              new H5::H5File(outfile_name, H5F_ACC_TRUNC));
    dmx = std::unique_ptr<H5::DataSet>(new H5::DataSet(
          outfile->createDataSet("d_if", H5::PredType::NATIVE_DOUBLE, 
          H5::DataSpace(2, T_dimms))));
    dmx->write(&T[0], H5::PredType::NATIVE_DOUBLE);

    T.clear();

    buf_Li = buf_Lf;
    buf_Lf = 1 - buf_Lf;
    Li_sz = Lf_sz;
  }
  return 0;
}

DMX2e::DMX2e(std::string cpot, char gau, int L_max, std::vector<int> &N_max) {
  pot = cpot;
  gauge = gau;

  sort_L(L_max, N_max);

  calc_dmx(L_max, N_max);
}

DMX2e::DMX2e(std::string cpot, char gau, int L_max, std::string inp_dir) {
  pot = cpot;
  gauge = gau;

  sort_L(L_max, inp_dir);

  calc_dmx(L_max, inp_dir);
}