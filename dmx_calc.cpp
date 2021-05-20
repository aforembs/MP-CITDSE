#include "dmx_calc.h"
#include <iostream>
#include <iomanip>

en_L DMX2e::make_enL(uint L) {
  int sz = (L+2)/2;
  uint mn;
  auto lp = std::vector<l_ab>(sz);
  for(auto i=0; i<sz; ++i) {
    lp[i].l1=i;
    lp[i].l2=L-i;
  }
  auto en_vec = std::vector<en_data>();
  en_L Len = {L,mn,lp,en_vec};
  return Len;
}

int DMX2e::sort_L(en_L &Lif, uint N_sz) {
  uint l1, l2;
  std::string filename;
  std::string outfile_name;
  H5::DataSet e1, ei;
  H5::H5File file;
  H5::H5File *outfile;
  std::vector<double> en12(N_sz*2);
  std::vector<double> en_strd;
  en_strd.reserve(N_sz);
  hsize_t write_sz[] = {N_sz};
  en_data en_d;
  int num_l_pairs = Lif.l_pair.size();

  uint t_sz=0;
  for(auto i=0; i<num_l_pairs; ++i) {
    if (l1!=l2) {
      t_sz+=N_sz*N_sz;
    } else {
      t_sz+=N_sz*(N_sz+1)/2;
    }
  }
  Lif.en_dat.reserve(t_sz);

  for(auto i=0; i<num_l_pairs; ++i) {
    l1=Lif.l_pair[i].l1;
    l2=Lif.l_pair[i].l2;

    if(l1==l2) {
      filename = pot + std::to_string(l1)+std::to_string(l1+1) + gauge + ".h5";
      file.openFile(filename, H5F_ACC_RDONLY);
      e1 = file.openDataSet("e_i");
      e1.read(&en12[0], H5::PredType::NATIVE_DOUBLE);
      file.close();

      double ent;
      for(uint n1=0; n1<N_sz; ++n1) {
        for(uint n2=n1; n2<N_sz; ++n2) {
          ent = en12[n1] + en12[n2];
          en_d = {ent,n1,l1,n2,l2};
          Lif.en_dat.emplace_back(en_d);
        }
      }

    } else {
      Lif.en_dat.reserve(N_sz*N_sz);
      filename = pot + std::to_string(l1)+std::to_string(l1+1) + gauge + ".h5";
      file.openFile(filename, H5F_ACC_RDONLY);
      e1 = file.openDataSet("e_i");
      e1.read(&en12[0], H5::PredType::NATIVE_DOUBLE);
      file.close();

      filename = pot + std::to_string(l2)+std::to_string(l2+1) + gauge + ".h5";
      file.openFile(filename, H5F_ACC_RDONLY);
      e1 = file.openDataSet("e_i");
      e1.read(&en12[N_sz], H5::PredType::NATIVE_DOUBLE);
      file.close();

      double ent;
      for(uint n1=0; n1<N_sz; ++n1) {
        for(uint n2=0; n2<N_sz; ++n2) {
          ent = en12[n1] + en12[n2];
          en_d = {ent,n1,l1,n2,l2};
          Lif.en_dat.emplace_back(en_d);
        }
      }
    }
  }

  std::sort(std::execution::par_unseq, Lif.en_dat.begin(), Lif.en_dat.end(), 
        [](en_data const &a, en_data const &b) { return a.en < b.en; });

  Lif.en_dat.erase(Lif.en_dat.begin()+N_sz, Lif.en_dat.end());

  auto max = std::max_element(std::execution::par_unseq, Lif.en_dat.begin(), Lif.en_dat.end(),
    [](en_data const &a, en_data const &b) { return a.n1 < b.n1; });

  Lif.max_n = std::max(max->n1, std::max_element(std::execution::par_unseq, Lif.en_dat.begin(), Lif.en_dat.end(),
    [](en_data const &a, en_data const &b) { return a.n2 < b.n2; })->n2);

  for(auto i : Lif.en_dat) {
    en_strd.emplace_back(i.en);
    //std::cout << std::setprecision(15) << i.en << ", " << i.n1 << ", " << i.l1 << ", " << i.n2 << ", " << i.l2 << "\n";
  }

  outfile_name = pot + "2_" + std::to_string(Lif.L) + std::to_string(Lif.L+1) + gauge + ".h5";
  outfile = new H5::H5File(outfile_name, H5F_ACC_TRUNC);
  ei = outfile->createDataSet("e_i", H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, write_sz));
  ei.write(&en_strd[0], H5::PredType::NATIVE_DOUBLE);
  outfile->close();

  return 0;
}

DMX2e::DMX2e(std::string cpot, char gau, uint L_max, uint N_max) {
  pot = cpot;
  gauge = gau;
  en_L L_struct;

  for(uint Li=0; Li<=L_max; ++Li) {
    L_struct = make_enL(Li);
    sort_L(L_struct, N_max);
    L_struct.en_dat.clear();
  }
}

// int calculate_2edmx(uint L_max, uint Ni_sz, uint Nf_sz) {
//   uint la=0,lb=0, lc=0, ld=1;
//   uint L=la+lb;
//   double min_sq=0;
//   double L_coeff_ac=0;
//   double L_coeff_bd=0;
//   int dsz=Ni_sz*Nf_sz;
//   char gauge='v';
//   std::string pot="he";
//   std::string filename;
//   std::string outfile_name;
//   std::string setname;

//   H5::DataSet dmx;
//   H5::DataSpace dmx_space;
//   hsize_t offset[2], count[2], stride[2], block[2];
//   hsize_t dimms[2];
//   offset[0]=0;
//   offset[1]=0;
//   count[0] =Nf_sz; // no. rows
//   count[1] =Ni_sz;  // no. coumns
//   stride[0]=1;
//   stride[1]=1;
//   block[0] =1;
//   block[1] =1;
//   dimms[0] =Nf_sz;
//   dimms[1] =Ni_sz;
//   H5::DataSpace memspace(2, dimms, NULL);
//   std::vector<double> Dac(dsz);
//   std::vector<double> Dbd(dsz);
//   std::vector<double> T;
//   T.reserve();
//   H5::H5File file;
//   H5::H5File *outfile;
//   H5::DataSet T_abcd;

//   min_sq = 2*pow(-1,L+1)*sqrt(L+1); // 2 for antisymmetrisation

//   // by L
//   filename = pot + std::to_string(lb) + std::to_string(ld) + gauge + ".h5";
//   file.openFile(filename, H5F_ACC_RDONLY);
//   dmx = file.openDataSet("d_if");
//   dmx_space = dmx.getSpace();
//   dmx_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
// 	dmx.read(&Dbd[0], H5::PredType::NATIVE_DOUBLE, memspace, dmx_space);

//   for(L=1; L<L_max; ++L) {
//     if ()
//     for (lb=0; lb<)
//     for (la=0; la<=lb; ++la) {

//     }
//   }




//   // need to iterate over L not la lb
//   for(la=0; la<la_max; ++la) {
//     lb=la;
//     L=la+lb;
//     lc=la;
//     ld=lb+1;
//     min_sq = 2*pow(-1,L+1)*sqrt(L+1); // 2 for antisymmetrisation

//     // Read Dbd
//     filename = pot + std::to_string(lb) + std::to_string(ld) + gauge + ".h5";
//     file.openFile(filename, H5F_ACC_RDONLY);
//     dmx = file.openDataSet("d_if");
//     dmx_space = dmx.getSpace();
//     dmx_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
// 		dmx.read(&Dbd[0], H5::PredType::NATIVE_DOUBLE, memspace, dmx_space);

//     // wigner
//     L_coeff_bd = min_sq*sqrt((2*lb+1)*(2*ld+1)/ld)*wigner_6j_2e(L,ld,lb,la);
//     // dbd only
//     for(int n=0; n<dsz; ++n) {
//       T[n] = L_coeff_bd*Dbd[n];
//     }

//     // save Tab;cd to file
//     setname = pot + std::to_string(la) + std::to_string(lb) + 
//               "-" + std::to_string(lc) + std::to_string(ld);
//     outfile_name = setname + ".h5";
//     outfile = new H5::H5File(outfile_name, H5F_ACC_TRUNC);
//     T_abcd = outfile->createDataSet(setname, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, count));
//     T_abcd.write(&T[0], H5::PredType::NATIVE_DOUBLE);
//     outfile->close();

//     for(lb=la+1; lb<lb_max; ++lb) {
//       lc=la;
//       ld=lb+1;
//       L=la+lb;
//       min_sq = 2*pow(-1,L+1)*sqrt(L+1);

//       // Read Dbd
//       filename = pot + std::to_string(lb) + std::to_string(ld) + gauge + ".h5";
//       file.openFile(filename, H5F_ACC_RDONLY);
//       dmx = file.openDataSet("d_if");
//       dmx_space = dmx.getSpace();
//       dmx_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
// 		  dmx.read(&Dbd[0], H5::PredType::NATIVE_DOUBLE, memspace, dmx_space);

//       // wigner
//       L_coeff_bd = min_sq*sqrt((2*lb+1)*(2*ld+1)/ld)*wigner_6j_2e(L,ld,lb,la);
//       // dbd only
//       for(int n=0; n<dsz; ++n) {
//         T[n] = L_coeff_bd*Dbd[n];
//       }

//       // save Tab;cd to file
//       setname = pot + std::to_string(la) + std::to_string(lb) + 
//                 "-" + std::to_string(lc) + std::to_string(ld);
//       outfile_name = setname + ".h5";
//       outfile = new H5::H5File(outfile_name, H5F_ACC_TRUNC);
//       T_abcd = outfile->createDataSet(setname, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, count));
//       T_abcd.write(&T[0], H5::PredType::NATIVE_DOUBLE);
//       outfile->close();

//       lc=la+1;
//       ld-=1;
//       // Read Dac
//       filename = pot + std::to_string(la) + std::to_string(lc) + gauge + ".h5";
//       file.openFile(filename, H5F_ACC_RDONLY);
//       dmx = file.openDataSet("d_if");
//       dmx_space = dmx.getSpace();
//       dmx_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
// 		  dmx.read(&Dac[0], H5::PredType::NATIVE_DOUBLE, memspace, dmx_space);

//       // wigner
//       L_coeff_ac = min_sq*sqrt((2*la+1)*(2*lc+1)/lb)*wigner_6j_2e(L,lc,la,lb);
//       // dac + dbd
//       for(int n=0; n<dsz; ++n) {
//         T[n] = L_coeff_ac*Dac[n];
//       }

//       // save Tab;cd to file
//       setname = pot + std::to_string(la) + std::to_string(lb) + 
//                 "-" + std::to_string(lc) + std::to_string(ld);
//       outfile_name = setname + ".h5";
//       outfile = new H5::H5File(outfile_name, H5F_ACC_TRUNC);
//       T_abcd = outfile->createDataSet(setname, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, count));
//       T_abcd.write(&T[0], H5::PredType::NATIVE_DOUBLE);
//       outfile->close();
//     }
//   }
//   return 0;
// }