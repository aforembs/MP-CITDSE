#include "dmx_calc.h"
#include <iostream>

int calculate_2edmx(unsigned int lb_max) {
  unsigned int la=0,lb=0, lc=0, ld=1;
  unsigned int L=0;
  double min_sq=0;
  double L_coeff_ac=0;
  double L_coeff_bd=0;
  int ni=10, nf=100;
  int dsz=ni*nf;
  unsigned int la_max = lb_max-1;
  char gauge='v';
  std::string pot="he";
  std::string filename;
  std::string outfile_name;
  std::string setname;

  H5::DataSet dmx;
  H5::DataSpace dmx_space;
  hsize_t offset[2], count[2], stride[2], block[2];
  hsize_t dimms[2];
  offset[0]=0;
  offset[1]=0;
  count[0] =nf; // no. rows
  count[1] =ni;  // no. coumns
  stride[0]=1;
  stride[1]=1;
  block[0] =1;
  block[1] =1;
  dimms[0] =nf;
  dimms[1] =ni;
  H5::DataSpace memspace(2, dimms, NULL);
  std::vector<double> Dac(dsz);
  std::vector<double> Dbd(dsz);
  std::vector<double> T(dsz);
  H5::H5File file;
  H5::H5File *outfile;
  H5::DataSet T_abcd;

  for(la=0; la<la_max; ++la) {
    lb=la;
    L=la+lb;
    lc=la;
    ld=lb+1;
    min_sq = pow(-1,L+1)*sqrt(L+1);

    // Read Dbd
    filename = pot + std::to_string(lb) + std::to_string(ld) + gauge + ".h5";
    file.openFile(filename, H5F_ACC_RDONLY);
    dmx = file.openDataSet("d_if");
    dmx_space = dmx.getSpace();
    dmx_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
		dmx.read(&Dbd[0], H5::PredType::NATIVE_DOUBLE, memspace, dmx_space);

    // wigner
    L_coeff_bd = min_sq*sqrt((2*lb+1)*(2*ld+1)/ld)*wigner_6j_2e(L,ld,lb,la);
    // dbd only
    for(int n=0; n<dsz; ++n) {
      T[n] = L_coeff_bd*Dbd[n];
    }

    // save Tab;cd to file
    setname = "T_" + std::to_string(la) + std::to_string(lb) + 
              "-" + std::to_string(lc) + std::to_string(ld);
    outfile_name = setname + ".h5";
    outfile = new H5::H5File(outfile_name, H5F_ACC_TRUNC);
    T_abcd = outfile->createDataSet(setname, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, count));
    T_abcd.write(&T[0], H5::PredType::NATIVE_DOUBLE);
    outfile->close();

    for(lb=la+1; lb<lb_max; ++lb) {
      lc=la;
      ld=lb+1;
      L=la+lb;
      min_sq = pow(-1,L+1)*sqrt(L+1);

      // Read Dbd
      filename = pot + std::to_string(lb) + std::to_string(ld) + gauge + ".h5";
      file.openFile(filename, H5F_ACC_RDONLY);
      dmx = file.openDataSet("d_if");
      dmx_space = dmx.getSpace();
      dmx_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
		  dmx.read(&Dbd[0], H5::PredType::NATIVE_DOUBLE, memspace, dmx_space);

      // wigner
      L_coeff_bd = min_sq*sqrt((2*lb+1)*(2*ld+1)/ld)*wigner_6j_2e(L,ld,lb,la);
      // dbd only
      for(int n=0; n<dsz; ++n) {
        T[n] = L_coeff_bd*Dbd[n];
      }

      // save Tab;cd to file
      setname = "T_" + std::to_string(la) + std::to_string(lb) + 
                "-" + std::to_string(lc) + std::to_string(ld);
      outfile_name = setname + ".h5";
      outfile = new H5::H5File(outfile_name, H5F_ACC_TRUNC);
      T_abcd = outfile->createDataSet(setname, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, count));
      T_abcd.write(&T[0], H5::PredType::NATIVE_DOUBLE);
      outfile->close();

      lc=la+1;
      ld-=1;
      // Read Dac
      filename = pot + std::to_string(la) + std::to_string(lc) + gauge + ".h5";
      file.openFile(filename, H5F_ACC_RDONLY);
      dmx = file.openDataSet("d_if");
      dmx_space = dmx.getSpace();
      dmx_space.selectHyperslab(H5S_SELECT_SET, count, offset, stride, block);
		  dmx.read(&Dac[0], H5::PredType::NATIVE_DOUBLE, memspace, dmx_space);

      // wigner
      L_coeff_ac = min_sq*sqrt((2*la+1)*(2*lc+1)/lb)*wigner_6j_2e(L,lc,la,lb);
      // dac + dbd
      for(int n=0; n<dsz; ++n) {
        T[n] = L_coeff_ac*Dac[n];
      }

      // save Tab;cd to file
      setname = "T_" + std::to_string(la) + std::to_string(lb) + 
                "-" + std::to_string(lc) + std::to_string(ld);
      outfile_name = setname + ".h5";
      outfile = new H5::H5File(outfile_name, H5F_ACC_TRUNC);
      T_abcd = outfile->createDataSet(setname, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(2, count));
      T_abcd.write(&T[0], H5::PredType::NATIVE_DOUBLE);
      outfile->close();
    }
  }
  return 0;
}