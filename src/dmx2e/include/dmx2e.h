#ifndef DMX2E_H_
#define DMX2E_H_

#include <cmath>
#include <vector>
#include <string>
#include <execution>
#include <algorithm>
#include <fstream>
#include <yaml-cpp/yaml.h>
#include "H5Cpp.h"
#include "dmx_typ.h"
#include "wigxjpf.h"
#include "cfg_in.h"

namespace dmx2e {
  int ReadConfig(std::string file, std::string &pot, 
                int &L_max, int &l_max, char &gauge);

  /* Function for sorting the n1,l1+n2,l2 LN states by ascending energies
   * these sorted energies are stored in the he2_<L><L+1>v.h5 files.
   * The n1,l1,n2,l2 indices are stored in the he<L>idx.h5 files
   * @param in L_max - The maximum total orbital angular momentum L
   * @param in dir
   */
  int SortL(std::string cpot, int L_max, char gauge, std::string dir);

  /* Function for calculating the 2e dipole matrix elements.
   * First all of the required subsections of the 1e dmx's are read into a vector.
   * Then for each L, L+1 pair the indices of both L's are read from the he<L>.idx files.
   * Finally the 2e dmx elements are calculated and saved in the he2_<L><L+1>v.h5 files.
   * @param in L_max - The maximum total orbital angular momentum L
   * @param in dir
   */
  int GenDipole(std::string cpot, int L_max, int l_m, char gauge, std::string dir);
};

#endif // DMX2E_H_