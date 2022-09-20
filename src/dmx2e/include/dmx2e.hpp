#ifndef DMX2E_H_
#define DMX2E_H_

#include "cfg_in.hpp"
#include "dmx_typ.hpp"
#include "wigxjpf.h"
#include <H5Cpp.h>
#include <algorithm>
#include <cmath>
#include <execution>
#include <fstream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace dmx2e {
int ReadConfig(std::string file, std::string &pot, int &L_max, int &l_max,
               char &gauge);

/* Function for sorting the n1,l1+n2,l2 LN states by ascending energies
 * these sorted energies are stored in the he2_<L><L+1>v.h5 files.
 * The n1,l1,n2,l2 indices are stored in the he<L>idx.h5 files
 * @param in L_max - The maximum total orbital angular momentum L
 * @param in dir
 */
int SortL(std::string cpot, int L_max, char gauge, std::string dir);

/* Function for calculating the 2e dipole matrix elements.
 * First all of the required subsections of the 1e dmx's are read into a vector.
 * Then for each L, L+1 pair the indices of both L's are read from the he<L>.idx
 * files. Finally the 2e dmx elements are calculated and saved in the
 * he2_<L><L+1>v.h5 files.
 * @param in L_max - The maximum total orbital angular momentum L
 * @param in dir
 */
int GenDipole(std::string cpot, int L_max, int l_m, char gauge,
              std::string dir);
}; // namespace dmx2e

#endif // DMX2E_H_