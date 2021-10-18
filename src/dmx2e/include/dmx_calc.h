#ifndef DMX_CALC_H_
#define DMX_CALC_H_

#include <cmath>
#include <vector>
#include <string>
#include <execution>
#include <algorithm>
#include <fstream>
#include "H5Cpp.h"
#include "dmx_typ.h"
#include "wigner_sym.h"

class DMX2e {
  private:
    /* Class Variables */
    char gauge;
    std::string pot;

    /* Function for calculating which l1,l2 combinations belong to which L
     * i.e. for L=0 l1,l2 = 0,0 1,1 2,2 ...
     * @param in L       - The 2e total orbital angular momentum L
     * @param in l1e_max - The maximum 1e angular momentum l1 and l2
     */ 
    en_L make_enL(uint L, uint l1e_max);

    /* Function for sorting the n1,l1+n2,l2 LN states by ascending energies
     * these sorted energies are stored in the he2_<L><L+1>v.h5 files.
     * The n1,l1,n2,l2 indices are stored in the he<L>idx.h5 files
     * @param in L_max - The maximum total orbital angular momentum L
     * @param in N_sz  - Vector of the total number of 1e states used from each
     *                   1e angular momentum l1 or l2. The number of elements must 
     *                   equal L_max+1
     */
    int sort_L(uint L_max, std::vector<uint> &N_sz);

    /* Function for calculating the 2e dipole matrix elements.
     * First all of the required subsections of the 1e dmx's are read into a vector.
     * Then for each L, L+1 pair the indices of both L's are read from the he<L>.idx files.
     * Finally the 2e dmx elements are calculated and saved in the he2_<L><L+1>v.h5 files.
     * @param in L_max - The maximum total orbital angular momentum L
     * @param in N_sz  - Vector of the total number of 1e states used from each
     *                   1e angular momentum l1 or l2. The number of elements must 
     *                   equal L_max+1
     */
    int calc_dmx(uint L_max, std::vector<uint> &N_max);

  public:
    /* Constructor for the DMX2e class, calculates and saves the sorted L states
     * and the 2e dmx elements for each L->L+1 transition.
     * @param in cpot - the type of the coulomb potential e.g. "he" for helium
     * @param in gau  - the gauge of the system
     * @param in N_sz  - Vector of the total number of 1e states used from each
     *                   1e angular momentum l1 or l2. The number of elements must 
     *                   equal L_max+1
     */
    DMX2e(std::string cpot, char gau, uint L_max, std::vector<uint> &N_max);

};

// int calculate_2edmx(uint L_max, uint Ni_sz, uint Nf_sz);

#endif // DMX_CALC_H_