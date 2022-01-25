#include <unistd.h>
#include <cstdlib>
#include <filesystem>
#include "V12.h"

namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
  std::string opt_file;

  for(;;) {
    switch(getopt(argc, argv, "hf:")) {
      case 'h':
        std::cout << "Program for calculating the inter-electronic interaction\n"
                  << "coefficients <n1l1;n2l2|r_12|n'1l'1;n'2l'2>\n"
                  << "-f <path> yaml input file with the input settings\n";
        return -1;
      case 'f':
        opt_file = optarg; 
        continue;
    }
    break;
  }

  

  return 0;
} 