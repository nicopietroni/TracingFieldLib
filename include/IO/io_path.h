#ifndef IO_PATH
#define IO_PATH

#include <IO/file_path_functions.h>
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <point3.h>
#include <sstream>
#include <string>
#include <vector>


inline bool WritePath(const std::string &filename,
                      const std::vector<std::vector<int> > &Path,
                      bool writeInfo = true) {
  if (writeInfo)
    std::cout << "Saving:" << filename.c_str() << std::endl;

  std::ofstream outfile(filename);
  if (!outfile) {
    std::cerr << "Error: Could not open file " << filename.c_str() << std::endl;
    return false; // Failed to open the file
  }

  for (size_t i=0;i<Path.size();i++) {
    outfile << "P " << i << "\n";
    if (!outfile) {
      return false; // Write operation failed
    }
    for (size_t j=0;j<Path[i].size();j++)
      outfile << Path[i][j] << "\n";
  }

  outfile.close();
  return true;
}


inline bool WritePath(const std::string &filename,
                      const std::vector<std::vector<int> > &Path,
                      const Geo::TracingGraph<ScalarType> &TGraph) 
{
  // if (writeInfo)
  //   std::cout << "Saving:" << filename.c_str() << std::endl;

  // std::ofstream outfile(filename);
  // if (!outfile) {
  //   std::cerr << "Error: Could not open file " << filename.c_str() << std::endl;
  //   return false; // Failed to open the file
  // }

  // for (size_t i=0;i<Path.size();i++) {
  //   outfile << "P " << i << "\n";
  //   if (!outfile) {
  //     return false; // Write operation failed
  //   }
  //   for (size_t j=0;j<Path[i].size();j++)
  //     outfile << Path[i][j] << "\n";
  // }

  // outfile.close();
  // return true;
}

#endif
