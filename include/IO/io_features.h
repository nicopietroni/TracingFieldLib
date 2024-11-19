#ifndef IO_FEATURES
#define IO_FEATURES

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

inline bool LoadClist(const std::string &filename,
                      std::vector<std::pair<int, int>> &Features,
                      bool writeInfo = true) {
  if (writeInfo)
    std::cout << "Loading:" << filename.c_str() << std::endl;
  std::ifstream file(filename.c_str());

  if (!file.is_open()) {
    std::cerr << "Error: Could not open file " << filename.c_str() << std::endl;
    return false;
  }

  std::string line;

  // Read each line from the file
  int CurrV = 0;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::vector<int> numbers;
    int NextV;

    // Split the line by spaces and store each number in the vector
    while (iss >> NextV) {
      Features.push_back(std::pair<int, int>(CurrV, NextV));
      // numbers.push_back(num);
    }

    // Process or output the numbers as needed (example: print the numbers)
    // std::cout << "Read numbers: ";
    // for (const int& number : numbers) {
    //     //std::cout << number << " ";
    // }
    // std::cout << std::endl;
    CurrV++;
  }

  file.close();
  return true;
}

inline bool WriteFeatures(const std::string &filename,
                          const std::vector<std::pair<int, int>> &Features,
                          bool writeInfo = true) {
  if (writeInfo)
    std::cout << "Saving:" << filename.c_str() << std::endl;

  std::ofstream outfile(filename);
  if (!outfile) {
    std::cerr << "Error: Could not open file " << filename.c_str() << std::endl;
    return false; // Failed to open the file
  }

  for (const auto &feature : Features) {
    outfile << feature.first << " " << feature.second << "\n";
    if (!outfile) {
      return false; // Write operation failed
    }
  }

  outfile.close();
  return true;
}

inline bool ReadFeatures(const std::string &filename,
                         std::vector<std::pair<int, int>> &Features,
                         bool writeInfo = true) {

  if (writeInfo)
    std::cout << "Reading:" << filename.c_str() << std::endl;

  std::ifstream infile(filename);
  if (!infile) {
    std::cerr << "Error: Could not open file " << filename.c_str() << std::endl;
    return false; // Failed to open the file
  }

  Features.clear(); // Clear the vector to ensure it's empty before reading new data

  int first, second;
  while (infile >> first >> second) {
    Features.emplace_back(first, second);
  }
  
  std::sort(Features.begin(), Features.end());
  auto last = std::unique(Features.begin(), Features.end());
  Features.erase(last, Features.end());

  if (!infile.eof()) {
    // An error occurred before reaching the end of the file
    return false;
  }

  infile.close();
  return true;
}

#endif
