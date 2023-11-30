#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <fstream>
#include <random>
#include <vector>

#include "globals.hpp"

template <typename T>
float popSum(T*& v, size_t start_idx = 0, size_t end_idx = 0)
{
    if (end_idx == 0) {
        end_idx = N;
    }
    if (start_idx >= end_idx) {
        return 0.0;  // or throw an exception
    }

    T sum = 0.0;
    for(size_t i=start_idx; i<end_idx; ++i)
      sum += v[i];

    return (float) sum ;
}

template <typename T>
float popMean(T*& v, size_t start_idx = 0, size_t end_idx = 0)
{
    if (end_idx == 0) {
        end_idx = N;
    }
    if (start_idx >= end_idx) {
        return 0.0;  // or throw an exception
    }

    T sum = 0.0;
    for(size_t i=start_idx; i<end_idx; ++i)
      sum += v[i];

    return (float) sum / (float) (end_idx - start_idx);
}

template<typename T>
void saveArrayToFile(std::ofstream& outFile, T* arr, size_t len) {
  
  if (!outFile.good()) {
    std::cerr << "Error: could not write to file" << std::endl;
    return;
  }

  outFile.write(reinterpret_cast<char*>(arr), len * sizeof(T));

  // for (size_t i = 0; i < len; ++i)
  //   outFile << arr[i] << " ";
  // outFile << std::endl;
}

template<typename T>
void loadArrayFromFile(std::ifstream& inFile, T* arr, size_t len) {

  if (!inFile.good()) {
    std::cerr << "Error: could not read from file" << std::endl;
    return;
  }

  inFile.read(reinterpret_cast<char*>(arr), len * sizeof(T));

  // Alternatively, if you want to read from file line by line:
  // for (size_t i = 0; i < len; ++i) {
  //   inFile >> arr[i];
  // }
}

template<typename T>
void saveVectorToFile(std::ofstream& outFile, std::vector<T>& vec) {
  
  if (!outFile.good()) {
    std::cerr << "Error: could not write to file" << std::endl;
    return;
  }
  
  outFile.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(T));
}

template<typename T>
void loadVectorFromFile(std::ifstream& inFile, std::vector<T>& vec, size_t len) {

  inFile.seekg(0, std::ios::end);
  size_t fileSize = inFile.tellg();
  inFile.seekg(0, std::ios::beg);
  
  // Calculate the number of elements of type T to read
  size_t numElements = fileSize / sizeof(T);
  
  vec.resize(numElements);
  
  if (!inFile.good()) {
    std::cerr << "Error: could not read from file" << std::endl;
    return;
  }
  
  inFile.read(reinterpret_cast<char*>(vec.data()), fileSize);
}

template<typename T>
void saveSpikes(std::ofstream& file, std::vector<std::pair<int, T>>& spikes) {
  if (file.is_open()) {
    for (const auto& spike : spikes) {
      file << spike.first << "," << spike.second << "\n";
    }
  }
  file.close();
}

template<typename T1, typename T2>
void saveVectorOfPairsToFile(std::ofstream& outFile, const std::vector<std::pair<T1, T2>>& vec) {
  if (!outFile.good()) {
    std::cerr << "Error: could not write to file" << std::endl;
    return;
  }

  for (const auto& element : vec) {
    outFile.write(reinterpret_cast<const char*>(&element.first), sizeof(element.first));
    outFile.write(reinterpret_cast<const char*>(&element.second), sizeof(element.second));
  }
}

#endif
