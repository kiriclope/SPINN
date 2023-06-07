#include <fstream>
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
void saveArrayToFile(std::ofstream& outFile, T& arr, size_t len) {

  if (!outFile.good()) {
    std::cerr << "Error: could not write to file" << std::endl;
    return;
  }

  std::cout << len << std::endl;

  for (size_t i = 0; i < len; ++i)
    outFile << arr[i] << " ";
  outFile << std::endl;
}

template<typename T>
void loadArrayFromFile(std::ifstream& inFile, T& arr, size_t len) {

  if (!inFile.good()) {
    std::cerr << "Error: could not read from file" << std::endl;
    return;
  }

  for (size_t i = 0; i < len; ++i)
    inFile >> arr[i];
}
