#include <fstream>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>
#include <stdexcept>

template <typename T>
void writeVectorToFile(const std::vector<T>& v, const std::string& filename, bool append = false)
{
  static std::ofstream file;
  if (!file.is_open() || file.good() == false || file.rdstate()) {
    // Open the file in the appropriate mode
    if (append) {
      file.open(filename, std::ios_base::app);
    } else {
      file.open(filename);
    }

    if (!file.is_open()) {
      throw std::runtime_error("Could not open file " + filename + " for writing!");
    }
  }

  // Write the vector to the file using std::copy
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(file, "\n"));

  // Check error status and close the file if it needs to be closed
  if (!append || file.rdstate()) {
    file.close();
  }
}

template <typename T>
std::vector<T> generateGaussianVector(std::size_t size, T mean, T variance)
{
    // Create a normal distribution with the specified mean and variance
    std::normal_distribution<T> distribution(mean, variance);

    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    // Generate a vector of random values using std::transform
    std::vector<T> result(size);
    std::transform(result.begin(), result.end(), result.begin(),
                   [&distribution, &gen](T& elem) { return distribution(gen); });

    return result;
}

template<typename T, typename Func>
void transformVector(std::vector<T>& vec, const T& constant, Func op)
{
  std::transform(vec.begin(), vec.end(), vec.begin(),
                 [constant, op](const T& elem) { return op(elem, constant); });
}

template <typename T>
void sumVectors(const std::vector<T>& vec1, const std::vector<T>& vec2, std::vector<T>& result)
{
  // std::size_t size = vec1.size();
  // result.resize(size);

  std::transform(vec1.begin(), vec1.end(),
                 vec2.begin(), result.begin(),
                 [](const T& a, const T& b) { return a + b; });
}

template <typename T, typename U>
void dotProduct(const std::vector<std::vector<T>>& matrix, const std::vector<U>& vec, std::vector<U>& result)
{
    // Check that dimensions of matrix and vec are compatible
    // if (matrix[0].size() != vec.size()) {
    //     throw std::invalid_argument("Dimensions of matrix and vector are not compatible!");
    // }

    // Compute the dot product
    // std::size_t size = matrix.size();
    // result.resize(size);

    for (std::size_t i = 0; i < matrix.size(); ++i) {
        U sum = std::inner_product(matrix[i].begin(), matrix[i].end(), vec.begin(), U());
        result[i] = sum;
    }
}

template <typename T>
void transferFunc(const std::vector<T>& x, T thresh, std::vector<T>& result)
{
    // std::vector<T> result(x.size());
    std::transform(x.begin(), x.end(), result.begin(), [=](T val) {
        return thresh * (1.0 + std::erf(val / std::sqrt(2.0)));
    });
    // return result;
}
