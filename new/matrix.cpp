#include "matrix.hpp"

std::vector<std::vector<double>> generateMatrix(int N) {
    double J0 = 1.0; // modify as desired
    double J1 = 1.0; // modify as desired

    std::vector<double> cosines(N/2+1);
    for (int k = 0; k <= N/2; k++) {
        double angle = 2.0 / M_PI / k / N;
        cosines[k] = std::cos(angle);
    }

    std::vector<std::vector<double>> C(N, std::vector<double>(N, 0.0));
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            int k = std::abs(j-i);
            if (k > N/2) {
                k = N - k;
            }
            C[i][j] = J0/N / (1.0 + J1 / cosines[k]);
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            C[i][j] = C[j][i];
        }
    }
    return C;
}
