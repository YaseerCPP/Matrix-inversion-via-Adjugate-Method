#include <iostream>
#include <vector>
#include <iomanip>

// Function to calculate the determinant of a matrix
double determinant(const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    if (n == 1) return matrix[0][0];
    if (n == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    double det = 0.0;
    for (int p = 0; p < n; ++p) {
        std::vector<std::vector<double>> submatrix(n - 1, std::vector<double>(n - 1));
        for (int i = 1; i < n; ++i) {
            int subcol = 0;
            for (int j = 0; j < n; ++j) {
                if (j == p) continue;
                submatrix[i - 1][subcol++] = matrix[i][j];
            }
        }
        det += matrix[0][p] * determinant(submatrix) * (p % 2 == 0 ? 1 : -1);
    }
    return det;
}

// Function to calculate the cofactor matrix
std::vector<std::vector<double>> cofactorMatrix(const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    std::vector<std::vector<double>> cofactor(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::vector<std::vector<double>> submatrix(n - 1, std::vector<double>(n - 1));
            int subi = 0;
            for (int k = 0; k < n; ++k) {
                if (k == i) continue;
                int subj = 0;
                for (int l = 0; l < n; ++l) {
                    if (l == j) continue;
                    submatrix[subi][subj++] = matrix[k][l];
                }
                ++subi;
            }
            cofactor[i][j] = determinant(submatrix) * ((i + j) % 2 == 0 ? 1 : -1);
        }
    }
    return cofactor;
}

// Function to transpose a matrix
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    std::vector<std::vector<double>> transposed(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            transposed[i][j] = matrix[j][i];
        }
    }
    return transposed;
}

// Function to calculate the inverse of a matrix using the adjugate method
bool inverseMatrix(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& invMatrix) {
    double det = determinant(matrix);
    if (det == 0.0) return false; // Singular matrix

    std::vector<std::vector<double>> cofactorMat = cofactorMatrix(matrix);
    std::vector<std::vector<double>> adjugateMat = transpose(cofactorMat);

    int n = matrix.size();
    invMatrix = std::vector<std::vector<double>>(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            invMatrix[i][j] = adjugateMat[i][j] / det;
        }
    }
    return true;
}

// Function to print the matrix
void printMatrix(const std::vector<std::vector<double>>& matrix) {
    int n = matrix.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << std::setw(10) << std::setprecision(4) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main() {
    // Example matrix (3x3)
    std::vector<std::vector<double>> A = {
        {4, 3, 2},
        {3, 7, 1},
        {2, 5, 3}
    };

    std::cout << "Original Matrix:" << std::endl;
    printMatrix(A);

    std::vector<std::vector<double>> invA;
    if (inverseMatrix(A, invA)) {
        std::cout << "Inverse Matrix:" << std::endl;
        printMatrix(invA);
    } else {
        std::cout << "Matrix is singular and cannot be inverted." << std::endl;
    }

    return 0;
}
