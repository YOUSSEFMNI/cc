#ifndef MATRIX_H
#define MATRIX_H

#include <list>
#include <iostream>

// Macro pour l'import/export de la DLL
#ifdef MATRIX_EXPORTS
#define MATRIX_API __declspec(dllexport)
#else
#define MATRIX_API __declspec(dllimport)
#endif


class Matrix {
private:
    std::list<std::list<double>> mat;
    int rows;
    int cols;

public:
    Matrix(int rows, int cols);
    void input();
    void display() const;
    Matrix inverseGauss() const;
    Matrix inverseCramer() const;
    double determinant(const std::list<std::list<double>>& m, int n) const;
    std::list<std::list<double>> adjoint() const;
    std::list<std::list<double>> getCofactor(const std::list<std::list<double>>& m, int p, int q, int n) const;
    Matrix multiply(const Matrix& other) const;
    bool isIdentity() const;
    bool isSquare() const { return rows == cols; }
    void roundSmallValues(double threshold = 1e-8);
    Matrix transpose() const;
    Matrix pseudoInverse() const;
    int getRows() const { return rows; }
    int getCols() const { return cols; }
};
extern "C" {
    __declspec(dllexport) bool isIdentity(double* matrix, int size);
    __declspec(dllexport) bool inverseCramer(double* matrix, int size, double* result);
    __declspec(dllexport) bool inverseGauss(double* matrix, int size, double* result);
    __declspec(dllexport) void pseudoInverse(double* matrix, int rows, int cols, double* result);
    __declspec(dllexport) double determinant(double* matrix, int size);
}
#endif // MATRIX_H