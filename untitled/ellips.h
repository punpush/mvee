#ifndef ELLIPS_H
#define ELLIPS_H

#include <limits>
#include <vector>
#include <algorithm>

class Ellips
{
public:
    Ellips(std::vector<std::vector<double>> &data, double tolerance);

    static void gaussian(std::vector<std::vector<double>> &data, std::vector<int> &index);
    static std::vector<std::vector<double>> inv2(std::vector<std::vector<double> > &data);
    std::vector<std::vector<double>> getBoundingCoordinates(int numPoints);
    std::vector<std::vector<double>> &matrix();
    std::vector<double> &center();

private:
    std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &data);
    void rotateMatrix(std::vector<std::vector<double>> &data);
    int setDiff(const std::vector<double> &data, int x);
    std::vector<double> sin(const std::vector<double> &data);
    std::vector<double> cos(const std::vector<double> &data);
    std::vector<double> sqrt(const std::vector<double> &data);
    std::vector<double> multiply(double n, const std::vector<double> &data);
    std::vector<double> diag(const std::vector<std::vector<double>> &data);
    std::vector<std::vector<double>> diag(const std::vector<double> &data);
    std::vector<std::vector<double>> minus(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B);
    std::vector<std::vector<double>> multiply(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B);
    std::vector<std::vector<double>> multiply(const std::vector<std::vector<double>> &A, const std::vector<double> &B);
    std::vector<std::vector<double>> multiply(double n, const std::vector<std::vector<double>> &data);
    std::vector<std::vector<double>> merge(const std::vector<std::vector<double>> &A, const std::vector<std::vector<double>> &B);
    std::vector<std::vector<double>> ones(int rows, int cols);
    double ssd(const std::vector<double> &A, const std::vector<double> &B);
    double max(const std::vector<double> &data);
    int findMaximumValueLocation(const std::vector<double> &data, double max);
    std::vector<double> linspace(double min, double max, int points);
    void reorderEigenValues(std::vector<std::vector<double>> &data);
    void reorderEigenVectors(std::vector<std::vector<double>> &data);

    static std::vector<std::vector<double>> minor(const std::vector<std::vector<double>> &data, int row, int column);
    static double determinant(const std::vector<std::vector<double>> &data);
    std::vector<std::vector<double>> inv1(const std::vector<std::vector<double> > &data);
    std::vector<std::vector<double>> inv(std::vector<std::vector<double>> &data);

    std::vector<double> m_center {};
    std::vector<std::vector<double>> m_data {};
    double smalSemi {0.0};
    double bigSemi {0.0};
    double thu {std::numeric_limits<double>::min()};
};

#endif // ELLIPS_H
