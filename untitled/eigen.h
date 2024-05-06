#ifndef EIGEN_H
#define EIGEN_H

#include <vector>

class Eigen
{
public:
    Eigen(const std::vector<std::vector<double>> &data);

    /** Return the eigenvector matrix
   @return     V
   */
    std::vector<std::vector<double>> getV();

    /** Return the block diagonal eigenvalue matrix
   @return     D
   */
    std::vector<std::vector<double>> getD();
private:
    void tred2();
    void tql2();
    void orthes();
    void hqr2();
    void cdiv(double xr, double xi, double yr, double yi);

    double cdivr, cdivi;
    /** Row and column dimension (square matrix).
   @serial matrix dimension.
   */
    int n;

    /** Symmetry flag.
   @serial internal symmetry flag.
   */
    bool isSymmetric;

    /** Arrays for internal storage of eigenvalues.
   @serial internal storage of eigenvalues.
   */
    std::vector<double> d {};
    std::vector<double> e {};

    /** Working storage for nonsymmetric algorithm.
   @serial working storage for nonsymmetric algorithm.
   */
    std::vector<double> ort {};

    /** Array for internal storage of eigenvectors.
   @serial internal storage of eigenvectors.
   */
    std::vector<std::vector<double>> V {};

    /** Array for internal storage of nonsymmetric Hessenberg form.
   @serial internal storage of nonsymmetric Hessenberg form.
   */
    std::vector<std::vector<double>> H {};
};

#endif // EIGEN_H
