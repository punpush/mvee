#include "ellips.h"
#include <random>
#include "eigen.h"

Ellips::Ellips(std::vector<std::vector<double> > &data, double tolerance)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(1, 10);

    for (auto &r : data) {
        for (auto &el : r) {
            el = dis(gen);
        }
    }

    // Dimension of the points
    double d = 2.0;

    size_t N = data.size();

    auto P = transpose(data);

    std::vector<std::vector<double>> Q = merge(P, ones(1, N));

    int count = 1;
    double err = 1;

    std::vector<double> u(N, (1 / static_cast<double>(N)));

     // Khachiyan Algorithm
    while (err > tolerance) {
        // Matrix multiplication:
        // diag(u) : if u is a vector, places the elements of u
        // in the diagonal of an NxN matrix of zeros
        //X = Q*diag(u)*Q'; // Q' - transpose of Q
        auto X = multiply(multiply(Q,diag(u)), transpose(Q));

        // inv(X) returns the matrix inverse of X
        // diag(M) when M is a matrix returns the diagonal vector of M
        //M = diag(Q' * inv(X) * Q); // Q' - transpose of Q
        auto M = diag(multiply(multiply(transpose(Q), inv(X)), Q));

        //Find the value and location of the maximum element in the vector M
        double maximum = max(M);
        int j = findMaximumValueLocation(M, maximum);


        // Calculate the step size for the ascent
        double step_size = (maximum - d - 1) / ((d+1) * (maximum - 1));


        // Calculate the new_u:
        // Take the vector u, and multiply all the elements in it by (1-step_size)
        auto new_u = multiply((1 - step_size), u);


        // Increment the jth element of new_u by step_size
        new_u[j] = new_u[j] + step_size;


        // Calculate error by taking finding the square root of the SSD
        // between new_u and u
        err = std::sqrt(ssd(new_u, u));


        // Increment count and replace u
        count = count + 1;
        u = new_u;
    }

    auto c = multiply(P, u);
    m_center = transpose(c)[0];
    auto U = diag(u);

    auto pup = multiply(multiply(P, U) , transpose(P));
    auto pupu = multiply((multiply(P, u)), transpose(multiply(P, u)));
    auto pup_pupu = minus(pup, pupu);
    auto a = inv(pup_pupu);
    m_data = multiply((1 / d), inv(pup_pupu));

    auto eig = Eigen(inv(m_data));
    auto Ve = eig.getV(); //eigenvalues
    auto De = eig.getD(); //right eigenvectors
    reorderEigenVectors(De);
    reorderEigenValues(Ve);

    auto v = sqrt(diag(De));
    smalSemi = max(v);
    int Ie = findMaximumValueLocation(v, smalSemi);

    std::vector<double> veig(Ve.size());
    for (int i=0; i<veig.size(); i++){
        veig[i] = Ve[Ie][i];
    }

    thu = std::atan2(veig[1], veig[0]);
    std::vector vec{0,1};
    bigSemi = v[setDiff({0, 1}, Ie)];
}

void Ellips::gaussian(std::vector<std::vector<double> > &data, std::vector<int> &index)
{
    int n = index.size();
    std::vector<double> c(n);

    for (int i = 0; i < n; ++i)
        index[i] = i;

    for (int i = 0; i < n; ++i) {
        double c1 = 0;
        for (int j = 0; j < n; j++) {
            double c0 = std::fabs(data[i][j]);
            if (c0 > c1)
                c1 = c0;
        }
        c[i] = c1;
    }

    int k = 0;

    for (int j = 0; j < n - 1; ++j) {
        double pi1 = 0;
        for (int i = j; i < n; ++i) {
            double pi0 = std::fabs(data[index[i]][j]);
            pi0 /= c[index[i]];
            if (pi0 > pi1) {
                pi1 = pi0;
                k = i;
            }
        }

        int itmp = index[j];
        index[j] = index[k];
        index[k] = itmp;

        for (int i = j + 1; i < n; ++i) {
            double pj = data[index[i]][j] / data[index[j]][j];
            data[index[i]][j] = pj;

            for (int l = j + 1; l < n; ++l)
                data[index[i]][l] -= pj * data[index[j]][l];
        }
    }
}

std::vector<std::vector<double> > Ellips::inv2(std::vector<std::vector<double> > &data)
{
    int n = data.size();
    std::vector<std::vector<double>> x(n, std::vector<double>(n));
    std::vector<std::vector<double>> b(n, std::vector<double>(n));
    std::vector<int> index(n);

    for (int i = 0; i < n; ++i)
        b[i][i] = 1;

    gaussian(data, index);

    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                b[index[j]][k] -= data[index[j]][i]*b[index[i]][k];
            }
        }
    }

    for (int i = 0; i < n; ++i){
        x[n - 1][i] = b[index[n - 1]][i] / data[index[n - 1]][n - 1];
        for (int j = n - 2; j >= 0; --j){
            x[j][i] = b[index[j]][i];
            for (int k = j + 1; k < n; ++k){
                x[j][i] -= data[index[j]][k] * x[k][i];
            }
            x[j][i] /= data[index[j]][j];
        }
    }

    return x;
}

std::vector<std::vector<double> > Ellips::inv1(const std::vector<std::vector<double> > &data)
{
    std::vector<std::vector<double>> inverse(data.size(), std::vector<double>(data.size()));

    // minors and cofactors
    for (int i = 0; i < data.size(); i++)
        for (int j = 0; j < data[i].size(); j++)
            inverse[i][j] = std::pow(-1, i + j) * determinant(minor(data, i, j));

    // adjugate and determinant
    double det = 1.0 / determinant(data);
    for (int i = 0; i < inverse.size(); i++) {
        for (int j = 0; j <= i; j++) {
            double temp = inverse[i][j];
            inverse[i][j] = inverse[j][i] * det;
            inverse[j][i] = temp * det;
        }
    }

    return inverse;
}

std::vector<std::vector<double> > Ellips::inv(std::vector<std::vector<double> > &data)
{
    // try{
    //     return inv1(data);
    // }
    // catch(int e){
    //     try{
    //         return inv2(data);
    //     }
    //     catch(int ex){
    //         (void)(ex);
    //     }
    // }

    auto res = inv1(data);

    bool hasError {false};

    for (auto &i : res) {
        for (auto j : i) {
            if (j != j) {
                hasError = true;
                break;
            }
        }
    }

    if (!hasError)
        return res;

    res.clear();

    res = inv2(data);

    hasError = false;

    for (auto &i : res) {
        for (auto j : i) {
            if (j != j) {
                hasError = true;
                break;
            }
        }
    }

    if (!hasError)
        return res;

    return {};
}

std::vector<std::vector<double>> Ellips::transpose(const std::vector<std::vector<double> > &data)
{
    if (data.empty())
        return {};

    std::vector<std::vector<double>> res {data[0].size(), std::vector<double>(data.size())};

    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[i].size(); j++) {
            res[j][i] = data[i][j];
        }
    }

    return res;
}

void Ellips::rotateMatrix(std::vector<std::vector<double>> &data)
{
    int N = data[0].size();

    for (int x = 0; x < N / 2; ++x) {
        for (int y = x; y < N - x - 1; ++y) {
            // store current cell in temp variable
            double temp = data[x][y];

            // move values from right to top
            data[x][y] = data[y][N-1-x];

            // move values from bottom to right
            data[y][N-1-x] = data[N-1-x][N-1-y];

            // move values from left to bottom
            data[N-1-x][N-1-y] = data[N-1-y][x];

            // assign temp to left
            data[N-1-y][x] = temp;
        }
    }
}

int Ellips::setDiff(const std::vector<double> &data, int x)
{
    for (auto i : data) {
        if (i != x)
            return i;
    }

    return 0; //???
}

std::vector<double> Ellips::sin(const std::vector<double> &data)
{
    std::vector<double> res((data.size()), 0);
    for (int i = 0; i < res.size(); ++i) {
        res[i] = std::sin(data[i]);
    }
    return res;
}

std::vector<double> Ellips::cos(const std::vector<double> &data)
{
    std::vector<double> res((data.size()), 0);
    for (int i = 0; i < res.size(); ++i) {
        res[i] = std::sin(data[i]);
    }
    return res;
}

std::vector<double> Ellips::sqrt(const std::vector<double> &data)
{
    std::vector<double> res((data.size()), 0);
    for (int i = 0; i < res.size(); ++i) {
        res[i] = std::sqrt(data[i]);
    }
    return res;
}

std::vector<double> Ellips::multiply(double n, const std::vector<double> &data)
{
    std::vector<double> res(data.size());
    for (int i = 0; i < data.size(); ++i){
        res[i] = data[i] * n;
    }
    return res;
}

std::vector<double> Ellips::diag(const std::vector<std::vector<double> > &data)
{
    std::vector<double> diag(data.size());

    for (int i = 0; i < data.size(); ++i)
        diag[i] = data[i][i];

    return diag;
}

std::vector<std::vector<double> > Ellips::diag(const std::vector<double> &data)
{
    std::vector<std::vector<double>> diag(data.size());
    for (int i = 0; i < data.size(); ++i) {
        std::vector<double> row(data.size());
        for (int j = 0; j < data.size(); ++j) {
            if (i == j)
                row[j] = data[i];
            else
                row[j] = 0;
        }
        diag[i] = row;
    }

    return diag;
}

std::vector<std::vector<double> > Ellips::minus(const std::vector<std::vector<double> > &A, const std::vector<std::vector<double> > &B)
{
    std::vector<std::vector<double>> res = {A.size(), std::vector<double>()};
    for (int i = 0; i < res.size(); ++i) {
        std::vector<double> row((A[i].size()));
        for (int j = 0; j < row.size(); ++j) {
            row[j] = A[i][j] - B[i][j];
        }
        res[i] = row;
    }

    return res;
}

std::vector<std::vector<double> > Ellips::multiply(const std::vector<std::vector<double> > &A, const std::vector<std::vector<double> > &B)
{
    int m1ColLength = A[0].size(); // m1 columns length
    int m2RowLength = B.size();    // m2 rows length
    if(m1ColLength != m2RowLength) return {}; // matrix multiplication is not possible

    int mRRowLength = A.size();    // m result rows length
    int mRColLength = B[0].size(); // m result columns length

    std::vector<std::vector<double>> mResult(mRRowLength, std::vector<double>(mRColLength));
    for(int i = 0; i < mRRowLength; i++) {         // rows from m1
        for(int j = 0; j < mRColLength; j++) {     // columns from m2
            for(int k = 0; k < m1ColLength; k++) { // columns from m1
                mResult[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return mResult;
}

std::vector<std::vector<double> > Ellips::multiply(const std::vector<std::vector<double> > &A, const std::vector<double> &B)
{
    std::vector<std::vector<double>> m2(B.size());
    for (int i = 0; i < m2.size(); ++i){
        std::vector<double> row(1);
        row[0] = B[i];
        m2[i] = row;
    }

    return multiply(A, m2);
}

std::vector<std::vector<double> > Ellips::multiply(double n, const std::vector<std::vector<double> > &data)
{
    std::vector<std::vector<double>> res(data.size());
    for (int i = 0; i < data.size(); ++i){
        std::vector<double> row = data[i];
        std::vector<double> r(row.size());
        for (int j = 0; j < row.size(); ++j){
            r[j] = row[j] * n;
        }
        res[i] = r;
    }
    return res;
}

std::vector<std::vector<double> > Ellips::merge(const std::vector<std::vector<double> > &A, const std::vector<std::vector<double> > &B)
{
    int x = 0;
    std::vector<std::vector<double> > res(A.size() + B.size());
    for (int i = 0; i < A.size(); ++i){
        res[x] = A[i];
        x++;
    }
    for (int i = 0; i < B.size(); ++i){
        res[x] = B[i];
        x++;
    }
    return res;
}

std::vector<std::vector<double> > Ellips::ones(int rows, int cols)
{
    std::vector<std::vector<double>> arr(rows);
    for (int i = 0; i < arr.size(); ++i){
        std::vector<double> row(cols);
        for (int j = 0; j < row.size(); ++j){
            row[j] = 1;
        }
        arr[i] = row;
    }
    return arr;
}

double Ellips::ssd(const std::vector<double> &A, const std::vector<double> &B)
{
    double ssd = 0;
    for (int i = 0; i < A.size(); ++i){
        ssd += std::pow(A[i] - B[i], 2);
    }
    return ssd;
}

double Ellips::max(const std::vector<double> &data)
{
    double max = data[0];
    for (double d : data){
        max = std::fmax(d, max);
    }
    return max;
}

int Ellips::findMaximumValueLocation(const std::vector<double> &data, double max)
{
    for (int i = 0; i < data.size(); i++){
        if (data[i] == max) return i;
    }
    return 0;
}

std::vector<double> Ellips::linspace(double min, double max, int points)
{
    std::vector<double> d(points);

    for (int i = 0; i < points; i++){
        d[i] = min + i * (max - min) / (points - 1);
    }
    return d;
}

void Ellips::reorderEigenValues(std::vector<std::vector<double> > &data)
{
    rotateMatrix(data);
    for (int i = 0; i < data.size(); ++i){
        for (int j = 0; j < data[i].size(); ++j){
            data[i][j] = -data[i][j];
        }
    }
}

void Ellips::reorderEigenVectors(std::vector<std::vector<double> > &data)
{
    rotateMatrix(data);
    rotateMatrix(data);
}

std::vector<std::vector<double> > Ellips::getBoundingCoordinates(int numPoints)
{
    return {};
    // //tq=linspace(-pi,pi,50);
    // auto tq = linspace(-M_PI, M_PI, numPoints);


    // //U=[cos(thu) -sin(thu);sin(thu) cos(thu)]*[l1*cos(tq);l2*sin(tq)];
    // auto U = multiply(
    //     std::vector<std::vector<double>>(
    //         createVector(std::cos(thu), -std::sin(thu)),
    //         createVector(std::sin(thu), std::cos(thu))
    //         ),
    //     new double[][]{
    //         multiply(l1, cos(tq)),
    //         multiply(l2, sin(tq))
    //     }
    //     );
    // //System.out.println(toString(transpose(U)));



    // double[][] coords = transpose(U);
    // for (int i=0; i<coords.length; i++){
    //     double x = coords[i][0] + center[0];
    //     double y = coords[i][1] + center[1];

    //     coords[i][0] = x;
    //     coords[i][1] = y;
    // }

    // return coords;
}

std::vector<std::vector<double> > &Ellips::matrix()
{
    return m_data;
}

std::vector<double> &Ellips::center()
{
    return m_center;
}

std::vector<std::vector<double> > Ellips::minor(const std::vector<std::vector<double> > &data, int row, int column)
{
    std::vector<std::vector<double>> minor(data.size() - 1, std::vector<double>(data.size() - 1));

    for (int i = 0; i < data.size(); i++)
        for (int j = 0; i != row && j < data[i].size(); j++)
            if (j != column)
                minor[i < row ? i : i - 1][j < column ? j : j - 1] = data[i][j];

    return minor;
}

double Ellips::determinant(const std::vector<std::vector<double> > &data)
{
    if (data.empty())
        return std::numeric_limits<double>::quiet_NaN();

    if (data.size() != data[0].size())
        return std::numeric_limits<double>::quiet_NaN();

    if (data.size() == 2)
        return data[0][0] * data[1][1] - data[0][1] * data[1][0];

    double det = 0;
    for (int i = 0; i < data[0].size(); i++)
        det += std::pow(-1, i) * data[0][i]
               * determinant(minor(data, 0, i));
    return det;
}
