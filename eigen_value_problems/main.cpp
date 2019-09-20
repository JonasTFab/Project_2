#include <iostream>
#include <armadillo>
#include <cmath>
using namespace std;



int main(){
    std::cout << M_PI << endl;
    unsigned int n = 5;
    int diag = 5;
    int semi_diag = 3;
    arma::Mat <double> A = arma::mat(n,n); A.zeros();

    for(unsigned int i=0; i<n-1; i++){
        A(i,i) = diag;
        A(i+1,i) = semi_diag;
        A(i,i+1) = semi_diag;

    }

    A(n-1,n-1) = diag;
    // Finds eigenvalues of matrix A stored in vector "eigval".
    // Then the eigenvectors in a matrix P. At last creates a
    // diagonal matrix D with the eigenvalues at the diagonal.
    arma::Col<double> eigval;
    arma::Mat<double> P;
    eig_sym(eigval, P, A);                      // calculates eigenvalues and eigenvectors using armadillo
    arma::Mat<double> P_t = P.t();              // transpose the matrix P
    arma::Mat<double> D = arma::mat(n,n); D.zeros();

    for(unsigned int i=0; i<n; i++){
        D(i,i) = eigval(i);
    }

    arma::Mat<double> A_2 = P * D * P_t;

    //cout << "Transpose matrix P:\n" << P_t << endl << "Eigenvalues:\n"<< D << endl << P << endl << A_2 << endl;

    return 0;
}
