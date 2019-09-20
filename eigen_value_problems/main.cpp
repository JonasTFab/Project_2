#include <iostream>
#include <armadillo>
#include <cmath>

double lambda_analytical(int N, double j,double d,double a){
    double lambda = d + 2*a* std::cos(j*M_PI/(N+1));
    return lambda;
}

int main(){
    unsigned int n = 5;
    int d = 5;
    int a = 3;
    arma::Mat<double> A = arma::mat(n,n); A.zeros();


    for(unsigned int i=0; i<n-1; i++){
        A(i,i) = d;
        A(i+1,i) = a;
        A(i,i+1) = a;
    }
    A(n-1,n-1) = d;

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

    arma::Col<double> lambs(n);
    for(unsigned int j=0; j<n; j++){
        lambs(j) = lambda_analytical(n, j+1, d, a);
    }


    //std::cout << P_t << std::endl << D << std::endl << P << std::endl << A_2 << std::endl;
    //std::cout << "Comparing the armadillos version of eigenvalues with the analytical solution:" << std::endl
    //     << "Armadillo:" << std::endl << eigval << std::endl << "Analytic:" << std::endl << lambs << std::endl;


    return 0;
}
