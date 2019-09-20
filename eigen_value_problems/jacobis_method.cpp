#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

// the offdiag function, using Armadillo
void non_diag(mat A, int &k, int &l, int n){
      double tol = 10e-8;
      for (int i = 0; i < n; i++){
        for ( int j = i+1; j < n; j++){
            double aij = abs(A(i,j));
            if ( aij > tol){
              tol = aij; k = i; l = j;
    }
    }
    }
    double tou = (A(l,l) - A(k,k))/(2*tol);
    double t_1 = -tou + sqrt(1+pow(tou,2));
    double t_2 = -tou - sqrt(1+pow(tou,2));
    if (t_1 < t_2){
      t = t_1
    }
    else{
      t = t_2
    }

    //cout << tou << pow(tou,2) << endl;
    }


int main(){
    unsigned int n = 5;
    int diag = 5;
    int semi_diag = 3;
    //double tol = 10e-8;
    arma::Mat <double> A = arma::mat(n,n); A.zeros();

    for(unsigned int i=0; i<n-1; i++){
        A(i,i) = diag;
        A(i+1,i) = semi_diag;
        A(i,i+1) = semi_diag;
      }

    A(n-1,n-1) = diag;

    cout << A << endl;
    //non_diag(A,float k,double l,n);

    return 0;




}
