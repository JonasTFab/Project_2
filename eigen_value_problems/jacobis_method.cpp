#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

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

    return A;
}
//cout << A << endl;

// the offdiag function
int non_diag(mat A,mat R, int k, int l, int n){
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
    double t;
    if (t_1 < t_2){ //chose smalles value of t
      t = t_1;
    }
    else{
      t = t_2;
    }
    double c = 1/sqrt(1+pow(t,2));
    double s = t*c;

    //Dette er i stor grad en kopi av morten sin kode
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    //change the matrix elements index k,l
    A(k,k) = pow(c,2) * a_kk - 2*c*s*A(k,l) + pow(s,2)*a_ll;
    A(l,l) = pow(s,2) * a_kk + 2*c*s*A(k,l) + pow(c,2)*a_ll;
    A(k,l) = 0.0;
    A(l,k) = 0.0;
    for(int i = 0; i < n; i++){
      if (i != k && i != l){
        a_ik = A(i,k);
        a_il = A(i,l);
        A(i,k) = c*a_ik - s*a_il;
        A(k,i) = A(i,k);
        A(i,l) = c*a_il + s*a_ik;
        A(l,i) = A(i,l);
      }
    //compute eigenvectors
    r_ik = R(i,k);
    r_il = R(i,l);
    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;
    }
    //slutt paa plagiat
    return 0;
    }
