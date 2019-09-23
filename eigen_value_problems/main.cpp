#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
//using namespace arma;


double lambda_analytical(int N, double j,double d,double a){
    double lambda = d + 2*a* std::cos(j*M_PI/(N+1));
    return lambda;
}


double jacobis_method(arma::mat A,int n){
    double tol = 10e-6;
    double max = 2*tol;
    int k,l,iter,iter_max;
    double c, s, tou, t_1, t_2, t,a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    iter = 0;
    iter_max = n*n*n;
    // Setting up the eigenvector matrix
    arma::Mat<double> R = arma::mat(n,n); R.eye();         // R.eye() automatically sets identity matrix
    while (max > tol && iter < iter_max){
        iter++;
        max = tol;
        for (int i = 0; i < n; i++){           // finds max value of off-diagonal elements assuming symmetric matrix
            for (int j = i+1; j < n; j++){
                double aij = abs(A(i,j));
                if (aij > max){
                    max = aij; k = i; l = j;
                    }
        }
        }
        if (A(k,l) != 0){
            tou = (A(l,l) - A(k,k))/(2*A(k,l));
            t_1 = -tou + sqrt(1+pow(tou,2));
            t_2 = -tou - sqrt(1+pow(tou,2));
            if (t_1 < t_2){ //chose smallest value of t
                t = t_1;
        }
            else{
                t = t_2;
        }
        c = 1/sqrt(1+pow(t,2));
        s = t*c;}
        else{
            c = 1;
            s = 0;
        }


        //Dette er i stor grad en kopi av morten sin kode
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
    }       // end of while loop

    arma::Col<double> e_val = arma::vec(n);
    for (int i=0; i<n; i++){
        e_val(i) = A(i,i);
    }

    std::cout << "Iterations out of max: " << iter << ":" << iter_max << endl << e_val << endl;
    //std::cout << "Iterations out of max: " << iter << ":" << iter_max << endl << R << endl << A << endl;
    //slutt paa plagiat



    //test 1: The largest non-diagonal value is correct (precise) for small matrices (max 5x5)
    //int n_test = 5;
    //double max2 = 1e-10;
    //for (int i=0; i<n_test; i++){
    //    for (int j=0; j<n_test; j++){
    //        if (A(i,j)>max){
    //            max2 = A(i,j);
    //        }
    //    }
    //}           // end of test 1


    //test 2: sum over all non-diagonal elements is less then a small error
    //double err = 10e-6;
    //double sum;
    //for (int i=0; i<n; i++){
    //    for (int j=i+1; j<n; j++){
    //        sum = sum + 2*abs(A(i,j));
    //}
    //}
    //if (sum > err){
    //   std::cout  << "The sum over the non-diagonal elements is greater than the biggest possible error!" << endl;
    //}
    //else{
    //    std::cout  << "The non-diagonal elements are all good!" << endl;    // maybe delete this
    //}           // end of test 2


    return 0;
    }






int main(){
    int n;
    std::cin >> n;
    //unsigned int n = 5;
    int d = 5;
    int a = 3;
    //double tol = 10e-8;
    arma::Mat <double> A = arma::mat(n,n); A.zeros();


    // Setting up the tridiagonal matrix A
    for(unsigned int i=0; i<n-1; i++){
        A(i,i) = d;
        A(i+1,i) = a;
        A(i,i+1) = a;
      }
    A(n-1,n-1) = d;


    // armadillo and analytical calculation of eigenvalues
    arma::Col<double> eig_val;
    arma::eig_sym(eig_val, A);                      // armadillo calculations
    arma::Col<double> lambs(n);
    for(unsigned int j=0; j<n; j++){                // analytical calculations
        lambs(j) = lambda_analytical(n, j+1, d, a);
    }
    //std::cout << "Comparing the armadillos version of eigenvalues with the analytical solution:" << std::endl
    //     << "Armadillo:" << std::endl << eig_val << std::endl << "Analytic:" << std::endl << lambs << std::endl;


    jacobis_method(A,n);


    std::cout << "Analytic eigenvalues: " << endl << eig_val << endl;
    return 0;
}
