#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <ctime>


using namespace std;
//using namespace arma;


double lambda_analytical(int N, double j,double d,double a){
    double lambda = d + 2*a* std::cos(j*M_PI/(N+1));
    return lambda;
}


int jacobis_method(arma::mat A,double d, double a, int n){
    double tol = 10e-8;
    double max = 2*tol;                // just to make max a little bigger than the tolerance before running the while loop
    int k,l,iter,iter_max;
    double c, s, tou, t_1, t_2, t,a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    iter = 0;
    iter_max = n*n*n;

    // Setting up the eigenvector matrix (R.eye() automatically sets identity matrix)
    arma::Mat<double> R = arma::mat(n,n); R.eye();

    // finds max value of off-diagonal elements assuming symmetric matrix
    while (max > tol && iter < iter_max){
        iter++;
        max = tol;
        for (int i = 0; i < n; i++){
            for (int j = i+1; j < n; j++){
                double aij = abs(A(i,j));
                if (aij > max){
                    max = aij; k = i; l = j;
                    }
        }
        }
        // defines tou, tan, cos and sin
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
        // test 1: check if orthogonality is preserved
        arma::Mat<double> A_tr = A.t(); //transpose A
        for (int i=0; i<n; i++){
            double dots = arma::norm_dot(A_tr.row(i),A.col(i));
            if (abs(dots - 1) > tol){
              std::cout << "Orthogonality is not preserved!!!" << endl;
            }
        }       // end of test 1
    }       // end of while loop

    arma::Col<double> e_val = arma::vec(n);
    for (int i=0; i<n; i++){
        e_val(i) = A(i,i);
    }


    // test 2: comparing the analytical and numerical eigenvalues
    arma::Col<double> sorted_e_val, e_val_analytic(n), sorted_e_val_analytic, diff;
    sorted_e_val = arma::sort(e_val);
    for(int i=0; i<n; i++){                // analytical calculations
        e_val_analytic(i) = lambda_analytical(n, i+1, d, a);
    }
    sorted_e_val_analytic = arma::sort(e_val_analytic);
    for (int i=0; i<n; i++){
        double diff = sorted_e_val(i) - sorted_e_val_analytic(i);
        if (diff > tol){
            std::cout << "The difference between the eigenvalues is bigger than the tolerance!" << endl <<
                         "Index: " << i << "     Analytic eigenvalue: " << sorted_e_val_analytic(i) <<
                         "     Numerical eigenvalue: " << sorted_e_val(i) << "\n\n" ;
        }
    } // end of test 2


    //OUTPUT!!!!!!!!!!!!!!!!
    //std::cout << "Number of iterations: " << iter << endl ;
    //std::cout << "Theoretical number of rotations between 3n^2-5n^2: " << 3*pow(n,2)<< "-" << 5*pow(n,2)<< endl;


    return iter;
    }




// This function creates a file which stores the number of iterations
// with respect to the dimension of the matrix using Jacobi's method
double iterations_per_dimension(arma::mat A, double d, double a, int n){
    arma::Col<double> dimensions = arma::vec(int (n/5));
    arma::Col<double> iterations = arma::vec(int (n/5));
    std::ofstream data;
    data.open("iter_per_dim.txt");
    for (int i=0; i<int (n/5); i++){
        dimensions(i) = 5*(i+1);
        iterations(i) = int(jacobis_method(A,d,a, dimensions(i)));
        std::cout << "Dimension: " << dimensions(i) << "     Iterations: " << iterations(i) << endl;
        data << dimensions(i) << " " << iterations(i) << endl;
    }
    data.close();

    return 0;
}





int main(){
    int n;
    std::cout << "Size of matrix: ";
    std::cin >> n;
    //int n = 20;
    int d = 5;
    int a = 3;
    arma::Mat <double> A = arma::mat(n,n); A.zeros();


    // Setting up the tridiagonal matrix A
    for(int i=0; i<n-1; i++){
        A(i,i) = d;
        A(i+1,i) = a;
        A(i,i+1) = a;
      }
    A(n-1,n-1) = d;


    // armadillo and analytical calculation of eigenvalues
    arma::Col<double> eig_val;
    time_t time_arma_0 = time(NULL);
    arma::eig_sym(eig_val, A);                      // armadillo calculations
    time_t time_arma_1 = time(NULL);
    arma::Col<double> lambs(n);
    for(int j=0; j<n; j++){                         // analytical calculations
        lambs(j) = lambda_analytical(n, j+1, d, a);
    }
    //std::cout << "Comparing the armadillos version of eigenvalues with the analytical solution:" << std::endl
    //     << "Armadillo:" << std::endl << eig_val << std::endl << "Analytic:" << std::endl << lambs << std::endl;


    // running the Jacobi's rotational algorithm to find eigenvalues
    time_t time_jaco_0 = time(NULL);
    jacobis_method(A,d,a,n);
    time_t time_jaco_1 = time(NULL);

    double t_arma, t_jaco;
    t_arma = time_arma_1-time_arma_0;
    t_jaco = time_jaco_1-time_jaco_0;
    std::cout << "Time in seconds using armadillo and Jacobi's method with n = " << n << endl
              << "Armadillo time: " << t_arma << endl << "Jacobi time: " << t_jaco << endl;

    // creating a file for plotting in python
    //iterations_per_dimension(A,d,a,n);

    return 0;
}
