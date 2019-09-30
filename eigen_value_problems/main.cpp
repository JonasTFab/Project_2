#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <ctime>
#include <chrono>

using namespace std;
//using namespace arma;


double lambda_analytical(int N, double j,double d,double a){
    double lambda = d + 2*a* std::cos(j*M_PI/(N+1));
    return lambda;
}


int jacobis_method(arma::mat A, arma::Col<double> d, arma::Col<double> a, int n){
    double tol = 10e-8;
    double max = 2*tol;                // just to make max a little bigger than the tolerance before running the while loop
    int k,l,iter,iter_max;
    double c, s, tou, t_1, t_2, t,a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    iter = 0;
    iter_max = n*n*n;

    // Setting up the eigenvector matrix (R.eye() automatically sets identity matrix)
    arma::Mat<double> R = arma::mat(n,n).eye();

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

    // Puts the diagonal elements (eigenvalues) in matrix A into a column vector
    arma::Col<double> e_val = arma::vec(n);
    arma::Col<double> sorted_e_val, e_val_analytic(n),sorted_e_val_analytic, diff;
    for (int i=0; i<n; i++){
        e_val(i) = A(i,i);
    }


    // test 2a: comparing the analytical and numerical eigenvalues (works only for a Töplitz matrix)
    if (d(0)==d(1) && d(0)==d(n-1) && a(0)==a(1) && a(0)==a(n-1)){
        sorted_e_val = arma::sort(e_val);
        for(int i=0; i<n; i++){                // analytical calculations
            e_val_analytic(i) = lambda_analytical(n, i+1, d(0), a(0));
        }
        sorted_e_val_analytic = arma::sort(e_val_analytic);
        for (int i=0; i<n; i++){
            double diff = sorted_e_val(i) - sorted_e_val_analytic(i);
            if (abs(diff) > tol){
                std::cout << "The difference between the eigenvalues are bigger than the tolerance!" << endl <<
                         "Index: " << i << "     Analytic eigenvalue: " << sorted_e_val_analytic(i) <<
                         "     Numerical eigenvalue: " << sorted_e_val(i) << "\n\n" ;
        }
        }
    } // end of test 2a
    // test 2b: comparing the analytical and numerical eigenvalues (works for Schrödinger equation)
    else if (d(0)!=d(1) && d(0)!=d(n-1) && a(0)==a(1) && a(0)==a(n-1)){
      double  tol_4 = 10e-4; //4 leading digits
        sorted_e_val = arma::sort(e_val);
        for (int i=0; i<n; i++){
            e_val_analytic(i) = 3 + i*4;       // assuming there exists a pattern for the eigenvalues
        }
        for (int i=0; i<4; i++){
            double diff = sorted_e_val(i) - e_val_analytic(i);
            if (abs(diff) > tol_4){
                std::cout << "The difference between the eigenvalues are bigger than the tolerance!" << endl <<
                         "Index: " << i << "     Analytic eigenvalue: " << e_val_analytic(i) <<
                         "     Numerical eigenvalue: " << sorted_e_val(i) << "\n\n" ;
        }
        }
    } // end of test 2b

    else {
        arma::cout << "Missing an analytic solution to compare with!" << endl;
        return 0;
    }
    //OUTPUT!!!!!!!!!!!!!!!!
    //std::cout << "The calculated eigenvalues is: " << endl << e_val << "\n\n";
    //std::cout << "The analytical eigenvalues is: " << endl << e_val_analytic << "\n\n";
    //std::cout << "Number of iterations: " << iter << "\n\n" ;
    //std::cout << "Theoretical number of rotations between 3n^2-5n^2: " << 3*pow(n,2)<< "-" << 5*pow(n,2)<< "\n\n";


    return iter;
    }




// This function creates a file which stores the number of iterations
// with respect to the dimension of the matrix using Jacobi's method
double iterations_per_dimension(arma::mat A, arma::Col<double> d, arma::Col<double> a, int n){
    arma::Col<double> dimensions = arma::vec(int (n/5));
    arma::Col<double> iterations = arma::vec(int (n/5));
    std::ofstream data;
    data.open("iter_per_dim.txt");
    for (int i=0; i<int (n/5); i++){
        dimensions(i) = 5*(i+1);
        iterations(i) = jacobis_method(A,d,a, dimensions(i));
        std::cout << "Dimension: " << dimensions(i) << "     Iterations: " << iterations(i) << endl;
        data << dimensions(i) << " " << iterations(i) << endl;
    }
    data.close();

    return 0;
}



// Solving the eigenvalues for schödinger equation using Jacobi's method
double schrodinger_solver(int N){
    // some initial conditions
    double rho_min, rho_N, h;
    rho_min = 0;
    //rho_N = 3.55;
    std::cout << "Chose rho_N: ";
    std::cin >> rho_N;
    h = (rho_N-rho_min)/(N-1);
    //cout << "h: " << h << endl;
    arma::Col<double> rho = arma::linspace(rho_min,rho_N,N);
    arma::Col<double> V = pow(rho,2);
    arma::Col<double> d = (2/pow(h,2))+V;
    arma::Col<double> e = arma::vec(N).ones()*((-1)/pow(h,2));
    //cout << "rho: " << endl << rho << endl;

    // setting up the Schrödinger matrix
    arma::Mat<double> S = arma::mat(N,N);
    for (int i=0; i < N-1; i++){
        S(i,i) = d(i);
        S(i+1,i) = e(i);
        S(i,i+1)= e(i);
    }
    S(N-1,N-1) = d(N-1);
    jacobis_method(S,d,e,N);



    return 0;


}





int main(){
    int n;
    std::cout << "Size of matrix: ";
    std::cin >> n;
    //int n = 20;
    arma::Col<double> d = arma::vec(n).ones()*7;
    arma::Col<double> a = arma::vec(n).ones()*-3;
    arma::Mat <double> A = arma::mat(n,n).zeros();


    // Setting up the tridiagonal matrix A
    for(int i=0; i<n-1; i++){
        A(i,i) = d(i);
        A(i+1,i) = a(i);
        A(i,i+1) = a(i);
      }
    A(n-1,n-1) = d(n-1);


    // armadillo and analytical calculation of eigenvalues
    arma::Col<double> eig_val;
    auto start = std::chrono::high_resolution_clock::now();
    arma::eig_sym(eig_val, A);                      // armadillo calculations
    auto finish = std::chrono::high_resolution_clock::now();

    arma::Col<double> lambs(n);
    for(int j=0; j<n; j++){                         // analytical calculations
        lambs(j) = lambda_analytical(n, j+1, d(0), a(0));
    }
    //std::cout << "Comparing the armadillos version of eigenvalues with the analytical solution:" << std::endl
    //     << "Armadillo:" << std::endl << eig_val << std::endl << "Analytic:" << std::endl << lambs << std::endl;


    // running the Jacobi's rotational algorithm to find eigenvalues
    auto start2 = std::chrono::high_resolution_clock::now();
    jacobis_method(A,d,a,n);
    auto finish2 = std::chrono::high_resolution_clock::now();

    // Calculating the time taken using armadillo's solver and our algorithm
    std::chrono::duration<double> elapsed_arma = finish - start;
    std::chrono::duration<double> elapsed_jacobi = finish2 - start2;


    std::cout << "Time in seconds using armadillo and Jacobi's method with n = " << n << endl
              << "Armadillo time: " << elapsed_arma.count() << endl << "Jacobi time: " << elapsed_jacobi.count() << endl;

    // creating a file for plotting in python
    //iterations_per_dimension(A,d,a,n);

    // solves Schrödinger equation
    time_t time_SE_0 = time(NULL);
    schrodinger_solver(n);
    time_t time_SE_1 = time(NULL);
    t_SE = time_SE_1 - time_SE_0;
    std::cout << "Time in seconds solving the Schrödinger equations with n = " << n << ": " << t_SE << "seconds" << endl;



    return 0;
}
