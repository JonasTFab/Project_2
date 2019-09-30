# Project 2: Eigenvalue problem

In this project we wanted to code in C++ for the first time. It has been a lot of educational experience when it comes to understanding a new programming language and will probably use it more in the future.

The C++ code file can be found in the folder named eigen_value_problems named main.cpp. We have now reduced the number of several codes to just one in unlike what we did in project 1. The code includes solving the eigenvalue problem for both an Töplitz matrix and also a tridiagonal matrix where the diagonal elements are changing.

By unslashing line 242 in the file, you will run the code that solves the Töplitz matrix. 

By unslashing line 254, you will activate a function that runs the algorithm throught increasing values of the size of the matrix N. It will then store the number of iteration, and all values of N, in a text file so we can plot them using python. The text file is located in the folder build-eigen_value_problems-Desktop_Qt_5_13_1_MinGW_64_bit-Debug. You will also find the python code which plots these values.

By unslashing line 258, you will run the algorithm with an matrix that corresponds to the Schrödinger equation.
