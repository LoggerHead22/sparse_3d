#pragma once

#include <cmath>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>
#include <cstring>
#include <pthread.h>
using namespace std;


#define EPS 1e-14

#define LOG(...) std::cout<< #__VA_ARGS__<< " = " <<(__VA_ARGS__)<< "\n"
#define LN       std::cout << "\n";

static pthread_barrier_t barrier;


void coord_trans( double x, double y, double &x_ , double &y_);
void matr_mult_vector(double *A, int *I, int M, double *x, double *b, int k, int p);
int get_k(int nx , int ny , int i , int j);
void get_ij(int nx, int ny, int k, int &i , int &j);
int get_num_offdiag(int nx, int ny, int k);
int get_non_zeros(int nx, int ny);
int get_offdiag_elem( int nx , int ny , int k , double *a_diag , double *a , int *I);
int allocate_MSR_matrix(int nx, int ny , double *&p_a , double *&p_I);
void build_MSR_matrix(int nx , int ny, double*a, int *I, int p , int k);
void* msl_approx(void *in_arg);
void reduce_sum(int p, int *sum = 0);
void print_vector(double* m, int size);
void print_MSR_matrix(double * A , int * I, int N );








struct Arg{
        int nx, ny ,m, p, thr_ind ,*error , *I;
        double * A,*b , *x, time_thr, fulltime, n_err , *xs , *ys;
		double hx , hy;
};
