#pragma once

#define _USE_MATH_DEFINES
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


struct parral{
public:
	double l1 = 1 , l2 = 1, alpha = 45./(180./M_PI) , k = 0.2 , hx, hy;
	int nx = 10 , ny = 10;
	double *xs = nullptr, *ys = nullptr;
	double l1_new , l2_new = l2;
	parral(){};	
	parral(double l1_ , double l2_ , double alpha_ ,  double k_ , int nx_ , int ny_){
		l1 = l1_; l2 = l2_;alpha = alpha_ /(180./M_PI);
		k = k_; nx = nx_; ny = ny_;
		l1_new = sin(alpha)*l1; l2_new = l2;
		hx = l2_new / nx; hy = l1_new / ny;
		xs = new double[nx + 1];
		ys = new double[ny + 1];
		compute_grid();
		//cout<<"HH"<<hx<<" "<<hy<<endl;
	}
	
	void compute_grid(){
		for(int i = 0 ; i<=nx; i++){
			xs[i] = hx*i;
		}
		for(int i = 0 ; i<=ny; i++){
			ys[i] = hy*i;
		}
		
	}
	
	
	double f_par(double x, double y , double (*f)(double, double) ){
		return f(x + cos(alpha) / sin(alpha)* y , y);
	}
	
	~parral(){
		delete[] xs;
		delete[] ys;
		xs = nullptr;
		ys = nullptr;
	};
	
	
};


void coord_trans( double x, double y, double &x_ , double &y_);
void matr_mult_vector(double *A, int *I, int M, double *x, double *b, int k, int p);
int get_k(int nx , int ny , int i , int j);
void get_ij(int nx, int ny, int k, int &i , int &j);
int get_num_offdiag(int nx, int ny, int k);
int get_non_zeros(int nx, int ny);
int get_offdiag_elem( int nx , int ny , int k , double *a_diag , double *a , int *I);
int allocate_MSR_matrix(int nx, int ny , double *&p_a , double *&p_I);
int get_offdiag_elem( int nx , int ny , int k , double *a_diag , double *a , int *I , double * b ,double hx, double hy , parral &par , double (*f) (double,double));
void build_MSR_matrix(int nx , int ny, double*a, int *I, double *b, int p , int k , parral &par , double (*f) (double,double));
void* msl_approx(void *in_arg);
void reduce_sum(int p, int *sum = 0);
void print_vector(double* m, int size);
void print_MSR_matrix(double * A , int * I, int N );











struct Arg{
        int nx, ny ,m, p, thr_ind ,*error , **I;
        double **A,**b , *x, time_thr, fulltime, n_err , *xs , *ys;
		double hx , hy;
		double (*f)(double,double);
		parral* par;
};



