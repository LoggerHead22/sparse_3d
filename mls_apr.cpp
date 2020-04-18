#include "mls.h"




const double tetta = 0.7;

void coord_trans( double x, double y, double &x_ , double &y_){
	
	x_ = x + tetta *y;
	y_ = y;
	
}


void print_MSR_matrix(double * A , int * I, int N ){
	int row_ind= 0;
	int row_diff;
	for(int i = 0 ; i < N - 1; i++){
		row_ind = I[i];
		row_diff = I[i + 1] - I[i];
		for(int k = 0; k < N; k++){
			if (i==k)  cout<<A[i]<<"\t";
			
			if(k == I[I[i] + I[i+1] - I[i] - row_diff]){
				cout<<A[row_ind + k]<<"\t";
				row_diff--;
			}else{ cout<< 0<<"\t"; }
		
		}
		cout<<endl;
	}
	row_ind = I[N-1];
	row_diff = I[N -1];
	cout<<A[N-1]<<"\t";
	
	for( int k = 0; k < N; k++){
		if(k == I[row_diff]){
			cout<<A[row_ind + row_diff]<<"\t";
			row_diff++;
		}else{
			cout<<0<<"\t";
		}
	
	}
}

void print_matrix (int N, double *a, int *I)
{
  int pos;
  int q;
  int pr_n = N > 16? 16: N;
  for (int i = 0; i < pr_n; i++)
    {
      pos = I[i];
      int raw_size = I[i+1] - I[i];
      for (int j = 0; j < pr_n; j++)
        {
          if (i != j)
            {
              for (q = 0; q < raw_size; q++)
                {
                  if (I[pos + q] == j)
                    {
                      printf ("%.3f \t", a[pos + q]);
                      break;
                    }
                }


              if (q == raw_size)
                {
                  printf ("%0.f \t", 0.);
                }
            }
          else
            {
              printf ("%.3f \t", a[i]);
            }
        }
      printf ("\n");
    }
}


void print_vector(double* m, int size) {
	for (int y = 0; y < ((size < 10) ? size : 10); y++) {
		cout << m[y] << "  ";
	}
	cout << endl;
}



void matr_mult_vector(double *A, int *I, int M, double *x, double *b, int k, int p)
{
    int begin = k*(M/p) + (k < M%p? k : 0);
    int end = begin + M/p + (k < M%p? 1 : 0);
    for(int i = begin; i < end; i++)
    {
//        printf("%lf\n", A[i]*x[i]);
        b[i] = A[i]*x[i];
        for(int j = I[i]; j < I[i+1]; j++)
            b[i] += A[j]*x[I[j]];
    }
}

int get_k(int nx , int ny , int i , int j){
	return i*(ny + 1) + j;
}	


void get_ij(int nx, int ny, int k, int &i , int &j){
	i = k / (ny + 1);
	j = k - i*(ny + 1);
}

int get_num_offdiag(int nx, int ny, int k){
	int i,j;
	get_ij(nx,ny,k,i,j);

	if( i >=1 && i <=nx -1 && j>=1 && j<=ny - 1){
		return 6;
	}
	else if( (i >=1 && i <=nx -1) || (j >=1 && j <=ny-1)){
		return  4;
	}else if( (i ==0 && j ==0) || (i == nx && j ==ny)){
		return  3;
	}else if( (i ==0 && j ==ny) || (i == nx && j == 0)){
		return 2;
	}
	cout<<"get_num_offdiag"<<endl;
	
	abort();
	return -1000;
}

	
int get_non_zeros(int nx, int ny){
	int K = (nx + 1)*(ny + 1);
	
	int nz = 0;
	
	for(int k = 0; k < K; k++){
		nz += get_num_offdiag(nx,ny ,k);
		//cout<<"nz "<<nz<<endl;
	}
	return nz;
}

//total_size = K + 1 + nz 

int get_offdiag_elem( int nx , int ny , int k , double *a_diag , double *a , int *I){
	int i, j; 
	
	get_ij(nx,ny,k,i,j);
	
	
	#define I(x,y) get_k(nx , ny, x , y)
	
	if( i >=1 && i <=nx -1 && j>=1 && j<=ny - 1){
		*a_diag = 0.5;
		 I[0] = I(i,j-1);   a[0] = 1./12;
		 I[1] = I(i,j+1);   a[1] = 1./12;
		 I[2] = I(i-1,j);   a[2] = 1./12;
		 I[3] = I(i+1,j);   a[3] = 1./12;
		 I[4] = I(i-1,j-1); a[4] = 1./12;
		 I[5] = I(i+1,j+1); a[5] = 1./12;
		 return 6;
	}
	else if( (i >=1 && i <=nx -1 && j==0)){
		*a_diag = 0.25;
		I[0] = I(i,j+1);   a[0] = 1./12;
		I[1] = I(i-1,j);   a[1] = 1./24;
		I[2] = I(i+1,j);   a[2] = 1./24;
		I[3] = I(i+1,j+1); a[3] = 1./12;
		return 4;
	}else if( (i >=1 && i <=nx -1 && j==ny)){
		*a_diag = 0.25;
		I[0] = I(i,j-1);   a[0] = 1./12;
		I[1] = I(i-1,j);   a[1] = 1./24;
		I[2] = I(i+1,j);   a[2] = 1./12;
		I[3] = I(i-1,j-1); a[3] = 1./24;
		return 4;
	}else if( (i ==0 && j >=1 && j<=ny-1)){
		*a_diag = 0.25;
		I[0] = I(i,j-1);  a[0] = 1./24;
		I[1] = I(i,j+1);  a[1] = 1./12;
		I[2] = I(i+1,j);  a[2] = 1./12;
		I[3] = I(i+1,j+1);a[3] = 1./24;
		return 4; 
	}else if( (i ==nx && j>=1 && j <= ny - 1)){
		*a_diag = 0.25;
		I[0] = I(i,j-1);   a[0] = 1./24;
		I[1] = I(i,j+1);   a[1] = 1./24;
		I[2] = I(i-1,j);   a[2] = 1./12;
		I[3] = I(i-1,j-1); a[3] = 1./12;
		return 4;
	}else if (i==0 && j==0){
		*a_diag = 1./6;
		I[0] = I(i+1,j);   a[0] = 1./24;
		I[1] = I(i,j+1);   a[1] = 1./24;
		I[2] = I(i+1,j+1); a[2] = 1./12;
		return 3;	
	}else if (i==nx && j==ny){
		*a_diag = 1./6;
		I[0] = I(i,j-1);   a[0] = 1./24; 
		I[1] = I(i-1,j);   a[1] = 1./24;
		I[2] = I(i-1,j-1); a[2] = 1./12;
		return 3;
	}else if (i==nx && j==0){
		*a_diag = 1./12;
		I[0] = I(i,j+1); a[0] = 1./24;
		I[1] = I(i-1,j); a[1] = 1./24;
		return 2;
	}else if (i==0 && j==ny){
		*a_diag = 1./12;
		I[0] = I(i,j-1);   a[0] = 1./24;
		I[1] = I(i + 1,j); a[1] = 1./24;
		return 2;
	}
		
		
	cout<<"get_offdiag_elem"<<endl;
	abort();
	  
	return -1000;
		
}

int allocate_MSR_matrix(int nx, int ny , double *&p_a , int *&p_I){
	double *a;
	int *I;
	
	
	int N = (nx + 1)*(ny+1);
	int k,l,s = 0;
	int nz = get_non_zeros(nx , ny);
	
	cout<<"nz "<<nz<< " "<<nx <<" "<<ny<<endl;
	
	int len = N + 1 +nz;
	
	a = new double[len];
	if(!a) return -1;
	I = new int[len];
	if(!I) {
		delete []a; return -2;
	}
	I[0] = N + 1;
	
	for(k =1 ; k <=N; k++){
		l = get_num_offdiag(nx , ny , k - 1);
		I[k] = I[k -1] + l;
		s+=l;
		//cout<<"L " << l <<" k"<<k -1<<endl;
	}
	cout<<"COUNT "<<nz <<" "<<s<<" "<<I[N]<<" "<<len;
	assert(nz == s);
	assert(I[N] == len);
	p_a = a;
	p_I = I;
	
	
	return 0;
} 


void build_MSR_matrix(int nx , int ny, double*a, int *I, int p , int k){
	int k1 , k2 , s , sum = 0;
	int N = (nx + 1) * (ny + 1);
	k1 = k*N / p;
	k2 = (k + 1)* N / p ;
	
	for(int  l = k1; l < k2; l++){
		s = get_offdiag_elem(nx , ny , l, a + l , a + I[l] , I + I[l]);
		sum+=s;
		cout<<"L " << s <<" k"<< l<<endl;
	}
	cout<<"N "<<N<<" "<<sum<<" "<<k1<<" "<<k2<<endl;

	//reduce_sum( p , &sum);
	pthread_barrier_wait(&barrier);
	
	cout<<"ASSSERT "<<sum<<" "<<I[N]<<" "<<N + 1 + sum<<endl;
	//assert( N + 1 + sum == I[N]);
	
}

void reduce_sum(int p, int *sum){
    static pthread_mutex_t mu=PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in=PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out=PTHREAD_COND_INITIALIZER;

    static int t_in=0;
    static int t_out=0;

	static int sums = -1;
	/*if(*sum !=0){
		if(sums == -1){
			sums = *sum;
		}else{
			sums+=*sum;
		}
	}*/
    pthread_mutex_lock(&mu);

    t_in++;
    if(t_in>=p){
            t_out=0;
            pthread_cond_broadcast(&c_in);
    }else while(t_in<p) pthread_cond_wait(&c_in,&mu);
	
	/*if(*sum !=0){
		*sum = sums;
	}*/
	
    t_out++;
	
    if(t_out>=p){
			t_in = 0;
            pthread_cond_broadcast(&c_out);
    }else while(t_out<p) pthread_cond_wait(&c_out,&mu);
    pthread_mutex_unlock(&mu);
}










void* msl_approx(void *in_arg) {
    Arg& arg=*((Arg *) in_arg);
  //  static pthread_mutex_t mu=PTHREAD_MUTEX_INITIALIZER;
   // pthread_mutex_lock(&mu);

	/*
    int n_cpus = get_nprocs();

//   cpu_set_t& cpu = *arg.cpu;
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    CPU_SET(n_cpus - arg.ind -1,&cpu);
    pthread_setaffinity_np(pthread_self(),sizeof(cpu),&cpu);

	*/
	
	
  //  pthread_mutex_unlock(&mu);
    const int nx = arg.nx;
    const int ny = arg.ny;
    const int N = (nx + 1)*(ny + 1);
    const int p = arg.p;
    const int thr_ind = arg.thr_ind;

    double *A = arg.A , *b = arg.b;
    int* err = arg.error , *I = arg.I;
	
    cout<<p<<endl;
	int l = thr_ind;
	if(thr_ind == 0) {
		cout<<"Im in "<<endl;
		if(allocate_MSR_matrix(nx, ny, A, I)!= 0){
			cout<<"very bad"<<endl;
			*err = - 1;
		}
		cout<<"Im out"<<endl;
	}
	
	
	//reduce_sum(p);
	
	pthread_barrier_wait(&barrier);
	cout<<"Thread: "<<thr_ind<<endl;
	
	cout<<l<<endl;
	
	/*if(*err == -1){
		cout<<"ERRORRRR"<<endl;
		return 0;
	}*/
	
	cout<<"Ready to built"<<endl;
	
	build_MSR_matrix(nx , ny , A , I, p, thr_ind);
	
	reduce_sum(p);
	if(thr_ind == 0){
		//print_MSR_matrix(A,I,N);
		print_matrix(N , A, I);
		cout<< l<<endl;
	}
	
	
	return 0 ;
}

/*
double get_time(){
        struct rusage buf;
        getrusage(RUSAGE_THREAD,&buf);
        return (double)buf.ru_utime.tv_sec+(double)buf.ru_utime.tv_usec/1000000.;
}

double get_full_time(){
        struct timeval buf;
        gettimeofday(&buf,0);
        return (double)buf.tv_sec+(double)buf.tv_usec/1000000.;
}

*/