#include "mls.h"



double f_1(double , double){
	return 1;
}
double f_2(double x, double){
	return x;
}
double f_3(double , double y){
	return y;
}
double f_4(double x, double y){
	return x + y;
}

double f_5(double x, double y){
	return x*y + x*x ;
}

double f_6(double x, double y){
	return exp(x)*exp(y);
}




int main( int argc , char** argv){
	int nx = 1024 , ny = 1024, p = 1, *error = new int(0);
	parral par(3 , 3 , 90, 0.9 , nx, ny);

	double *a, *b, *x , *v , *u , *r , *buf ;
	int * I;
    pthread_t *tids;
    Arg *args;

    

	if (!(argc == 2) || ((p = stoi(argv[1])) <= 0)) {
            printf("usage: %s p \n ", argv[0]);
	    return -1;
	}
	
	
	int N = (nx + 1)*(ny + 1) - par.nx_rect*par.ny_rect;

	cout<<par.nx<<" "<<par.ny<<" "<<par.nx_rect<<" "<<par.ny_rect<<endl;
	
    tids= new pthread_t[p];
    args= new Arg[p];
	buf = new double[p];
	x = new double[N];
	u = new double[N];
	v = new double[N];
	r = new double[N];
	
	memset(x , 0 , N*sizeof(double));
	memset(u , 0 , N*sizeof(double));
	memset(v , 0 , N*sizeof(double));
	memset(r , 0 , N*sizeof(double));
	memset(buf , 0 , p*sizeof(double));

	// pthread_barrier_init (&barrier, nullptr, p);

    for(int i=0;i<p;i++){
        args[i].p=p;
        args[i].thr_ind=i;
        args[i].A=&a;
        args[i].b=&b;
		args[i].I=&I;
        args[i].x=x;
        args[i].u=u;
        args[i].v=v;
        args[i].r=r;
		args[i].buf = buf;
        args[i].error=error;
		args[i].nx = par.nx;
		args[i].ny = par.ny;
		args[i].f = f_4;
		args[i].par = &par;
	}

        //double TIME = get_full_time();

	cout<<"Lets GO"<<endl;
    for(int i=0;i<p;i++){
           if(pthread_create(tids+i,0,&msl_approx,(void *) (args+i))){
                    printf("Cannot create thread %d\n",i);
					delete[] tids;
					delete[] args;
					delete[] x;
					delete[] u;
					delete[] v;
					delete[] r;
					delete error;
                    abort();
            }
    }

	for(int i=0;i<p;i++)pthread_join(tids[i],0);
	/*
    
    TIME= get_full_time() - args[0].fulltime ;
	
	
	
    for(int i = 0; i<p;i++){
        cout<<"TIME of "<<args[i].ind<<" potok "<<args[i].time_thr<<endl;
    }
    LN;

    printf("FULL TIME IS %f\n",TIME);

	*/
	//cout<<"bue"<<endl;
	
	

    //delete[] x;
	delete[] tids;
	delete[] args;
	delete[] x;
	delete[] u;
	delete[] v;
	delete[] r;
	delete error;
	
		
	
	return 0;
}