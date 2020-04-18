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



int main( ){
	parral par(2 , 3 , 30 , 0.2, 3, 3);

	double *a, *b, *x;
	int * I;
    pthread_t *tids;
    Arg *args;

    int nx = 3 , ny = 3 , p = 16, *error;


	
    tids= new pthread_t[p];
    args= new Arg[p];


	pthread_barrier_init (&barrier, nullptr, p);

    for(int i=0;i<p;i++){
        args[i].p=p;
        args[i].thr_ind=i;
        args[i].A=&a;
        args[i].b=&b;
		args[i].I=&I;
        args[i].x=x;
        args[i].error=error;
		args[i].nx = par.nx;
		args[i].ny = par.ny;
		args[i].f = f_3;
		args[i].par = &par;
	}

        //double TIME = get_full_time();

	cout<<"Lets GO"<<endl;
    for(int i=0;i<p;i++){
           if(pthread_create(tids+i,0,&msl_approx,(void *) (args+i))){
                    printf("Cannot create thread %d\n",i);
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
	cout<<"bue"<<endl;
	
	

    //delete[] x;
	delete[] tids;
	delete[] args;
		
		
	
	return 0;
}