#include "mls.h"





int main( ){
	double xs[] =  {0 , 1 , 2 , 3};
	double ys[] =  {0 , 1 , 2 , 3};
	
	double x_new[4];
	double y_new[4];
	
	for( int i = 0 ; i < 4; i ++){
		coord_trans(xs[i] , ys[i] , x_new[i] , y_new[i]);
		cout<<x_new[i]<<" "<<y_new[i]<<endl;
	}
	

	double *a, *b, *x;
	int * I;
    pthread_t *tids;
    Arg *args;

    int nx = 3 , ny = 3 , p =3, *error;


	
    tids= new pthread_t[p];
    args= new Arg[p];


	pthread_barrier_init (&barrier, nullptr, p);

    for(int i=0;i<p;i++){
        args[i].p=p;
        args[i].thr_ind=i;
        args[i].A=a;
        args[i].b=b;
		args[i].I=I;
        args[i].x=x;
        args[i].error=error;
		args[i].nx = nx;
		args[i].ny = ny;
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
	
	
    delete[] a;
    delete[] b;
    //delete[] x;
	delete[] tids;
	delete[] args;
		
		
	
	return 0;
}