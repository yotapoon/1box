#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define N 720
#define n ((int)ceil(log2(N)))
#define q ((int)pow(2,n)-N)
#define p ((N-q)/2)
#define DEFAULT_INPUT 0.0

//scaleはすべてcm
double a = 0.01;//radius of disk
double V0 = 0.5;
double e = 0.95,e_wall = 1.0;
double g = 1.0;
int N_cell_x = 64,N_cell_y = 12;
double Xmin = -3.2,Xmax = 3.2,X0 = 0.0;
double Ymin = 0.0,Ymax = 1.0;
double ell = 0.05;
double U = 0.149;//箱の運動を決めるパラメータ
double T = 3000.0;
double T_demon = 0.0;
double epsilon = 0.000001;
double p_demon = 0.00;



typedef struct node{
	int number;
	double time;
	struct node *left;
	struct node *right;
	struct node *parent;
} Node;

struct event{
	double time;
	int number_particle;
	int number_col;
};

struct PARTICLE_CELL{
	int first;
	int last;
};

int getcell_x(double x,double cell_length_x);//ellはcell_lengthのこと
int getcell_y(double y,double cell_length_y);
void initialize(struct PARTICLE_CELL particle_cells[N_cell_x+1][N_cell_y+1],int nextof[N]);

void status_initialize(double x[N],double y[N],double u[N],double v[N],double tau[N],int C[N]);
double distance(double x1,double y1,double x2,double y2);
double Uniform(void);
int set(double x[N],double y[N],int i);//setに成功していれば1,失敗していれば0を返す
double rand_normal( double mu, double sigma );
int intpow(int a,int b);
void CBT_initialize(Node *entry[n+1][2*p+2*q],struct event L[N]);
void CBT_update(Node *entry[n+1][2*p+2*q],double time_new,int i_new);
void initialize_LM(double x[N],double y[N],double u[N],double v[N],struct event *L,int i,int C[N]);
void Free_evolution(double *x,double *y,double *u,double *v,double t);
void G1(double *x,double *y,double *u,double *v,int j,double h);
void G2(double *x1,double *y1,double *u1,double *v1,double *x2,double *y2,double *u2,double *v2);
void Predictions(double x[N],double y[N],double u[N],double v[N],struct event *L,int i,double t,struct PARTICLE_CELL particle_cells[N_cell_x+1][N_cell_y+1],int nextof[N+1],double cell_length_x,double cell_length_y,double tau[N]);
double T_DDC(double x1,double y1,double u1,double v1,double x2,double y2,double u2,double v2,double t);
double T_DWC(double x,double y,double u,double v,double t,int j);
double m_cal(double x[N]);
double Vmax(double u[N],double v[N]);
double nonamin(int x,int y);
double nonamin(int x,int y);

int main(int argc,char* argv[]){
	double h = 2 <= argc ? atof(argv[1]) : DEFAULT_INPUT;
	
	double x[N],y[N];
	double u[N],v[N];
	double tau[N];//粒子の固有時間を記録、,DelayedStateAlgorithm(DSA)による高速化のために必要
	int C[N];//粒子のこれまでの衝突回数を記録、eventがinvalidか判定するために必要
	double t=0.0,dt=0.01,trec_start = T*0.5,trec=trec_start,dtrec = (double)T/10000.0,t0,t_old;
	double t_print=0.0,dt_print = T/100.0;
	double t_cell=0.0,dt_cell,t_cell_old;
	double pi = atan(1.0)*4.0,theta;
	int i,j,k=0,loop_judge,i_current,j_current;
	double unit;
	int nextof[N];
	struct PARTICLE_CELL particle_cells[N_cell_x+1][N_cell_y+1];
	int cell_x,cell_y,lastPrev,c1,c2;
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x;
	double cell_length_y =(Ymax-Ymin)/(double)N_cell_y;
	double m,m_max,m_min,m_average,v_max;
	//double h_start=0.8,h_end=1.0,dh=0.02,h;
	int counter=1,counter_end = 10;
	//printf("cell_length = %lf %lf\n",cell_length_x,cell_length_y);
	initialize(particle_cells,nextof);
	
	//realization
	Node *entry[n+1][2*p+2*q],hogee[n+1][2*p+2*q];
	for(i=0;i<=n;i++){
		for(j=0;j<2*p+2*q;j++){
			entry[i][j] = &hogee[i][j];
		}
	}
	
	srand((unsigned) time(NULL));
	struct event L[N];
	
	
	FILE *fp_m;
	char name_m[256],name_mplot[256];
	sprintf(name_m,"m(p=%lf,h=%lf).txt",p_demon,h);
	if((fp_m = fopen(name_m,"w"))==NULL){
		printf("file_m open error\n");
	}
	m_average = 0.0;
	m_max = 0.0;
	m_min = 1.0;
	counter = 1;
	while(counter <= counter_end){
		printf("h=%lf (%d times)\n",h,counter);
		initialize(particle_cells,nextof);
		//初期条件の設定
		t_old = 0.0;
		t_cell = 0.0;
		t = 0.0;
		//初期条件の設定
		status_initialize(x,y,u,v,tau,C);
		for(i=0;i<N;i++){
			//セルへの登録
			cell_x = getcell_x(x[i],cell_length_x);
			cell_y = getcell_y(y[i],cell_length_y);
			
			lastPrev = particle_cells[cell_x][cell_y].last;
			particle_cells[cell_x][cell_y].last = i;
			
			if(lastPrev == -1){
				particle_cells[cell_x][cell_y].first = i;
			}else{
				nextof[lastPrev] = i;
			}
		}
		v_max = Vmax(u,v);

		dt_cell = (cell_length_y-2.0*a)/(2.0*v_max);
		t_cell = (cell_length_y-2.0*a)/(2.0*v_max);
		t_cell_old = 0.0;
		
		for(i=0;i<N;i++){
			Predictions(x,y,u,v,&L[i],i,t,particle_cells,nextof,cell_length_x,cell_length_y,tau);
		}
		
		CBT_initialize(entry,L);
		//反射とinitialize_LMは必要,update_scheduleはwhileの中で
		printf("set up ok\n");
		trec = trec_start;
		m = 0.0;
		t_old = t;
		t_print = dt_print;
		
		//時間発展
		while(t<=T){
			if(t > t_print){
				printf("t = %lf:m = %lf\n",t,m_cal(x));
				t_print += dt_print;
			}
			//NEXT EVENTの検索
			i_current = entry[0][0]->number;
			j_current = L[i_current].number_particle;
			
			C[i_current] += 1;
			t = L[i_current].time;
			if(j_current >= 0){
				C[j_current] += 1;
			}
			Free_evolution(&x[i_current],&y[i_current],&u[i_current],&v[i_current],t-tau[i_current]);
			tau[i_current] = t;
			if(j_current >= 0){//DDC
				Free_evolution(&x[j_current],&y[j_current],&u[j_current],&v[j_current],t-tau[j_current]);
				tau[j_current] = t;
				G2(&x[i_current],&y[i_current],&u[i_current],&v[i_current],&x[j_current],&y[j_current],&u[j_current],&v[j_current]);
			}
			if(j_current < 0){//DWC
				G1(&x[i_current],&y[i_current],&u[i_current],&v[i_current],j_current,h);
				if(j_current == -3){
					if(v_max*v_max < u[i_current]*u[i_current]+v[i_current]*v[i_current]){
						v_max = sqrt(u[i_current]*u[i_current]+v[i_current]*v[i_current]);
					}
				}
				dt_cell = (cell_length_y-2.0*a)/(2.0*v_max);
				t_cell = t_cell_old+dt_cell;
			}
			Predictions(x,y,u,v,&L[i_current],i_current,t,particle_cells,nextof,cell_length_x,cell_length_y,tau);
			CBT_update(entry,L[i_current].time,i_current);
			if(j_current >= 0){//DDC
				Predictions(x,y,u,v,&L[j_current],j_current,t,particle_cells,nextof,cell_length_x,cell_length_y,tau);
				CBT_update(entry,L[j_current].time,j_current);
			}
			
			//fprintf(fp_check,"time:%lf\n",t);
			//fprintf(fp_check,"i=%d:%lf %lf %lf %lf,%lf %d\n",i_current,x[i_current],y[i_current],u[i_current],v[i_current],L[i_current].time,L[i_current].number_particle);
			if(j_current >= 0){
				//fprintf(fp_check,"i=%d:%lf %lf %lf %lf,%lf %d\n",j_current,x[j_current],y[j_current],u[j_current],v[j_current],L[j_current].time,L[j_current].number_particle);
			}
			
			//i_current,j_currentと衝突する粒子がいた場合はその粒子のeventはinvalidになってしまう
			//そのような粒子は同じマスク内にしか存在しないはずなのでその中で探索
			cell_x = getcell_x(x[i_current],cell_length_x);
			cell_y = getcell_y(y[i_current],cell_length_y);
			
			for(c1=-1;c1<=1;c1++){
				for(c2=-1;c2<=1;c2++){
					if((((cell_x+c1 >= 1) && (cell_x+c1 <= N_cell_x)) && (cell_y+c2 >= 1)) && (cell_y+c2 <= N_cell_y)){
						j = particle_cells[cell_x+c1][cell_y+c2].first;
						while(j >= 0){
							if(L[j].number_particle == i_current){
								Predictions(x,y,u,v,&L[j],j,t,particle_cells,nextof,cell_length_x,cell_length_y,tau);
								CBT_update(entry,L[j].time,j);
								//fprintf(fp_check,"(invalid)j=%d:%lf %lf %lf %lf,%lf %d\n",j,x[j],y[j],u[j],v[j],L[j].time,L[j].number_particle);
							}
							j = nextof[j];
						}
					}
				}
			}
			if(j_current >= 0){
				cell_x = getcell_x(x[j_current],cell_length_x);
				cell_y = getcell_y(y[j_current],cell_length_y);
				
				for(c1=-1;c1<=1;c1++){
					for(c2=-1;c2<=1;c2++){
						if((((cell_x+c1 >= 1) && (cell_x+c1 <= N_cell_x)) && (cell_y+c2 >= 1)) && (cell_y+c2 <= N_cell_y)){
							j = particle_cells[cell_x+c1][cell_y+c2].first;
							while(j >= 0){
								if(L[j].number_particle == j_current){
									Predictions(x,y,u,v,&L[j],j,t,particle_cells,nextof,cell_length_x,cell_length_y,tau);
									CBT_update(entry,L[j].time,j);
									//fprintf(fp_check,"(invalid)j=%d:%lf %lf %lf %lf,%lf %d\n",j,x[j],y[j],u[j],v[j],L[j].time,L[j].number_particle);
								}
								j = nextof[j];
							}
						}
					}
				}
			}
			//fprintf(fp_check,"\n\n");
			//EEPGM
			if(t >= t_cell){
				//printf("t_cell = %lf\n",t_cell);
				for(i=0;i<N;i++){
					Free_evolution(&x[i],&y[i],&u[i],&v[i],t-tau[i]);
					tau[i] = t;
				}
				initialize(particle_cells,nextof);
				for(i=0;i<N;i++){
					//セルへの登録
					cell_x = getcell_x(x[i],cell_length_x);
					cell_y = getcell_y(y[i],cell_length_y);
					
					lastPrev = particle_cells[cell_x][cell_y].last;
					particle_cells[cell_x][cell_y].last = i;
					
					if(lastPrev == -1){
						particle_cells[cell_x][cell_y].first = i;
					}else{
						nextof[lastPrev] = i;
					}
				}
				for(i=0;i<N;i++){
					Predictions(x,y,u,v,&L[i],i,t,particle_cells,nextof,cell_length_x,cell_length_y,tau);
				}
				CBT_initialize(entry,L);
				
				v_max = Vmax(u,v);
				t_cell_old = t_cell;
				dt_cell = (cell_length_y-2.0*a)/(2.0*v_max);
				t_cell = t_cell_old+dt_cell;
				
				for(i=0;i<N;i++){
					if(y[i] < Ymin+a-epsilon){
						//fprintf(fp_check,"(error)i=%d:%lf %lf %lf %lf,%lf %d\n",i,x[i],y[i],u[i],v[i],L[i].time,L[i].number_particle);
						printf("i=%d:error\n",i);
						printf("%lf %lf %lf %lf\n",x[i],y[i],u[i],v[i]);
						printf("%lf %d %d\n",L[i].time,L[i].number_particle,L[i].number_col);
						G1(&x[i],&y[i],&u[i],&v[i],-3,h);
						Predictions(x,y,u,v,&L[i],i,t,particle_cells,nextof,cell_length_x,cell_length_y,tau);
						CBT_update(entry,L[i].time,i);
						
						cell_x = getcell_x(x[i],cell_length_x);
						cell_y = getcell_y(y[i],cell_length_y);
						for(c1=-1;c1<=1;c1++){
							for(c2=-1;c2<=1;c2++){
								if((((cell_x+c1 >= 1) && (cell_x+c1 <= N_cell_x)) && (cell_y+c2 >= 1)) && (cell_y+c2 <= N_cell_y)){
									j = particle_cells[cell_x+c1][cell_y+c2].first;
									while(j >= 0){
										Predictions(x,y,u,v,&L[j],j,t,particle_cells,nextof,cell_length_x,cell_length_y,tau);
										CBT_update(entry,L[j].time,j);
										j = nextof[j];
									}
								}
							}
						}
					}
				}
				//fprintf(fp_check,"EEPGM\n");
				for(j=0;j<N;j++){
					//fprintf(fp_check,"(invalid)j=%d:%lf %lf %lf %lf,%lf %d\n",j,x[j],y[j],u[j],v[j],L[j].time,L[j].number_particle);
				}
				//fprintf(fp_check,"\n\n");
			}
			
			
			if(t > trec){
				
				for(i=0;i<N;i++){
					Free_evolution(&x[i],&y[i],&u[i],&v[i],t-tau[i]);
					tau[i] = t;
				}
				initialize(particle_cells,nextof);
				for(i=0;i<N;i++){
					//セルへの登録
					cell_x = getcell_x(x[i],cell_length_x);
					cell_y = getcell_y(y[i],cell_length_y);
					
					lastPrev = particle_cells[cell_x][cell_y].last;
					particle_cells[cell_x][cell_y].last = i;
					
					if(lastPrev == -1){
						particle_cells[cell_x][cell_y].first = i;
					}else{
						nextof[lastPrev] = i;
					}
				}
				for(i=0;i<N;i++){
					Predictions(x,y,u,v,&L[i],i,t,particle_cells,nextof,cell_length_x,cell_length_y,tau);
				}
				CBT_initialize(entry,L);
				
				v_max = Vmax(u,v);
				t_cell_old = t_cell;
				dt_cell = (cell_length_y-2.0*a)/(2.0*v_max);
				t_cell = t_cell_old+dt_cell;
				
				for(i=0;i<N;i++){
					if((y[i] < Ymin+a-epsilon)||(fabs(x[i]) > Xmax-a+epsilon)){
						//fprintf(fp_check,"(error)i=%d:%lf %lf %lf %lf,%lf %d\n",i,x[i],y[i],u[i],v[i],L[i].time,L[i].number_particle);
						printf("i=%d:error\n",i);
						printf("%lf %d %d\n",L[i].time,L[i].number_particle,L[i].number_col);
						printf("%lf %lf %lf %lf\n",x[i],y[i],u[i],v[i]);
						G1(&x[i],&y[i],&u[i],&v[i],-3,h);
						Predictions(x,y,u,v,&L[i],i,t,particle_cells,nextof,cell_length_x,cell_length_y,tau);
						CBT_update(entry,L[i].time,i);
						
						cell_x = getcell_x(x[i],cell_length_x);
						cell_y = getcell_y(y[i],cell_length_y);
						for(c1=-1;c1<=1;c1++){
							for(c2=-1;c2<=1;c2++){
								if((((cell_x+c1 >= 1) && (cell_x+c1 <= N_cell_x)) && (cell_y+c2 >= 1)) && (cell_y+c2 <= N_cell_y)){
									j = particle_cells[cell_x+c1][cell_y+c2].first;
									while(j >= 0){
										Predictions(x,y,u,v,&L[j],j,t,particle_cells,nextof,cell_length_x,cell_length_y,tau);
										CBT_update(entry,L[j].time,j);
										j = nextof[j];
									}
								}
							}
						}
					}
				}
				v_max = Vmax(u,v);
				m += m_cal(x)*(dtrec/(double)(T-trec_start));
				trec += dtrec;
			}
		}
		m = fabs(m);
		printf("m = %lf\n",m);
		m_average += m/(double)counter_end;
		m_max = ((m+m_max)+fabs(m-m_max))/2.0;
		m_min = ((m+m_min)-fabs(m-m_min))/2.0;
		counter += 1;
	}
	fprintf(fp_m,"%lf %lf %lf %lf\n",h,m_average,m_min,m_max);
	fclose(fp_m);
	return 0;
}

void status_initialize(double x[N],double y[N],double u[N],double v[N],double tau[N],int C[N]){
	double prob;
	int i;
	for(i=0;i<N;i++){
		prob = Uniform();
		x[i] = (Xmin+a)*prob+(Xmax-a)*(1-prob);
		prob = Uniform();
		y[i] = (Ymin+a)*prob+(0.5*Ymax-a)*(1-prob);
		while(set(x,y,i) == 0){
			prob = Uniform();
			x[i] = (Xmin+a)*prob+(Xmax-a)*(1-prob);
			prob = Uniform();
			y[i] = (Ymin+a)*prob+(Ymax-a)*(1-prob);
		}
		C[i] = 0;
		tau[i] = 0.0;
		u[i] = rand_normal(0.0,V0);
		v[i] = rand_normal(0.0,V0);
	}
}


double distance(double x1,double y1,double x2,double y2){
	double d;
	d = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
	return d;
}


double Uniform(void){
	return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}


int set(double x[N],double y[N],int i){//setに成功していれば1,失敗していれば0を返す
	int j,r=1;
	double d;
	
	if(fabs(x[i]) < a){
		r = 0;
	}
	for(j=1;j<=i-1;j++){
		d = sqrt(pow(x[j]-x[i],2.0)+pow(y[j]-y[i],2.0));
		if(d <= 2.0*a){
			r = 0;
			break;
		}
	}
	return r;
}

double rand_normal( double mu, double sigma ){
	double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
	return mu + sigma*z;
}



int intpow(int a,int b){
	return (int)pow(a,b);
}
void CBT_initialize(Node *entry[n+1][2*p+2*q],struct event L[N]){
	int i,n_index;
	//initialization for bottom nodes
	for(i=0;i<2*p+2*q;i++){
		if(i < 2*p){
			entry[n][i]->time = L[i].time;
			entry[n][i]->number = i;
		}else{
			entry[n][i]->time = L[2*p+(i-2*p)/2].time;
			entry[n][i]->number = 2*p+(i-2*p)/2;
		}
	}
	//tournament
	for(n_index=n-1;n_index>=0;n_index--){
		for(i=0;i<=intpow(2,n_index)-1;i++){
			entry[n_index][i]->left = entry[n_index+1][2*i];
			entry[n_index+1][2*i]->parent = entry[n_index][i];
			entry[n_index][i]->right = entry[n_index+1][2*i+1];
			entry[n_index+1][2*i+1]->parent = entry[n_index][i];
			if(entry[n_index+1][2*i]->time <= entry[n_index+1][2*i+1]->time){
				entry[n_index][i]->time = entry[n_index+1][2*i]->time;
				entry[n_index][i]->number = entry[n_index+1][2*i]->number;
			}else{
				entry[n_index][i]->time = entry[n_index+1][2*i+1]->time;
				entry[n_index][i]->number = entry[n_index+1][2*i+1]->number;
			}
		}
	}
}

void CBT_update(Node *entry[n+1][2*p+2*q],double time_new,int i_new){
	Node *entry_now,hoge_now;
	if(i_new < 2*p){
		entry[n][i_new]->time = time_new;
		entry_now = entry[n][i_new];
	}else{//practiceはここでミスっていた、修正すること
		entry[n][2*i_new-2*p]->time = time_new;
		entry[n][2*i_new-2*p+1]->time = time_new;
		entry_now = entry[n][2*i_new-2*p];
	}
	while(entry_now->parent != NULL){
		entry_now = entry_now->parent;
		if(entry_now->left->time < entry_now->right->time){
			entry_now->time = entry_now->left->time;
			entry_now->number = entry_now->left->number;
		}else{
			entry_now->time = entry_now->right->time;
			entry_now->number = entry_now->right->number;
		}
	}
}
void Predictions(double x[N],double y[N],double u[N],double v[N],struct event *L,int i,double t,struct PARTICLE_CELL particle_cells[N_cell_x+1][N_cell_y+1],int nextof[N+1],double cell_length_x,double cell_length_y,double tau[N]){
	int j,c1,c2;
	double d;
	double Cx,Cy;
	int i_col,j_col;
	int cell_x,cell_y;
	double r_relative,v_relative,b,hoge;
	double t_min = 2.0*T,t_temp;
	double xj,yj,uj,vj;
	
	for(j=-4;j<0;j++){
		t_temp = T_DWC(x[i],y[i],u[i],v[i],tau[i],j);//t = tau[]よりok
		if((t_temp > t) && (t_temp < t_min)){
			t_min = t_temp;
			j_col = j;
		}
	}
	
	cell_x = getcell_x(x[i],cell_length_x);
	cell_y = getcell_y(y[i],cell_length_y);
	
	for(c1=-1;c1<=1;c1++){
		for(c2=-1;c2<=1;c2++){
			if((((cell_x+c1 >= 1) && (cell_x+c1 <= N_cell_x)) && (cell_y+c2 >= 1)) && (cell_y+c2 <= N_cell_y)){
				j = particle_cells[cell_x+c1][cell_y+c2].first;
				while(j >= 0){
					xj = x[j];
					yj = y[j];
					uj = u[j];
					vj = v[j];
					Free_evolution(&xj,&yj,&uj,&vj,tau[i]-tau[j]);
					t_temp = T_DDC(x[i],y[i],u[i],v[i],xj,yj,uj,vj,tau[i]);
					if((t_temp > t) && (t_temp < t_min)){
						t_min = t_temp;
						j_col = j;
					}
					j = nextof[j];
				}
			}
		}
	}
	L->time = t_min;
	L->number_particle = j_col;
}



void Free_evolution(double *x,double *y,double *u,double *v,double t){
	*x += (*u)*t;
	*y += (*v)*t-0.5*g*t*t;
	*v += -g*t;
}

void G1(double *x,double *y,double *u,double *v,int j,double h){
	double temp;
	if((j == -1) || (j == -2)){//collision with R or L wall
		*u = -e_wall**u;
		if(j == -1){
			*x = Xmax-a-epsilon;
		}else{
			*x = Xmin+a+epsilon;
		}
	}else if(j == -3){//collision with Bottom wall
		*v = (1+e_wall)*U-e_wall**v;
		*y = Ymin+a+epsilon;
	}else if(j == -4){//collision with Cener wall
		if(*x > 0){
			*x = a;
		}else if(*x < -0){
			*x = -a;
		}
		if(fabs(*y-h) > ell){//really collision
			*u = -e_wall**u;
		}else{
			temp = Uniform();
			if((temp < p_demon)&&(*x > 0.0)){//demon prevent the particle from passing through right-to-left
				*u = -e_wall**u;
			}
		}
	}
}


	
	
void G2(double *x1,double *y1,double *u1,double *v1,double *x2,double *y2,double *u2,double *v2){
	double d,Xtemp,Ytemp,Utemp1,Utemp2,Vtemp1,Vtemp2,Cx,Cy;
	d = distance(*x1,*y1,*x2,*y2);
	Utemp1 = *u1;
	Vtemp1 = *v1;
	Utemp2 = *u2;
	Vtemp2 = *v2;
	Cx = (*x1-*x2)/d;
	Cy = (*y1-*y2)/d;
	*u1 = 0.5*(1+e)*((Utemp2-Utemp1)*Cx+(Vtemp2-Vtemp1)*Cy)*Cx+Utemp1;
	*v1 = 0.5*(1+e)*((Utemp2-Utemp1)*Cx+(Vtemp2-Vtemp1)*Cy)*Cy+Vtemp1;
	*u2 = 0.5*(1+e)*((Utemp1-Utemp2)*Cx+(Vtemp1-Vtemp2)*Cy)*Cx+Utemp2;
	*v2 = 0.5*(1+e)*((Utemp1-Utemp2)*Cx+(Vtemp1-Vtemp2)*Cy)*Cy+Vtemp2;
	
//	d = distance(*x1,*y1,*x2,*y2);
//	Utemp = *u1;
//	Vtemp = *v1;
//	Cx = (*x1-*x2)/d;
//	Cy = (*y1-*y2)/d;
//	*u1 = 0.5*(1+e)*((*u2-*u1)*Cx+(*v2-*v1)*Cy)*Cx+*u1;
//	*v1 = 0.5*(1+e)*((*u2-*u1)*Cx+(*v2-*v1)*Cy)*Cy+*v1;
//	*u2 = 0.5*(1+e)*((Utemp-*u2)*Cx+(Vtemp-*v2)*Cy)*Cx+*u2;
//	*v2 = 0.5*(1+e)*((Utemp-*u2)*Cx+(Vtemp-*v2)*Cy)*Cy+*v2;
//	Xtemp = *x1;
//	Ytemp = *y1;
//	*x1 = *x1*2.0*(a/d)+*x2*(1-2.0*a/d);
//	*y1 = *y1*2.0*(a/d)+*y2*(1-2.0*a/d);
//	*x2 = *x2*2.0*(a/d)+Xtemp*(1-2.0*a/d);
//	*y2 = *y2*2.0*(a/d)+Ytemp*(1-2.0*a/d);
}

double T_DDC(double x1,double y1,double u1,double v1,double x2,double y2,double u2,double v2,double t){
	double r_relative,v_relative,b,hoge;
	double tau = t;
	r_relative = distance(x1,y1,x2,y2);
	v_relative = distance(u1,v1,u2,v2);
	b = (x1-x2)*(u1-u2)+(y1-y2)*(v1-v2);
	hoge = b*b-v_relative*v_relative*(r_relative*r_relative-4.0*a*a);
	if(hoge > 0){
		tau += -(b+sqrt(hoge))/(v_relative*v_relative);
	}else{
		tau += T;
	}
	return tau;
}

double T_DWC(double x,double y,double u,double v,double t,int j){
	double tau = t;
	//collision with RIGHT wall(-1)
	if(j==-1){
		if(u>0.0){
			tau += (Xmax-a-x)/u;
		}else{
			tau += 2.0*T;
		}
	}else if(j==-2){
		if(u<0.0){
			tau += (Xmin+a-x)/u;
		}else{
			tau += 2.0*T;
		}
	}else if(j==-3){
		tau += (v+sqrt(v*v+2*g*(y-a)))/g;
	//	if(v<0.0){
	//		tau += (v+sqrt(v*v+2*g*(y-a)))/g;
	//	}else{
	//		tau += 2.0*T;
	//	}
	}else if(j==-4){
		if(x > a){
			tau += -(x-a)/u;
		}else if(x < -a){
			tau += (-x-a)/u;
		}
	}
	if(tau < t){
		return 2.0*T;
	}else{
		return tau;
	}
}


void initialize(struct PARTICLE_CELL particle_cells[N_cell_x+1][N_cell_y+1],int nextof[N]){
	int i,j;
	for(i=1;i<=N_cell_x;i++){
		for(j=1;j<=N_cell_y;j++){
			particle_cells[i][j].first = -1;
			particle_cells[i][j].last = -1;
		}
	}
	for(i=0;i<N;i++){
		nextof[i] = -1;
	}
}

int getcell_x(double x,double cell_length_x){
	
	if((x < Xmin+a)||(Xmax-a < x)){
//		printf("x is out of range\n");
	}
	if(x > 0){
		return N_cell_x/2+1+(int)(x/cell_length_x);
	}else{
		return N_cell_x/2+(int)(x/cell_length_x);
	}
	
}

int getcell_y(double y,double cell_length_y){
	if(y < Ymin){
		//printf("error:y<0(%lf)\n",y);
		//return -1;
		return 1;
	}else if(y>Ymax){
		return N_cell_y;
	}else{
		return (int)(y/cell_length_y)+1;
	}
}

double nonamax(int x,int y){
	if(x > y){
		return (double)x;
	}else{
		return (double)y;
	}
}

double nonamin(int x,int y){
	if(x < y){
		return (double)x;
	}else{
		return (double)y;
	}
}

double m_cal(double x[N]){
	int i,N_left=0,N_right;
	for(i=0;i<N;i++){
		if(x[i] < X0){
			N_left += 1;
		}
	}
	double m;
	N_right = N-N_left;
	m = (nonamax(N_left,N_right)-0.5*N)/(double)N;
	return m;
}

double Vmax(double u[N],double v[N]){
	int i;
	double v_max = 0.0;
	for(i=0;i<N;i++){
		if(v_max*v_max < u[i]*u[i]+v[i]*v[i]){
			v_max = sqrt(u[i]*u[i]+v[i]*v[i]);
		}
	}
	return v_max;
}