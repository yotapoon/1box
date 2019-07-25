#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <windows.h>
#define N 100
double a = 0.01;//radius of disk
double Xmin = -1.6,Xmax = 1.6,X0 = 0.0;
double Ymin = 0.0,Ymax = 1.0;
double V0 = 1.0;
int N_cell_x = 16,N_cell_y = 12;


struct PARTICLE{
	double x,y,u,v;
	double tau;//粒子の固有時間を記録，Delayed State Algorithm(DSA)による高速化のために必要
	double next;
};

struct CELL{
	int first[N_cell_x][N_cell_y];
	int last[N_cell_x][N_cell_y];
	double xlength = (Xmax-Xmin)/(double)N_cell_x,ylength = (Ymax-Ymin)/(double)N_cell_y;
};

void status_initialize(struct PARTICLE particle[N]);
void cell_register(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y]);
int set(struct PARTICLE partile[N],int i);
double distance(struct PARTICLE particle[N],int i,int j);
double Uniform(void);
double rand_normal( double mu, double sigma );
int getcell_x(double x,double cell_length_x);//ellはcell_lengthのこと
int getcell_y(double y,double cell_length_y);


int main(void){
	
	struct PARTICLE particle[N];
	struct CELL cell[N_cell_x+1][N_cell_y+1];
	
	status_initialize(particle);//位置や速度の初期化
	cell_register(particle,cell);//粒子をセルに登録する
	
	for(int i=0;i<N;i++){
		printf("%lf %lf\n",particle[i].x,particle[i].u);
	}
	return 0;
}

void status_initialize(struct PARTICLE particle[N]){
	double prob;
	int i;
	for(i=0;i<N;i++){
		prob = Uniform();
		particle[i].x = (Xmin+a)*prob+(Xmax-a)*(1-prob);
		particle[i].y = (Ymin+a)*prob+(0.5*Ymax-a)*(1-prob);
		while(set(particle,i) == 0){
			prob = Uniform();
			particle[i].x = (Xmin+a)*prob+(Xmax-a)*(1-prob);
			prob = Uniform();
			particle[i].y = (Ymin+a)*prob+(Ymax-a)*(1-prob);
		}
		particle[i].u = rand_normal(0.0,V0);
		particle[i].v = rand_normal(0.0,V0);
		particle[i].next = -1;
	}
}

void cell_register(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y]){
	int i,cell_x,cell_y,lastPrev;
	//initialize particle.next and cell
	for(i=0;i<N;i++){
		particle[i].next = -1;
	}
	for(cell_x = 1;cell_x <= N_cell_x;cell_x++){
		for(cell_y = 1;cell_y <= N_cell_y;cell_y++){
			cell[cell_x][cell_y].first = -1;
			cell[cell_x][cell_y].last = -1;
		}
	}
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	for(i=0;i<N;i++){
		cell_x = getcell_x(particle[i].x,cell_length_x);
		cell_y = getcell_y(particle[i].y,cell_length_y);
		lastPrev = cell[cell_x][cell_y].last;
		cell[cell_x][cell_y].last = i;
		
		if(lastPrev == -1){
			cell[cell_x][cell_y].first = i;
		}else{
			particle[lastPrev].next = i;
		}
	}
}


int set(struct PARTICLE particle[N],int i){//setに成功していれば1,失敗していれば0を返す
	int j,r=1;
	double d;
	
	if(fabs(particle[i].x) < a){
		r = 0;
	}
	for(j=1;j<=i-1;j++){
		d = distance(particle,i,j);
		if(d <= 2.0*a){
			r = 0;
			break;
		}
	}
	return r;
}

double distance(struct PARTICLE particle[N],int i,int j){
	double d;
	d = sqrt(pow(particle[i].x-particle[j].x,2.0)+pow(particle[i].y-particle[j].y,2.0));
	return d;
}

double Uniform(void){
	return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

double rand_normal( double mu, double sigma ){
	double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
	return mu + sigma*z;
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