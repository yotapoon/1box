#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <windows.h>
#define N 300
#define n ((int)ceil(log2(N)))
#define q ((int)pow(2,n)-N)
#define p ((N-q)/2)
double a = 0.02;//radius of disk
double e = 0.95,e_wall = 1.0;
double g = 1.0;
double Xmin = -1.0,Xmax = 1.0,X0 = 0.0;
double Ymin = 0.0,Ymax = 1.0;
double h_rec = 0.1;
double U = 0.149;//箱の運動を決めるパラメータ
double V0 = 1.0;
int N_cell_x = 32,N_cell_y = 12;
double T = 20.0;
double epsilon = 0.000001;
int step_gif = 400;//gifアニメーションのステップ数


struct PARTICLE{
	double x,y,u,v;
	double next;
};

struct CELL{
	int first;
	int last;
};

int intpow(int a,int b);
void status_initialize(struct PARTICLE particle[N]);
void cell_register(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y]);
int set(struct PARTICLE partile[N],int i);
double r_distance(struct PARTICLE particle1,struct PARTICLE particle2);
double v_distance(struct PARTICLE particle1,struct PARTICLE particle2);
double Uniform(void);
double rand_normal( double mu, double sigma );
int getcell_x(double x,double cell_length_x);//ellはcell_lengthのこと
int getcell_y(double y,double cell_length_y);
void Free_evolution(struct PARTICLE *particle,double t);
void G1(struct PARTICLE *particle,int j);
void G2(struct PARTICLE *particle1,struct PARTICLE *particle2);
double T_DWC(struct PARTICLE particle,double t,int j);
double T_DDC(struct PARTICLE particle1,struct PARTICLE particle2,double t);
double Vmax(struct PARTICLE particle[N]);
void DDC(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y]);
void DWC(struct PARTICLE particle[N]);

int main(void){
	FILE *fp_position,*fp_height,*fp_setting;
	char name_position[256];
	sprintf(name_position,"position(N=%d).txt",N);
	if((fp_position = fopen(name_position,"w"))==NULL){
		printf("file_check open error\n");
	}
	if((fp_setting = fopen("setting.txt","w"))==NULL){//gif生成時に必要なパラメータを格納する
		printf("file open error\n");
	}
	if((fp_height = fopen("height_tsd.txt","w"))==NULL){
		printf("file open error\n");
	}
	int i_current,j_current;
	double v_max = 0.0;
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y =(Ymax-Ymin)/(double)N_cell_y;
	double t=0.0,dt=0.01,trec=0.0,dtrec = (double)T/(double)step_gif,t0,t_old,Trec = 0.2*T;
	double t_cell=0.0,t_cell_old=0.0;
	double height;
	//gif生成に必要な情報を出力
	fprintf(fp_setting,"Xmin = %lf\nXmax = %lf\nYmin = %lf\nYmax = %lf\n\n",Xmin,Xmax,Ymin,Ymax);
	fprintf(fp_setting,"n1 = %d\ndt = %lf\n",step_gif,dtrec);
	fprintf(fp_setting,"file = \"position(N=%d)\"\n",N);

	//srand((unsigned) time(NULL));
	struct PARTICLE particle[N];
	struct CELL cell[N_cell_x+1][N_cell_y+1];
	status_initialize(particle);//位置や速度の初期化
	cell_register(particle,cell);//粒子をセルに登録する,nextofの初期化
	
	printf("set up ok\n");
	trec += dtrec;
	
	
	while(t <= T+epsilon){
		for(int i=0;i<N;i++){
			Free_evolution(&particle[i],dt);
		}
		DDC(particle,cell);
		DWC(particle);
		
		cell_register(particle,cell);
		if((t > trec)&&(t <= T)){
			printf("t = %lf\n",t);
			height = 0.0;
			for(int i=0;i<N;i++){
				fprintf(fp_position,"%lf %lf\n",particle[i].x,particle[i].y);
				height += particle[i].y/(double)N;
			}
			fprintf(fp_position,"\n\n");
			fprintf(fp_height,"%lf %lf\n",t,height);
			trec += dtrec;
		}
		t += dt;
	}
	fclose(fp_position);
	fclose(fp_setting);
	fclose(fp_height);
	system("gnuplot gif.txt");
	return 0;
}

int intpow(int a,int b){
	return (int)pow(a,b);
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
		d = r_distance(particle[i],particle[j]);
		if(d <= 2.0*a){
			r = 0;
			break;
		}
	}
	return r;
}

double r_distance(struct PARTICLE particle1,struct PARTICLE particle2){
	double d;
	d = sqrt(pow(particle1.x-particle2.x,2.0)+pow(particle1.y-particle2.y,2.0));
	return d;
}

double v_distance(struct PARTICLE particle1,struct PARTICLE particle2){
	double d;
	d = sqrt(pow(particle1.u-particle2.u,2.0)+pow(particle1.v-particle2.v,2.0));
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
		//printf("x = %lf is out of range\n",x);
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



void Free_evolution(struct PARTICLE *particle,double t){
	particle->x += (particle->u)*t;
	particle->y += (particle->v)*t-0.5*g*t*t;
	particle->v += -g*t;
}

void G1(struct PARTICLE *particle,int j){
	double temp;
	if((j == -1) || (j == -2)){//collision with R or L wall
		particle->u = -e_wall*particle->u;
		if(j == -1){
			particle->x = Xmax-a-epsilon;
		}else{
			particle->x = Xmin+a+epsilon;
		}
	}else if(j == -3){//collision with Bottom wall
		particle->v = (1+e_wall)*U-e_wall*particle->v;
		particle->y = Ymin+a+epsilon;
	}
}


void G2(struct PARTICLE *particle1,struct PARTICLE *particle2){
	double d,Xtemp1,Xtemp2,Ytemp1,Ytemp2,Utemp1,Utemp2,Vtemp1,Vtemp2,Cx,Cy,d_ij;
	d = r_distance(*particle1,*particle2);
	d_ij = 2.0*a-d;
	if(d < 2.0*a){
		Xtemp1 = particle1->x;
		Xtemp2 = particle2->x;
		Ytemp1 = particle1->y;
		Ytemp2 = particle2->y;
		particle1->x += -0.5*d_ij*(Xtemp2-Xtemp1)/d;
		particle2->x += 0.5*d_ij*(Xtemp2-Xtemp1)/d;
		particle1->y += -0.5*d_ij*(Ytemp2-Ytemp1)/d;
		particle2->y += 0.5*d_ij*(Ytemp2-Ytemp1)/d;
		Utemp1 = particle1->u;
		Vtemp1 = particle1->v;
		Utemp2 = particle2->u;
		Vtemp2 = particle2->v;
		Cx = (Xtemp1-Xtemp2)/d;
		Cy = (Ytemp1-Ytemp2)/d;
		particle1->u = 0.5*(1+e)*((Utemp2-Utemp1)*Cx+(Vtemp2-Vtemp1)*Cy)*Cx+Utemp1;
		particle1->v = 0.5*(1+e)*((Utemp2-Utemp1)*Cx+(Vtemp2-Vtemp1)*Cy)*Cy+Vtemp1;
		particle2->u = 0.5*(1+e)*((Utemp1-Utemp2)*Cx+(Vtemp1-Vtemp2)*Cy)*Cx+Utemp2;
		particle2->v = 0.5*(1+e)*((Utemp1-Utemp2)*Cx+(Vtemp1-Vtemp2)*Cy)*Cy+Vtemp2;
	}
}



double Vmax(struct PARTICLE particle[N]){
	double v_max = 0.0;
	for(int i=0;i<N;i++){
		if(v_max*v_max < particle[i].u*particle[i].u+particle[i].v*particle[i].v){
			v_max = sqrt(particle[i].u*particle[i].u+particle[i].v*particle[i].v);
		}
	}
	return v_max;
}


void DDC(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y]){
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	int cell_x,cell_y,j;
	for(int i=0;i<N;i++){
		cell_x = getcell_x(particle[i].x,cell_length_x);
		cell_y = getcell_y(particle[i].y,cell_length_y);
		for(int c1=-1;c1<1;c1++){
			for(int c2=-1;c2<=1;c2++){
				if((((cell_x+c1 >= 1) && (cell_x+c1 <= N_cell_x)) && (cell_y+c2 >= 1)) && (cell_y+c2 <= N_cell_y)){
					j = cell[cell_x+c1][cell_y+c2].first;
					while(j >= 0){
						if(j != i){
							G2(&particle[i],&particle[j]);
						}
						j = particle[j].next;
					}
				}
			}
		}
	}
}

void DWC(struct PARTICLE particle[N]){
	int j;
	for(int i=0;i<N;i++){
		j = 0;
		if(particle[i].x > Xmax-a){
			j = -1;
		}
		if(particle[i].x < -Xmax+a){
			j = -2;
		}
		if(particle[i].y < Ymin+a){
			j = -3;
		}
		G1(&particle[i],j);
	}
}