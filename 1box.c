#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 300//���q��
#define n ((int)ceil(log2(N)))//���S���t�؂̐[��
#define q ((int)pow(2,n)-N)
#define p ((N-q)/2)
double a = 0.02;//���q�̔��a
double e = 0.95,e_wall = 1.0;//���ꂼ�ꗱ�q���m�A�ǂƂ̔����W��
double g = 1.0;//�K�i�����ꂽ�d�͉����x
double Xmin = -1.0,Xmax = 1.0;//���E�̕ǂ̈ʒu
double Ymin = 0.0,Ymax = 1.0;//��ʂƃZ���̍ō��_�̈ʒu
double U = 0.149;//���ʂ̐U�����鑬�x
double V0 = 1.0;//���������ł̑��x���z�̕W���΍�
int N_cell_x = 12,N_cell_y = 8;//x,y�����̃Z���̕�����
double T = 20.0;//�V�~�����[�V�����I������
double epsilon = 0.000001;//���l�덷�̉e�����������߂ɓ���Ă���
int step_gif = 200;//gif�A�j���[�V�����̃X�e�b�v��


struct NODE{//���S���t�؂̃m�[�h�̍\����
	int number;//�Ή����闱�q�ԍ�
	double time;//�Ή����闱�q�̗\�����ꂽ�ŒZ�Փˎ���
	struct NODE *left;//��ƍ��E���Ȃ�
	struct NODE *right;
	struct NODE *parent;
};

struct EVENT{//���闱�q�̃C�x���g�̏ڍ�(�Փˎ����E����)���L�^
	double time;
	int number_particle;
};

struct PARTICLE{//���q�Ɋւ�������܂Ƃ߂�
	double x,y,u,v;//�ʒu�E���x
	double tau;//���q�̌ŗL���Ԃ��L�^�CDelayed State Algorithm(DSA)�ɂ�鍂�����̂��߂ɕK�v
	double next;//�����̎��ɓ����Z���ɓ��������q�̔ԍ�
	struct EVENT event;//���ɗ\�肳���C�x���g
};

struct CELL{
	int first;//�Z���ւ̓o�^���ɍł������o�^���ꂽ���q
	int last;
};

int intpow(int a,int b);
struct EVENT Predictions(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],double t,int i);
void CBT_build(struct NODE *node[n+1][2*p+2*q],struct PARTICLE particle[N]);
void CBT_update(struct NODE *entry[n+1][2*p+2*q],double time_new,int i_new);
void status_initialize(struct PARTICLE particle[N]);
void cell_register(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y]);
int set(struct PARTICLE partile[N],int i);
double r_distance(struct PARTICLE particle1,struct PARTICLE particle2);
double v_distance(struct PARTICLE particle1,struct PARTICLE particle2);
double Uniform(void);
double rand_normal( double mu, double sigma );
int getcell_x(double x,double cell_length_x);//ell��cell_length�̂���
int getcell_y(double y,double cell_length_y);
void Free_evolution(struct PARTICLE *particle,double t);
void G1(struct PARTICLE *particle,int j);
void G2(struct PARTICLE *particle1,struct PARTICLE *particle2);
double T_DWC(struct PARTICLE particle,double t,int j);
double T_DDC(struct PARTICLE particle1,struct PARTICLE particle2,double t);
double NextEvent(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],int i_current,int j_current);
double t_cell_update(struct PARTICLE particle,int j_current,double t_cell_old,double *v_max);
double Vmax(struct PARTICLE particle[N]);
void MaskUpdate(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],int i_current,double t);
double EEPGM(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],double t,double *v_max);


int main(void){
	FILE *fp_position,*fp_height,*fp_setting;//�t�@�C���̐���
	char name_position[256];
	sprintf(name_position,"position(N=%d).txt",N);//��莞�����Ƃɗ��q�̈ʒu��ۑ�����t�@�C��
	if((fp_position = fopen(name_position,"w"))==NULL){
		printf("file_check open error\n");
	}
	if((fp_height = fopen("height.txt","w"))==NULL){//���q�̍����̕��ϒl���L�^
		printf("file open error\n");
	}
	if((fp_setting = fopen("setting.txt","w"))==NULL){//gif�������ɕK�v�ȃp�����[�^���i�[����
		printf("file open error\n");
	}
	int i_current,j_current;//���ݒ��ڂ��Ă��闱�q�̃y�A,j_current<0:��,j_current>=0:���q
	double v_max = 0.0;//�ő呬�x��ۑ��A�Z���̍X�V�̂��߂ɕK�v
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y =(Ymax-Ymin)/(double)N_cell_y;//1�̃Z���̍����Ɖ���
	double t=0.0,dt=0.01,trec=0.0,dtrec = (double)T/(double)step_gif;//dtrec��gif�����̎��ԊԊu
	double t_cell=0.0,t_cell_old=0.0;//�Z���̍X�V����
	double height;
	//gif�����ɕK�v�ȏ����o��
	fprintf(fp_setting,"Xmin = %lf\nXmax = %lf\nYmin = %lf\nYmax = %lf\n",Xmin,Xmax,Ymin,Ymax);
	fprintf(fp_setting,"n1 = %d\ndt = %lf\n",step_gif,dtrec);
	fprintf(fp_setting,"file = \"position(N=%d)\"\n",N);
	
	srand((unsigned) time(NULL));
	struct PARTICLE particle[N];//���q��\���\����
	struct CELL cell[N_cell_x][N_cell_y];
	status_initialize(particle);//�ʒu�⑬�x�̏�����
	v_max = Vmax(particle);//���q�̍ő呬�x
	t_cell = (cell_length_y-2.0*a)/(2.0*v_max);//���̎��Ԃ܂łɃ}�X�N�O����̏Փ˂͂��肦�Ȃ�
	cell_register(particle,cell);//���q���Z���ɓo�^����,nextof�̏�����
	
	struct NODE *node[n+1][2*p+2*q];//���S���t��(���邢�̓g�[�i�����g)��\���\����
	//node�̎��̉��A����������Ƃ������ɂ�����
	for(int i=0;i<=n;i++){
		for(int j=0;j<2*p+2*q;j++){
			node[i][j] = (struct NODE *)malloc(sizeof(struct NODE));//�̈���m�ۂ���
		}
	}
	
	for(int i=0;i<N;i++){
		particle[i].event = Predictions(particle,cell,t,i);//���ꂼ��̗��q�̎��̃C�x���g��\��
	}
	CBT_build(node,particle);//Complete Binary Tree��g�ݗ��Ă�
	printf("set up ok\n");
	
	while(t <= T){
		//NEXT EVENT�̌���
		i_current = node[0][0]->number;//�����̃m�[�h�͍ŒZ�̎��ԂŏՓ˂��闱�q������
		j_current = particle[i_current].event.number_particle;//i_current�̏Փˑ���(����͕ǂ̉\��������)
		t = NextEvent(particle,cell,node,i_current,j_current);//NEXT EVENT��������t��particle,cell,node���X�V
		t_cell = t_cell_update(particle[i_current],j_current,t_cell_old,&v_max);//t_cell��v_max�̍X�V
		
		//i_current,j_current�ƏՓ˂���\�肾�������q�������ꍇ�͂��̗��q��event��invalid�ɂȂ��Ă��܂��̂ŐV����event�����
		//���̂悤�ȗ��q�͓����}�X�N���ɂ������݂��Ȃ��͂��Ȃ̂ł��̒��ŒT��
		MaskUpdate(particle,cell,node,i_current,t);//i_current�̎���̗��q��invalid�Ȃ��̂�����΃A�b�v�f�[�g
		if(j_current >= 0){//j�ɂ��Ă����l
			MaskUpdate(particle,cell,node,j_current,t);
		}
		
		//EEPGM �}�X�N�O�̗��q�Ƃ��Փ˂���\����������̂œo�^������
		if(t >= t_cell){
			t_cell_old = t;
			t_cell = EEPGM(particle,cell,node,t,&v_max);
			//���ɗ��q���߂荞��ł����炱�̃G���[��������
//			for(int i=0;i<N;i++){
//				if(particle[i].y < Ymin+a-epsilon){
//					printf("i=%d:error\n",i);
//					printf("%lf %lf %lf %lf\n",particle[i].x,particle[i].y,particle[i].u,particle[i].v);
//					printf("%lf %d\n",particle[i].event.time,particle[i].event.number_particle);
//					G1(&particle[i],-3);
//					particle[i].event = Predictions(particle,cell,t,i);
//					CBT_update(node,particle[i].event.time,i);
//					MaskUpdate(particle,cell,node,i,t);
//				}
//			}
			
		}
		//���q�̈ʒu�̏o��
		if((t > trec)&&(t < T)){
			t_cell_old = t;
			t_cell = EEPGM(particle,cell,node,t,&v_max);
			printf("t = %lf, v_max = %lf\n",t,v_max);
			height = 0.0;
			for(int i=0;i<N;i++){
				fprintf(fp_position,"%lf %lf\n",particle[i].x,particle[i].y);
				height += particle[i].y/(double)N;
				//���q���m���߂荞��ł���󋵂ɂȂ��Ă��Ȃ����m�F
//				for(int j = 0;j<i;j++){
//					if(r_distance(particle[i],particle[j]) < 2.0*a-epsilon){
//						printf("%d %d is too close!!\n",i,j);
//						printf("distance = %lf\n",r_distance(particle[i],particle[j]));
//						printf("particle1:%lf %lf %lf %lf\n",particle[i].x,particle[i].y,particle[i].u,particle[i].v);
//						printf("particle2:%lf %lf %lf %lf\n",particle[j].x,particle[j].y,particle[j].u,particle[j].v);
//					    break;
//					}
//				}
			}
			fprintf(fp_position,"\n\n");//gif�����̂��߂ɕK�v�ȓ�̋�
			fprintf(fp_height,"%lf %lf\n",t,height);
			trec += dtrec;
		}
	}
	
	fclose(fp_position);
	fclose(fp_height);
	fclose(fp_setting);
	system("gnuplot gif.txt\n");
	return 0;
}

int intpow(int a,int b){//pow()�̐�����,�g�[�i�����g���쐬����Ƃ��ɕK�v
	return (int)pow(a,b);
}

struct EVENT Predictions(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],double t,int i){//���qi�ɑ΂��čŒZ�Ő�����C�x���g��\�����ďo��
	double t_min = 2.0*T,t_temp;
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	int j_col,j;
	struct PARTICLE particle_j;
	struct EVENT L;
	
	for(j=-3;j<0;j++){//�ǂƂ̏Փˎ��Ԃ��m�F
		t_temp = T_DWC(particle[i],particle[i].tau,j);
		if((t_temp > t) & ( t_temp < t_min)){
			t_min = t_temp;
			j_col = j;
		}
	}
	
	int cell_x = getcell_x(particle[i].x,cell_length_x),cell_y = getcell_y(particle[i].y,cell_length_y);
	for(int c1=-1;c1<=1;c1++){//�ߕӂ̗��q�Ƃ̏Փˎ��Ԃ��m�F
		for(int c2=-1;c2<=1;c2++){
			if((((cell_x+c1 >= 0) && (cell_x+c1 < N_cell_x)) && (cell_y+c2 >= 0)) && (cell_y+c2 < N_cell_y)){
				j = cell[cell_x+c1][cell_y+c2].first;//�����N���X�g�\���𗘗p���ė��q��{��
				while(j >= 0){
					particle_j = particle[j];
					Free_evolution(&particle_j,particle[i].tau-particle[j].tau);//���藱�qj�̎��Ԃ�i�Ƃ��낦��
					t_temp = T_DDC(particle[i],particle_j,particle[i].tau);//�Փ˂ɂ����鎞�Ԃ��v�Z
					if((t_temp > t) && (t_temp < t_min)){//���ݎ������x���Ct_min�����������Ԃł����t_min�̍X�V
						t_min = t_temp;
						j_col = j;//���̂Ƃ��̑���j���L�^
					}
					j = particle[j].next;
				}
			}
		}
	}
	L.time = t_min;
	L.number_particle = j_col;
	return L;//���ԂƑ���̏����o��
}
void CBT_build(struct NODE *node[n+1][2*p+2*q],struct PARTICLE particle[N]){//CBT�����
	int i,n_index;
	//initialization for bottom nodes
	for(i=0;i<2*p+2*q;i++){
		if(i < 2*p){
			node[n][i]->time = particle[i].event.time;
			node[n][i]->number = i;
		}else{
			node[n][i]->time = particle[2*p+(i-2*p)/2].event.time;
			node[n][i]->number = 2*p+(i-2*p)/2;
		}
	}
	//tournament
	for(n_index=n-1;n_index>=0;n_index--){
		for(i=0;i<=intpow(2,n_index)-1;i++){
			node[n_index][i]->left = node[n_index+1][2*i];
			node[n_index+1][2*i]->parent = node[n_index][i];
			node[n_index][i]->right = node[n_index+1][2*i+1];
			node[n_index+1][2*i+1]->parent = node[n_index][i];
			if(node[n_index+1][2*i]->time <= node[n_index+1][2*i+1]->time){
				node[n_index][i]->time = node[n_index+1][2*i]->time;
				node[n_index][i]->number = node[n_index+1][2*i]->number;
			}else{
				node[n_index][i]->time = node[n_index+1][2*i+1]->time;
				node[n_index][i]->number = node[n_index+1][2*i+1]->number;
			}
		}
	}
}

void CBT_update(struct NODE *entry[n+1][2*p+2*q],double time_new,int i_new){//i_new�̏����X�V����
	struct NODE *entry_now,hoge_now;
	if(i_new < 2*p){
		entry[n][i_new]->time = time_new;
		entry_now = entry[n][i_new];
	}else{
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

void status_initialize(struct PARTICLE particle[N]){//���q�̏������������߂�
	double prob;
	int i;
	for(i=0;i<N;i++){
		prob = Uniform();
		particle[i].x = (Xmin+a)*prob+(Xmax-a)*(1-prob);
		prob = Uniform();
		particle[i].y = (Ymin+a)*prob+(0.5*Ymax-a)*(1-prob);
		while(set(particle,i) == 0){//�����d�Ȃ��Ă��闱�q���������Ƃ��͏d�Ȃ肪�Ȃ��Ȃ�܂œo�^������
			prob = Uniform();
			particle[i].x = (Xmin+a)*prob+(Xmax-a)*(1-prob);
			prob = Uniform();
			particle[i].y = (Ymin+a)*prob+(Ymax-a)*(1-prob);
		}
		particle[i].u = rand_normal(0.0,V0);
		particle[i].v = rand_normal(0.0,V0);
		particle[i].next = -1;
		particle[i].tau = 0.0;
		particle[i].event.time = 2.0*T;
		particle[i].event.number_particle = -1;
	}
}

void cell_register(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y]){//���ׂĂ̗��q���Z���ɓo�^������
	int i,cell_x,cell_y,lastPrev;
	//initialize particle.next and cell
	for(i=0;i<N;i++){
		particle[i].next = -1;
	}
	for(cell_x = 0;cell_x < N_cell_x;cell_x++){
		for(cell_y = 0;cell_y < N_cell_y;cell_y++){
			cell[cell_x][cell_y].first = -1;
			cell[cell_x][cell_y].last = -1;
		}
	}
	//�����N���X�g�\���̍쐬
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


int set(struct PARTICLE particle[N],int i){//�d�Ȃ�Ȃ����q��u�����Ƃɐ������Ă����1,���s���Ă����0��Ԃ�
	int j,r=1;
	double d;
	
	if(fabs(particle[i].x) < a){
		r = 0;
	}
	for(j=0;j<=i-1;j++){
		d = r_distance(particle[i],particle[j]);
		if(d <= 2.0*a){
			r = 0;
			break;
		}
	}
	return r;
}

double r_distance(struct PARTICLE particle1,struct PARTICLE particle2){//2�̗��q�̋������v�Z
	double d;
	d = sqrt(pow(particle1.x-particle2.x,2.0)+pow(particle1.y-particle2.y,2.0));
	return d;
}

double v_distance(struct PARTICLE particle1,struct PARTICLE particle2){//2�̗��q�̑��x�x�N�g���̍��̑傫�����v�Z
	double d;
	d = sqrt(pow(particle1.u-particle2.u,2.0)+pow(particle1.v-particle2.v,2.0));
	return d;
}

double Uniform(void){//0����1�̈�l�����𐶐�
	return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

double rand_normal( double mu, double sigma ){//����mu,�W���΍�sigma�̐��K���z�𐶐�
	double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
	return mu + sigma*z;
}

int getcell_x(double x,double cell_length_x){//x�Ƃ����ʒu�̃Z���̔ԍ���Ԃ�
	if((x < Xmin+a)||(Xmax-a < x)){
		printf("x is out of range\n");
	}
	return (int)((x-Xmin)/cell_length_x);
}

int getcell_y(double y,double cell_length_y){//y�Ƃ����ʒu�̃Z���̔ԍ���Ԃ�
	if(y < Ymin){
		printf("error:y<0(%lf)\n",y);
		return 0;
	}else if(y>Ymax){
		return N_cell_y-1;//Ymax���������ʒu�̗��q�͈�ԍ����Z���ɓo�^
	}else{
		return (int)(y/cell_length_y);
	}
}

void Free_evolution(struct PARTICLE *particle,double t){//���闱�q������t�������Ԕ��W������
	particle->x += (particle->u)*t;
	particle->y += (particle->v)*t-0.5*g*t*t;
	particle->v += -g*t;
	particle->tau += t;//�ŗL���Ԃ̍X�V���K�v�Ȃ��Ƃɒ���
}

void G1(struct PARTICLE *particle,int j){//���q�ƕǂ̏Փˏ������s��
	double temp;
	if((j == -1) || (j == -2)){//collision with R or L wall
		particle->u = -e_wall*particle->u;
		if(j == -1){
			particle->x = Xmax-a-epsilon;//����epsilon������getcell_x�̂Ƃ��ȂǂɕK�v�ɂȂ�
		}else{
			particle->x = Xmin+a+epsilon;
		}
	}else if(j == -3){//collision with Bottom wall
		particle->v = (1+e_wall)*U-e_wall*particle->v;
		particle->y = Ymin+a+epsilon;
	}
}


void G2(struct PARTICLE *particle1,struct PARTICLE *particle2){//���q���m�̏Փˏ���
	double d,Xtemp,Ytemp,Utemp1,Utemp2,Vtemp1,Vtemp2,Cx,Cy;
	d = r_distance(*particle1,*particle2);
	Utemp1 = particle1->u;
	Vtemp1 = particle1->v;
	Utemp2 = particle2->u;
	Vtemp2 = particle2->v;
	Cx = (particle1->x-particle2->x)/d;
	Cy = (particle1->y-particle2->y)/d;
	particle1->u = 0.5*(1+e)*((Utemp2-Utemp1)*Cx+(Vtemp2-Vtemp1)*Cy)*Cx+Utemp1;
	particle1->v = 0.5*(1+e)*((Utemp2-Utemp1)*Cx+(Vtemp2-Vtemp1)*Cy)*Cy+Vtemp1;
	particle2->u = 0.5*(1+e)*((Utemp1-Utemp2)*Cx+(Vtemp1-Vtemp2)*Cy)*Cx+Utemp2;
	particle2->v = 0.5*(1+e)*((Utemp1-Utemp2)*Cx+(Vtemp1-Vtemp2)*Cy)*Cy+Vtemp2;
}

double T_DDC(struct PARTICLE particle1,struct PARTICLE particle2,double t){//���q���m�̏Փ�(DDC)�̎��Ԃ��v�Z
	double r_relative,v_relative,b,hoge;
	double tau = t;
	double x1 = particle1.x ,x2 = particle2.x ,y1 = particle1.y , y2 = particle2.y;
	double u1 = particle1.u ,u2 = particle2.u ,v1 = particle1.v , v2 = particle2.v;
	r_relative = r_distance(particle1,particle2);
	v_relative = v_distance(particle1,particle2);
	b = (x1-x2)*(u1-u2)+(y1-y2)*(v1-v2);
	hoge = b*b-v_relative*v_relative*(r_relative*r_relative-4.0*a*a);
	if(hoge > 0.0){
		tau += -(b+sqrt(hoge))/(v_relative*v_relative);
	}else{
		tau += T;
	}
	return tau;
}

double T_DWC(struct PARTICLE particle,double t,int j){//���q�ƕǂ̏Փ˂̎��Ԃ��v�Z
	double tau = t;
	if(j==-1){//collision with RIGHT wall(-1)
		if(particle.u>0.0){
			tau += (Xmax-a-particle.x)/particle.u;
		}else{
			tau += 2.0*T;
		}
	}else if(j==-2){//collision with LEFT wall(-2)
		if(particle.u<0.0){
			tau += (Xmin+a-particle.x)/particle.u;
		}else{
			tau += 2.0*T;
		}
	}else if(j==-3){//collision with BOTTOM wall(-3)
		tau += (particle.v+sqrt(particle.v*particle.v+2*g*(particle.y-a)))/g;
	}
	if(tau < t){
		return 2.0*T;
	}else{
		return tau;
	}
}
double NextEvent(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],int i_current,int j_current){//i_current��j_current�̃C�x���g�����ۂɍs���C���̎�����Ԃ��֐�
	double t = particle[i_current].event.time;
	Free_evolution(&particle[i_current],t-particle[i_current].tau);//i_current�̎��Ԕ��W
	if(j_current >= 0){//Disk Disk Collision
		Free_evolution(&particle[j_current],t-particle[j_current].tau);//j_current�̎��Ԕ��W
		G2(&particle[i_current],&particle[j_current]);//���q���m�̏Փˏ���
		if(r_distance(particle[i_current],particle[j_current]) < 2.0*a-epsilon){
			printf("%d %d is too close!!\n",i_current,j_current);
			printf("distance = %lf\n",r_distance(particle[i_current],particle[j_current]));
		}
	}
	if(j_current < 0){//Disk Wall Collision
		G1(&particle[i_current],j_current);//�ǂƂ̏Փˏ���
	}
	particle[i_current].event = Predictions(particle,cell,t,i_current);//i_current�̃C�x���g�X�V
	CBT_update(node,particle[i_current].event.time,i_current);//i_current��node�A�b�v�f�[�g
	if(j_current >= 0){//j_current�ɂ��Ă����l
		particle[j_current].event = Predictions(particle,cell,t,j_current);
		CBT_update(node,particle[j_current].event.time,j_current);
	}
	return t;
}

double t_cell_update(struct PARTICLE particle,int j_current,double t_cell_old,double *v_max){//�Z���̍X�V�������v�Z
	double t_cell,dt_cell;
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	if(j_current == -3){
		if(*v_max*(*v_max) < pow(particle.u,2.0)+pow(particle.v,2.0)){
			*v_max = sqrt(pow(particle.u,2.0)+pow(particle.v,2.0));
		}
	}
	dt_cell = (cell_length_y-2.0*a)/(2.0*(*v_max));
	t_cell = t_cell_old+dt_cell;
	return t_cell;
}

double Vmax(struct PARTICLE particle[N]){//�ő呬�x���v�Z
	double v_max = 0.0;
	for(int i=0;i<N;i++){
		if(v_max*v_max < particle[i].u*particle[i].u+particle[i].v*particle[i].v){
			v_max = sqrt(particle[i].u*particle[i].u+particle[i].v*particle[i].v);
		}
	}
	return v_max;
}

void MaskUpdate(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],int i_current,double t){//�����}�X�N�Ɋ܂܂�闱�q�̃C�x���g���X�V
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	int cell_x = getcell_x(particle[i_current].x,cell_length_x) , cell_y = getcell_y(particle[i_current].y,cell_length_y),j;
	for(int c1=-1;c1<=1;c1++){
		for(int c2=-1;c2<=1;c2++){
			if((((cell_x+c1 >= 0) && (cell_x+c1 < N_cell_x)) && (cell_y+c2 >= 0)) && (cell_y+c2 < N_cell_y)){
				j = cell[cell_x+c1][cell_y+c2].first;
				while(j >= 0){
					if(particle[j].event.number_particle == i_current){
						particle[j].event = Predictions(particle,cell,t,j);
						CBT_update(node,particle[j].event.time,j);
					}
					j = particle[j].next;
				}
			}
		}
	}
}

double EEPGM(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],double t,double *v_max){//���ׂĂ̗��q�Ɋւ��Ď��Ԕ��W�������̂��Z���ɓo�^������
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	double dt_cell,t_cell;

	for(int i=0;i<N;i++){//���ݎ����܂Ŏ��Ԕ��W
		Free_evolution(&particle[i],t-particle[i].tau);
	}

	cell_register(particle,cell);//�S���q���Z���ɓo�^������
	for(int i=0;i<N;i++){//�S���q�ɂ���event���v�Z������
		particle[i].event = Predictions(particle,cell,t,i);
	}
	CBT_build(node,particle);//CBT���ŏ�����\��
	*v_max = Vmax(particle);
	dt_cell = (cell_length_y-2.0*a)/(2.0*(*v_max));
	t_cell = t+dt_cell;
	return t_cell;
}