#include "navi.h"
#include <stdlib.h>
#include<Arduino.h>
#include<TimerOne.h>


extern long volatile  milliss,t;
unsigned long long int delta,t0;

unsigned int prevMillis=0;
unsigned int currentMillis=0;

//Motors realted variables
extern int PWM_R,PWM_L;
extern int PWM_Max;
extern int PWM_R_Min,PWM_L_Min;
extern int PWM_R_Min_Rot,PWM_L_Min_Rot;
extern int PWM_R_sign_counter,PWM_L_sign_counter;

//Odometry related variables
extern volatile float total_right,total_left, total_centre;
extern volatile float current_x,current_y,current_phi_deg,current_phi_rad;
extern float ref_x, ref_y;
extern  float spacing_encoder,spacing_wheel/*,dec*/;
extern volatile double right_speed, left_speed;

//Robot Navi Related Variables
float goal_distance, target_angle;
float accel_dist, decel_dist;
int PWM_LB,PWM_RB;
int coef_correct_angle=100; //50 Correction rotate
int right_correction=0,left_correction=0;

//Speed Regulation
float speed_ref, ramp = 500,rampR=100,rampC=500;//max 1500
int sens;
double right_error=0,i_right_error=0;
double left_error=0,i_left_error=0;
float kp = 15, ki = 0.8;

//Move
int coef_correct_dist = 20;

//Trajectory
float target_x, target_y;
float target_x_prime, target_y_prime;
float right_target_speed, left_target_speed;

// curv
float remain_distC=0,goalC=0;
float speed_refR=0,speed_refL=0,prev_speed_refR=0,prev_speed_refL=0,new_speed_refR=0,new_speed_refL=0,speedC=0,speed_refC=0;
float kpL = 15.5, kiL = 1.75;
float kpR = 15.5, kiR = 1.75;
float corde=0,tetaC=0,phi_prim=0,corde_angle=0,Xc=0,Yc=0,Rayon=0,phi_target_rad=0,sens_de_mouvement=0;
float Distance_empietement=50;
float** matrix;
int flag1=0;
int flag2=0;



float constrainn( float x,float min,float Max)
{
	if (x<min) return min;
	if (x>Max) return Max;
	else return(x);
}

void initt (void)
{
	total_right=0;
	total_left=0;
	total_centre=0;
	PWM_L_sign_counter=0;
	PWM_R_sign_counter=0;
	i_right_error = 0;
	i_left_error = 0;
	right_error=0;
	left_error=0;
}
// intiailisation de la fonction multi_curv (pour accumuler l'Integrale)
void initMultiCurv (void)
{
	total_right=0;
	total_left=0;
	total_centre=0;
	PWM_L_sign_counter=0;
	PWM_R_sign_counter=0;
	right_error=0;
	left_error=0;
}

//allocation d une matrice [taille][3] : x/y/phi_target pour (taille) points
void allocation (int taille ) {	
	matrix=(float**) malloc (taille*(sizeof(float*)));
	for(int i=0;i<taille;i++)
		matrix[i]=(float*) malloc(3*sizeof(float));
}


void move_distance(float distance,float speed)
{
	initt();
	//Set accel/decel distance
	if (fabs(distance) < (speed*speed/ramp))
	{
		accel_dist = fabs(distance)/2;
		decel_dist = fabs(distance)/2;
		speed = sqrt (2*ramp*accel_dist);
	}
	else
	{
		accel_dist = (float)0.5*speed*speed/ramp;
		decel_dist = (float)0.5*speed*speed/ramp;
	}
	while(fabs((total_right+total_left)/2-distance)>1)
	{
		t0=t;
		//Accel/Decel Speed Set
		if (((total_right+total_left)/2 -distance)<0)
			sens = 1;
		else
			sens = -1;
		if (fabs((total_right+total_left)/2) < accel_dist)
			speed_ref = sens*50+sens*(constrainn(sqrt (ramp*fabs(total_right+total_left))-50,0,1000));
	
		else if (fabs((total_right+total_left)/2 -distance) < decel_dist)
			speed_ref = sens*10+sens*constrainn((sqrt(2*ramp*fabs((total_right+total_left)/2 -distance))-10),0,1000);//fabs((total_right+total_left)/2 -distance)
		else
			speed_ref = sens*speed;
		//Right wheel regulation
		right_error = speed_ref - right_speed;
		i_right_error += right_error;
		PWM_RB = kp * right_error + ki * i_right_error;
		if (PWM_RB>PWM_Max) PWM_RB = PWM_Max;
		if (PWM_RB<-PWM_Max) PWM_RB = -PWM_Max;
		//Left wheel regulation
		left_error = speed_ref - left_speed;
		i_left_error += left_error;
		PWM_LB = kp * left_error + ki * i_left_error;
		if (PWM_LB>PWM_Max) PWM_LB = PWM_Max;
		if (PWM_LB<-PWM_Max) PWM_LB = -PWM_Max;
		//Orientation Correction²
		left_correction = coef_correct_dist * (total_right-total_left);
		right_correction = - left_correction;
		PWM_R = PWM_RB + right_correction;
		PWM_L = PWM_LB + left_correction;
		//Execution
		run_motors();
		do delta=t-t0;
		while (delta<T);
	}
	stop_motors();
}

//calcul des donnees de la fonction robot locate curv
void Calcul_Curv(float x, float y,float phi_target) {

	corde=0;
	corde_angle=0;
	Rayon=0;
	Xc=0;
	Yc=0;
	tetaC=0;
	phi_target_rad=(double)PI*((phi_target)/180);
	
	//Calcul de la corde 
	corde=sqrtf((pow((x-current_x),2)) + (pow((y-current_y),2)));
	
	//Calcul des coordonnées du centre de la courbure
		
	if ((phi_target==0) || (phi_target==180) || (phi_target==-180)){ 
		if ((y==0) && (current_y==0))
				current_y=0.000000001; //Résolution du probleme de division par zero 
		Xc=x;
		Yc=(float)(((current_x-x)/(y-current_y))*Xc)+(0.5f*(y+current_y))-((pow(current_x,2)-pow(x,2))/(2*(y-current_y)));
	}
	else{
		Xc=(float)((cosf(phi_target_rad)/sinf(phi_target_rad))*x+(((current_x*current_x)-(x*x))/(2*(y-current_y)))-(0.5f*(current_y-y)))/((cosf(phi_target_rad)/sinf(phi_target_rad))+((current_x-x)/(y-current_y)));
		Yc=(float)((cosf(phi_target_rad)/sinf(phi_target_rad))*(x-Xc))+y;
	}
	
	//Calcul du rayon de la courbure
	Rayon=sqrtf((pow((x-Xc),2)) + (pow((y-Yc),2)));
	
	//Calcul de l'angle de la corde
	sens = (asinf((y-current_y)/sqrtf((current_x-x)*(current_x-x)+(current_y-y)*(current_y-y)))>0)? 1 : -1;
	corde_angle = sens * rad_to_deg(acosf((x-current_x)/sqrtf((current_x-x)*(current_x-x)+(current_y-y)*(current_y-y))));	
	
	//Calcul de la l'angle parcouru tout au long de la courbure [0,180] en valeur absolue
	tetaC=rad_to_deg (2*asin((double)(corde/(2*Rayon))));
	
	//Calcul du sens_de_mouvement {recucler ou avancer}
	sens_de_mouvement=1;	
	if((phi_target-corde_angle>=90) ||(phi_target- corde_angle<=-90))
		sens_de_mouvement=-1;
	
	//Calcul du sens de rotation {trigonometique ou horraire}
	if (sens_de_mouvement==1) {
		if((phi_target-corde_angle)<0) 
			tetaC=-1*tetaC; // horraire
	}
	else if (sens_de_mouvement==-1) {
		if((-phi_target+corde_angle)<0) 
			tetaC=-1*tetaC; // trigo
	}
	
	//Calcul de l'orientation initiale avant de commencer la courbure
	phi_prim=phi_target-tetaC;
	sens_de_mouvement=1; //ajoutéé
}

// Calcul des vitesses de chaque roue selon la courbure et le sens de mouvement
void Calcul_speed_Refs(float theta,float R,int sens_de_mouvement, float speed){
		if(sens_de_mouvement==1){
			if (theta>=0) {
				speed_refL=speed*((R-(spacing_wheel*0.5f))/R);
				speed_refR=speed*((R+(spacing_wheel*0.5f))/R);
			}
			else if (theta<0) {
				speed_refR=speed*((R-(spacing_wheel*0.5f))/R);
				speed_refL=speed*((R+(spacing_wheel*0.5f))/R);
				}
		}
		else{
			if (theta<0) {
				speed_refL=speed*((R-(spacing_wheel*0.5f))/R);
				speed_refR=speed*((R+(spacing_wheel*0.5f))/R);
			}
			else if (theta>=0) {
				speed_refR=speed*((R-(spacing_wheel*0.5f))/R);
				speed_refL=speed*((R+(spacing_wheel*0.5f))/R);
			}
		}
}
// robot locate curv
void Robot_locateCurv(float x, float y, float phi_target, float speed) {
	Calcul_Curv(x,y,phi_target);
  orientate(phi_prim,speed);
	curv(Rayon,tetaC,speed);
}
//curv
void curv(float R,float theta ,float speed)
{ initt();	 
	theta=(theta*PI)/180;
	goalC=fabs(R*theta); 
	remain_distC=goalC-fabs(total_centre);
	speedC=speed;
	if (goalC < (speedC*speedC/rampC)) {
		accel_dist = fabs(remain_distC)/2;
		decel_dist = fabs(remain_distC)/2;
		speedC = speedC* accel_dist / ((float)0.5*speedC*speedC/rampC);}
	else{
		accel_dist = (float)0.5*speedC*speedC/rampC;
		decel_dist = (float)0.5*speedC*speedC/((rampC));
	}	
	while((((remain_distC)>=2))||((remain_distC)<=-2)) { 
		t0=t;
		goalC=fabs(R*theta); 
	  remain_distC=goalC-fabs(total_centre);
	  if ((fabs(total_centre)) < accel_dist){
			speed_refC = sens_de_mouvement*50+sens_de_mouvement*(speedC-50)*fabs((total_centre)/(accel_dist));		
		}
		else if (fabs((remain_distC)) < decel_dist)
			speed_refC = sens_de_mouvement*10+sens_de_mouvement*(speedC-10)*fabs((remain_distC)/(decel_dist));
		else
			speed_refC = sens_de_mouvement*speedC;
		Calcul_speed_Refs( theta, R, sens_de_mouvement,speed_refC);
		//Right wheel regulation
		right_error = speed_refR - right_speed;
		i_right_error += right_error;
		PWM_R = kpR * right_error + kiR * i_right_error;
		//Left wheel regulation
		left_error = speed_refL - left_speed;
		i_left_error += left_error;
		PWM_L = kpL * left_error + kiL * i_left_error;
		//Execution
		run_motors();
		do delta=t-t0;
		while (delta<T);
	  }
  	stop_motors();
}
//curv function used in the robot_locate_multi_curv function
 void Multi_Curv(float R,float theta ,float speed,int i,int n){ 
	initMultiCurv(); // integrale accumulée
	theta=(theta*PI)/180;
	goalC=fabs(R*theta); 
	remain_distC=goalC-fabs(total_centre);
	speedC=speed;
	if (goalC < (speedC*speedC/rampC)) {
		accel_dist = fabs(remain_distC)/2;
		decel_dist = fabs(remain_distC)/2;
		speedC = speedC* accel_dist / ((float)0.5*speedC*speedC/rampC);
	}
	else{
		accel_dist = (float)0.5*speedC*speedC/rampC;
		decel_dist = (float)0.5*speedC*speedC/((rampC));
	}	
	
		if (i!=0) accel_dist=0; //Accélération du robot juste au début de la trajectoire
		if (i!=(n-1)) decel_dist=0;// Décèlération du robot seulement a la fin de la trajectoire
	
	while(fabs(total_centre)<fabs(goalC)) {  //(((remain_distC)>=2))||((remain_distC)<=-2) 
		t0=t;
	  remain_distC=goalC-fabs(total_centre);
	  if ((fabs(total_centre)) < accel_dist){
			speed_refC = sens_de_mouvement*50+sens_de_mouvement*(speedC-50)*fabs((total_centre)/(accel_dist));		
		}
		else if (fabs((remain_distC)) < decel_dist)
			speed_refC = sens_de_mouvement*10+sens_de_mouvement*(speedC-10)*fabs((remain_distC)/(decel_dist));
		else
			speed_refC = sens_de_mouvement*speedC;		
		Calcul_speed_Refs(theta,R,sens_de_mouvement,speed_refC);
		//Empietement {just ask us} Pour plus d'informations, contactez Nesrine ou Wassim
		if(i!=(n-1)){
			if ((remain_distC<=Distance_empietement)) {
				if (flag1==0) {
					prev_speed_refL=speed_refL;
					prev_speed_refR=speed_refR;
				
					Calcul_Curv(matrix[i+1][0],matrix[i+1][1],matrix[i+1][2]);
					speedC=speed;
					tetaC=(tetaC*PI)/180;
					if ((fabs(Rayon*tetaC)-Distance_empietement)< (speedC*speedC/rampC)) {
						speedC = speedC* (fabs((Rayon*tetaC)/2)-(0.5f*Distance_empietement)) / ((float)0.5f*speedC*speedC/rampC);
					}
					Calcul_speed_Refs(tetaC,Rayon,sens_de_mouvement,sens_de_mouvement*speedC);
					new_speed_refL=speed_refL;
					new_speed_refR=speed_refR;
					flag1=1;
				}
			speed_refR=(((1/Distance_empietement)*(((-0.5f)*(prev_speed_refR+new_speed_refR))+prev_speed_refR))*(remain_distC))+(0.5f*(prev_speed_refR+new_speed_refR));
			speed_refL=(((1/Distance_empietement)*(((-0.5f)*(prev_speed_refL+new_speed_refL))+prev_speed_refL))*(remain_distC))+(0.5f*(prev_speed_refL+new_speed_refL));
			} 
		}
		if ((total_centre<=Distance_empietement)&&(i!=0)) {
			if (flag2==0) {
					Calcul_Curv(matrix[i][0],matrix[i][1],matrix[i][2]);
					Calcul_speed_Refs(tetaC,Rayon,sens_de_mouvement,speed_refC);
					new_speed_refL=speed_refL;
					new_speed_refR=speed_refR;
					flag2=1;
				}
			speed_refR=(((1/Distance_empietement)*(((-0.5f)*(prev_speed_refR+new_speed_refR))+new_speed_refR))*(total_centre))+(0.5f*(prev_speed_refR+new_speed_refR));
			speed_refL=(((1/Distance_empietement)*(((-0.5f)*(prev_speed_refL+new_speed_refL))+new_speed_refL))*(total_centre))+(0.5f*(prev_speed_refL+new_speed_refL));
		}	
		//Right wheel regulation
		right_error = speed_refR - right_speed;
		i_right_error += right_error;
		PWM_R = kpR * right_error + kiR * i_right_error;
		//Left wheel regulation
		left_error = speed_refL - left_speed;
		i_left_error += left_error;
		PWM_L = kpL * left_error + kiL * i_left_error;
		//Execution
		run_motors();
		do delta=t-t0;
		while (delta<T);
	  }
	flag1=0;
	flag2=0;
	prev_speed_refL=speed_refL;
	prev_speed_refR=speed_refR;
}

// matrix colonne : x/y/phi_target
// in the matrix put the x,y and phi_target you want the robot to go to in order 

 void Robot_Locate_Multi_Curv( float** matrix , int n , int speed) {
 for(int i=0;i<n;i++) {
		Calcul_Curv(matrix[i][0],matrix[i][1],matrix[i][2]);
		if (i==0)
		{
			orientate(phi_prim,speed);
			initt();
		}
			Multi_Curv(Rayon,tetaC,speed,i,n);
 }
	 stop_motors();
}
		
void rotate(float angle, float speed)
{
	initt();
	//Set accel/decel distance
	goal_distance = angle * PI * spacing_encoder/ 180;
	if (fabs(goal_distance) < (2*speed*speed/rampR))
	{
		accel_dist = fabs(goal_distance)/2;
		decel_dist = fabs(goal_distance)/2;
		speed = sqrt(rampR*accel_dist);
	}
	else
	{
		accel_dist = (float)speed*speed/rampR;
		decel_dist = (float)speed*speed/rampR;
	}
	while( fabs(total_right-total_left)<fabs(goal_distance))
	{
		t0=t;
		//Accel/Decel Speed Set
		if (((total_right-total_left)-goal_distance)<0)
			sens = 1;
		else
			sens = -1;
		if (fabs((total_right-total_left)) < accel_dist)
			speed_ref = sens*constrainn(sqrt(rampR*fabs(total_right-total_left)),50,1000);
		else if (fabs((total_right-total_left)-goal_distance) < decel_dist)
			speed_ref = sens*constrainn(sqrt(rampR*fabs((total_right-total_left)-goal_distance)),10,1000);
		else
			speed_ref = sens*speed;
		//Right wheel regulation
		right_error = speed_ref - right_speed;
		i_right_error += right_error;
		PWM_RB = kp * right_error + ki * i_right_error;
		//Left wheel regulation
		left_error = - speed_ref - left_speed;
		i_left_error += left_error;
		PWM_LB = kp * left_error + ki * i_left_error;
		//Position Correction;
		left_correction = coef_correct_angle * (total_right + total_left);
		right_correction = - left_correction;
		PWM_R = PWM_RB + right_correction;
		PWM_L = PWM_LB - left_correction;
		//Execution
		run_motors();
		do delta=t-t0;
		while (delta<T);
	}
	stop_motors();
}

void orientate (float orientation, float speed)
{
	target_angle = orientation - current_phi_deg;
	if (target_angle>180) target_angle -= 360;
	if (target_angle<-180) target_angle += 360;
	rotate(target_angle,speed);
}

void Robot_Locate(float goal_x, float goal_y, float speed)
{
	sens = (asinf((goal_y-current_y)/sqrtf((current_x-goal_x)*(current_x-goal_x)+(current_y-goal_y)*(current_y-goal_y)))>0)? 1 : -1;
	target_angle = sens * rad_to_deg(acosf((goal_x-current_x)/sqrtf((current_x-goal_x)*(current_x-goal_x)+(current_y-goal_y)*(current_y-goal_y))));
	orientate(target_angle,speed);
	goal_distance = sqrtf((current_x-goal_x)*(current_x-goal_x)+(current_y-goal_y)*(current_y-goal_y));
	move_distance (goal_distance,speed);
}

// 3ADDEL BEHA
void move(float speed, int delay)
{
	i_right_error = 0;
	i_left_error = 0;
	t = 0;
	while (t<delay)
	{
		t0=t;
		speed_ref = ramp*t/1000;
		if (speed_ref>speed)
			speed_ref = speed;
		//Right wheel regulation
		right_error = speed_ref - right_speed;
		i_right_error += right_error;
		PWM_R = kp * right_error + ki * i_right_error;
		//Left wheel regulation
		left_error = speed_ref - left_speed;
		i_left_error += left_error;
		PWM_L = kp * left_error + ki * i_left_error;
		run_motors();
		do delta=t-t0;
		while (delta<T);
	}
	stop_motors();
}



								///////////////////////////////////////////////////////////////////////////
								///////////////////////// FOR THE NEXT GENERATION /////////////////////////
								///////////////////////// 			RABI Y3INKOM 			/////////////////////////
								///////////////////////////////////////////////////////////////////////////

void trajectory_genrator(void)
{
	target_y = 0; 
	target_y_prime = 0;
	if (t<2000)
	{
		target_x = /*dec+*/(float)0.5*ramp*t*t/1000000;
		target_x_prime = ramp*t/1000;
	}
	else if (t<2500)
	{
		target_x = /*dec + */1000 + (t-2000);
		target_x_prime = 1000;
		
	}
	else if (t<4500)
	{
		target_x = /*dec +*/ 1500 + (t-2500) - (float)0.5*ramp*(t-2500)*(t-2500)/1000000;
		target_x_prime = 1000 - ramp * (t-2500)/1000;
	}
	else
	{
		target_x = /*dec + */2500;
		target_x_prime = 0;
	}
}

void trajectory (void)
{
	initt();
	float kVx = 0.5, kVy = 0.5, kPx = 0, kPy = 0;
	double pass_matrix[2][2], correction_matrix[2];
	while(1)
	{
		t0=t;
		trajectory_genrator();
		//Calculate the correction matrix (x,y)
		correction_matrix[0] = kVx*target_x_prime + kPx*(target_x-ref_x);
		correction_matrix[1] = kVy*target_y_prime + kPy*(target_y-ref_y);
		//Calculate the passing matrix from coordinate system to speed system
		pass_matrix[0][0] = /*dec**/cos(current_phi_rad) + spacing_wheel*sin(current_phi_rad)/2;
		pass_matrix[0][1] = /*dec**/sin(current_phi_rad) - spacing_wheel*cos(current_phi_rad)/2;
		pass_matrix[1][0] = /*dec**/cos(current_phi_rad) - spacing_wheel*sin(current_phi_rad)/2;
		pass_matrix[1][1] = /*dec**/sin(current_phi_rad) + spacing_wheel*cos(current_phi_rad)/2;
		//Calculate the right and left speed
		left_target_speed = (pass_matrix[0][0]*correction_matrix[0] + pass_matrix[0][1]*correction_matrix[1])/*/dec*/;
		right_target_speed = (pass_matrix[1][0]*correction_matrix[0] + pass_matrix[1][1]*correction_matrix[1])/*/dec*/;
		//Speed regulation
		//Right wheel regulation
		right_error = right_target_speed - right_speed;
		i_right_error += right_error;
		PWM_R = kp * right_error + ki * i_right_error;
		//Left wheel regulation
		left_error = left_target_speed - left_speed;
		i_left_error += left_error;
		PWM_L = kp * left_error + ki * i_left_error;
		//Execute
		run_motors();
		do delta=t-t0;
		while (delta<T);
	}
}


void Time_interrupt(){
    update_position();
    // PWM_sign_change_counter();
   
}

void setup(){
Serial.begin(9600);
// Timer1.initialize(1000);
// Timer1.attachIn40);
set_dimentions(45,45,355,185);
set_motors(255,7,8,9,10);
set_right_encoder(1024,4,1,18,19);
set_left_encoder(1024,4,1,20,21);
set_PWM_min(80,80,80,80);
}


void loop(){
 	 update_position();
    // PWM_sign_change_counter();
    Serial.print("x= ");
    Serial.print(ref_x);
    Serial.print(" ");
    Serial.print("y= ");
    Serial.print(ref_y);
    Serial.print(" ");
    Serial.print("teta= ");
    Serial.println(current_phi_deg);
    delay(1);


}