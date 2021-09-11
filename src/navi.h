#include "Odometry.h"
#include "Motors.h"
#include "math.h"
#include <stdlib.h>
void move_distance(float distance,float speed); //distance (mm); speed (mm/s)
void rotate(float angle, float speed); //angle (degree) ; speed (mm/s) (Wheel linear speed)
void orientate (float orientation, float speed);
void Robot_Locate(float goal_x, float goal_y, float speed);
void trajectory (void);
void curv(float R,float theta ,float speed);
void Robot_locateCurv(float x, float y, float phi_target ,float speed); //phi_target(degree) phi_target==180 || ==-180 selon le sens souhaité
void Multi_Curv(float R,float theta ,float speed, int i,int n);
void Robot_Locate_Multi_Curv( float** matrix, int n , int speed);//!!!!speed==300 // matrix contains the (x,y) of the desired locations your want the robot to go to in order 

//test fonction
void move(float speed, int delay);
//
void Calcul_Curv(float x, float y,float phi_target);
//
void allocation (int taille);
//



