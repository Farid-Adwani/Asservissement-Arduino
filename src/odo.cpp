#include <stdlib.h>
#include<Arduino.h>
#include<TimerOne.h>
#include <ros.h>
#include <ros/time.h>
#include <tf/tf.h>
#include <tf/transform_broadcaster.h>

ros::NodeHandle  nh;

geometry_msgs::TransformStamped t;
tf::TransformBroadcaster broadcaster;

char base_link[] = "/base_link";
char odom[] = "/odom";

int encoder_rightA=20;      // encoder R
int encoder_rightB =21;
int encoder_leftA =18 ;  // encoder L
int encoder_leftB =19;

volatile long right_delta_ticks = 0;    // encoder 1
volatile long left_delta_ticks = 0;    // encoder 2

unsigned int right_ticks;
int right_resolution;
int right_precision;
int right_sens;
volatile long long current_right_count=0;
volatile long long last_right_count=0;
volatile long d_right;
volatile long total_right_count;
volatile float total_right=0,d_right_counter=0;
volatile double right_speed,right_encoder_speed;


int left_resolution;
int left_precision;
int left_sens;


// extern volatile  int right_delta_ticks;
// extern volatile  int left_delta_ticks;
volatile long long current_left_count;
volatile long long last_left_count;
volatile long d_left;
volatile long total_left_count;
volatile float total_left=0,d_left_counter=0,total_centre=0;
volatile double left_speed,left_encoder_speed;

volatile float dR,dL,dC;
volatile float left_radius,right_radius,spacing_encoder,spacing_wheel;
volatile float current_x=0;
volatile float current_y=0;
volatile double current_phi_rad=0;
volatile double phi_speed,d_phi_counter;
volatile float current_phi_deg;



float ref_x,ref_y;

void reset_position()
{
	current_x = 0;
	current_y = 0;
	current_phi_deg = 0;
	current_phi_rad = 0;
}

void initt (void)
{
	total_right=0;
	total_left=0;
	total_centre=0;
}

float rad_to_deg(double x)
{
	return (x*360/(2*PI));
}

void set_right_encoder(int resolution, int precision, int sens,int pinA,int pinB)
{
	
	encoder_rightA=pinA;
	encoder_rightB=pinB;
	right_resolution = resolution;
	right_precision = precision;
	right_sens=sens;
	right_delta_ticks = 0;
	
	
}

void set_left_encoder(int resolution, int precision, int sens,int pinA,int pinB)
{

	encoder_leftA=pinA;
	encoder_leftB=pinB;
	left_resolution = resolution;
	left_precision = precision;
	left_sens=sens;
	left_delta_ticks = 0;

}

void read_right_encoder()
{
	last_right_count = current_right_count;
	current_right_count = right_sens*right_delta_ticks;
	d_right = current_right_count - last_right_count;
  
}

void read_left_encoder()
{
	last_left_count = current_left_count;
	current_left_count = left_sens*left_delta_ticks;
	d_left = current_left_count - last_left_count;
}

void set_dimentions(float right_wheel_radius, float left_wheel_radius, float encoder_spacing, float wheels_spacing)
{
	right_radius = right_wheel_radius;
	left_radius = left_wheel_radius;
	spacing_encoder = encoder_spacing;
	spacing_wheel = wheels_spacing;
}

float ticks_to_distance(int x, float r, int resolution, int precision) 
{
	return (x*2*PI*r/(resolution*precision));
}

void update_position()
{
	read_right_encoder();
	read_left_encoder();
	dR = ticks_to_distance(d_right,right_radius,right_resolution,right_precision);
	total_right += dR;
	d_right_counter += dR;
	dL = ticks_to_distance(d_left,left_radius,left_resolution,left_precision);
	total_left += dL;
	d_left_counter += dL;
	dC = (dR+dL)/2;
	total_centre+=dC;
	current_x += dC*cos(current_phi_rad);
	current_y += dC*sin(current_phi_rad);
	current_phi_rad += ((dR-dL)/spacing_encoder);
	d_phi_counter += ((dR-dL)/spacing_encoder);
	while (current_phi_rad>PI)
	{
		current_phi_rad -= 2*PI;
	}
	while (current_phi_rad<-PI)
	{
		current_phi_rad += 2*PI;
	}
	current_phi_deg = rad_to_deg(current_phi_rad);
}


void Time_interrupt(){
    update_position(); 
}

void doEncoderA() {
  // look for a low-to-high on channel A
  if (digitalRead(encoder_rightA) == HIGH) {
    // check channel B to see which way encoder is turning
    if (digitalRead(encoder_rightB) == LOW) {
      right_delta_ticks = right_delta_ticks + 1;         // CW
    }
    else {
      right_delta_ticks = right_delta_ticks - 1;         // CCW
    }
  }
  else   // must be a high-to-low edge on channel A
  {
    // check channel B to see which way encoder is turning
    if (digitalRead(encoder_rightB) == HIGH) {
      right_delta_ticks = right_delta_ticks + 1;          // CW
    }
    else {
      right_delta_ticks = right_delta_ticks - 1;          // CCW
    }
  }

}

void doEncoderB() {

  // look for a low-to-high on channel B
  if (digitalRead(encoder_rightB) == HIGH) {
    // check channel A to see which way encoder is turning
    if (digitalRead(encoder_rightA) == HIGH) {
      right_delta_ticks = right_delta_ticks + 1;         // CW
    }
    else {
      right_delta_ticks = right_delta_ticks - 1;         // CCW
    }
  }
  // Look for a high-to-low on channel B
  else {
    // check channel B to see which way encoder is turning
    if (digitalRead(encoder_rightA) == LOW) {
      right_delta_ticks = right_delta_ticks + 1;          // CW
    }
    else {
      right_delta_ticks = right_delta_ticks - 1;          // CCW
    }
  }


}

void doEncoderC() {

  // look for a low-to-high on channel A
  if (digitalRead(encoder_leftA) == HIGH) {
    // check channel B to see which way encoder is turning
    if (digitalRead(encoder_leftB) == LOW) {
      left_delta_ticks = left_delta_ticks - 1;         // CW
    }
    else {
      left_delta_ticks = left_delta_ticks + 1;         // CCW
    }
  }
  else   // must be a high-to-low edge on channel A
  {
    // check channel B to see which way encoder is turning
    if (digitalRead(encoder_leftB) == HIGH) {
      left_delta_ticks = left_delta_ticks - 1;          // CW
    }
    else {
      left_delta_ticks = left_delta_ticks + 1;          // CCW
    }
  }

}

void doEncoderD() {

  // look for a low-to-high on channel B
  if (digitalRead(encoder_leftB) == HIGH) {
    // check channel A to see which way encoder is turning
    if (digitalRead(encoder_leftA) == HIGH) {
      left_delta_ticks = left_delta_ticks - 1;         // CW
    }
    else {
      left_delta_ticks = left_delta_ticks + 1;         // CCW
    }
  }
  // Look for a high-to-low on channel B
  else {
    // check channel B to see which way encoder is turning
    if (digitalRead(encoder_leftA) == LOW) {
      left_delta_ticks = left_delta_ticks - 1;          // CW
    }
    else {
      left_delta_ticks = left_delta_ticks + 1;          // CCW
    }
  }


}

void setup(){
  nh.getHardware()->setBaud(115200);
nh.initNode();
broadcaster.init(nh);
// Serial.println(9600);

initt();
reset_position();
Timer1.initialize(10000);
Timer1.attachInterrupt(Time_interrupt);
set_dimentions(40,40,245,190);
set_right_encoder(400,4,1,20,21);
set_left_encoder(400,4,1,18,19);

  pinMode(encoder_rightA,INPUT_PULLUP);
  pinMode(encoder_rightB,INPUT_PULLUP);
  pinMode(encoder_leftA,INPUT_PULLUP);
  pinMode(encoder_leftB,INPUT_PULLUP);
  attachInterrupt(digitalPinToInterrupt(encoder_leftA), doEncoderC, CHANGE);
  attachInterrupt(digitalPinToInterrupt(encoder_leftB), doEncoderD, CHANGE);
  attachInterrupt(digitalPinToInterrupt(encoder_rightA), doEncoderA, CHANGE);
  attachInterrupt(digitalPinToInterrupt(encoder_rightB), doEncoderB, CHANGE);
}


void loop(){
    // update_position();
  t.header.frame_id = odom;
  t.child_frame_id = base_link;
  
  t.transform.translation.x = current_x/10;
  t.transform.translation.y = current_y/10;
  
  t.transform.rotation = tf::createQuaternionFromYaw(current_phi_rad);

  t.header.stamp = nh.now();
  broadcaster.sendTransform(t);
  nh.spinOnce();
  delay(10);
    // Serial.print(current_x);
    // Serial.print(" , ");
    // Serial.print(current_y);
	  // Serial.print(" , ");
    // Serial.println(current_phi_deg);
    // delay(1);

}