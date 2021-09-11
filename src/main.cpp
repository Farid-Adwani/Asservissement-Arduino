// #include <Arduino.h>
// #include <ros.h>
// #include <std_msgs/String.h>
// #include <Motors.h>



// #define led_pin 3

// ros::NodeHandle  nh;
// std_msgs::String str_msg;
// ros::Publisher chatter("odometry", &str_msg); 

// void messageCb( const std_msgs::String& msg){
//   Serial.println(msg.data);
//    const char* msgg=msg.data;
//   if(!strcmp(msgg,"on")){
//    digitalWrite(led_pin, HIGH);  
//    str_msg.data="ON ok :)";
//   chatter.publish( &str_msg );  
//     }else
//     if(!strcmp(msgg,"off")){
//    digitalWrite(led_pin,LOW); 
//    str_msg.data="OFF ok :(";
//   chatter.publish( &str_msg );    
//     }else
//     if(strcmp(msgg,"blink")==0){
//    for(int i=0;i<10;i++){
//     digitalWrite(led_pin, HIGH);
//     delay(1000);
//     digitalWrite(led_pin, LOW);
//     delay(1000);
//     }}
//     else{
//      for(int i=0;i<5;i++){
//     digitalWrite(led_pin, HIGH);
//     delay(1000);
//     digitalWrite(led_pin, LOW);
//     delay(1000);
//     } 
//       }    
    


// }


// ros::Subscriber<std_msgs::String> sub("led", &messageCb );


// void setup()
// { 
//   Serial.begin(96000);
//   pinMode(LED_BUILTIN, OUTPUT);
//   pinMode(led_pin, OUTPUT);
//   nh.initNode();
//   nh.subscribe(sub);
//   nh.advertise(chatter);
  
// }

// void loop()
// {  
//   nh.spinOnce();
//   delay(1);
// }
