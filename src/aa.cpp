// #include<Arduino.h>

// #define encoder1PinA 20      // encoder R
// #define encoder1PinB 21

// #define encoder0PinA 18   // encoder L
// #define encoder0PinB 19

// unsigned long currentMillis;
// unsigned long previousArmMillis;
// unsigned long previousMillis;
// int ki=0;
// int kp=0;
// int kd=0;
// int error;
// int prverror;
// int pwm_B=60;
// int pwm_A=60;
// volatile long encoder0Pos = 0;    // encoder 1
// volatile long encoder1Pos = 0;    // encoder 2
// int correction;
// float x = 0.0;
// float y = 0.0;
// float teta = 0.0;
// float resolution = 4096.0000;
// float dia = 80.600;
// float dl = 0.0000;
// float dr = 0.0000;
// float deltaX = 0.0000;
// float entreAxe = 374.000;
// float ndc = 0.0000;
// float arctan = 0.00;
// float dc = 0.0000;
// float actualSpeed = 0.0000;
// float speedX = 0.0000;
// float speedY = 0.0000;
// float speedTheta = 0.0000;
// float vl = 0.0000;
// float vr = 0.0000;
// float angleTheta = 0.0000;

// double radToDeg(double rad) {
//   return (rad * 180.00) / PI;
// }


// void updatePosition() {
//   dl = PI * dia * (encoder0Pos / resolution);
//   dr = PI * dia * (encoder1Pos / resolution);
//   dc = (dr + dl) / 2.0000;
//   teta = (dr - dl) / entreAxe;
//   deltaX = dc - ndc;
//   y += deltaX * cos(teta);
//   x += deltaX * sin(teta) * -1.0000;
//   teta = atan2(sin(teta), cos(teta));
//   if (teta < 0) {
//     teta += 2 * PI;
//   }
//   ndc = dc;
// }


// int mr1 = A0;
// int mr2 = A1;
// int ml1 = A2;
// int ml2 = A3;

// void doEncoderA();
// void doEncoderB();
// void doEncoderC();
// void doEncoderD();

// void setup() {
//   pinMode(mr1, OUTPUT);
//   pinMode(mr2, OUTPUT);
//   pinMode(ml1, OUTPUT);
//   pinMode(ml2, OUTPUT);
//   pinMode(encoder1PinA,INPUT_PULLUP);
//   pinMode(encoder1PinB,INPUT_PULLUP);
//   pinMode(encoder0PinA,INPUT_PULLUP);
//   pinMode(encoder0PinB,INPUT_PULLUP);
//   attachInterrupt(digitalPinToInterrupt(encoder0PinA), doEncoderA, CHANGE);
//   attachInterrupt(digitalPinToInterrupt(encoder0PinB), doEncoderB, CHANGE);

//   attachInterrupt(digitalPinToInterrupt(encoder1PinA), doEncoderC, CHANGE);
//   attachInterrupt(digitalPinToInterrupt(encoder1PinB), doEncoderD, CHANGE);
//   Serial.begin(9600);
// }



// void loop() {
//   currentMillis = millis();   // bookmark the time
//   updatePosition();
//   if (currentMillis - previousMillis >= 10) {  // start timed loop for everything else 
//     previousMillis = currentMillis;
//     Serial.print(encoder0Pos);
//     Serial.print(" , ");
//     Serial.println(encoder1Pos);
//    /* Serial.print(" || ");
//     Serial.print(x);
//     Serial.print(" , ");
//     Serial.print(y);
//     Serial.print(" , ");
//     Serial.println(radToDeg(teta));
//   }
// error=encoder0Pos-encoder1Pos;
// correction=kp*error+ki*(error+prverror)+kd*(error-prverror);
// pwm_A+=correction;
// pwm_B-=correction;
// analogWrite(mr1,pwm_A); 
// analogWrite(ml1,pwm_B);*/





// /* if(millis()<10000000)
//   {

//     if(encoder0Pos<encoder1Pos) digitalWrite(ml1,150);
//     else if(encoder0Pos>encoder1Pos) digitalWrite(mr1,150);
// }
// else{
//   digitalWrite(mr1,0);
//   digitalWrite(ml1,0);
// }*/
// }}










// void doEncoderA() {

//   // look for a low-to-high on channel A
//   if (digitalRead(encoder0PinA) == HIGH) {
//     // check channel B to see which way encoder is turning
//     if (digitalRead(encoder0PinB) == LOW) {
//       encoder0Pos = encoder0Pos + 1;         // CW
//     }
//     else {
//       encoder0Pos = encoder0Pos - 1;         // CCW
//     }
//   }
//   else   // must be a high-to-low edge on channel A
//   {
//     // check channel B to see which way encoder is turning
//     if (digitalRead(encoder0PinB) == HIGH) {
//       encoder0Pos = encoder0Pos + 1;          // CW
//     }
//     else {
//       encoder0Pos = encoder0Pos - 1;          // CCW
//     }
//   }

// }

// void doEncoderB() {

//   // look for a low-to-high on channel B
//   if (digitalRead(encoder0PinB) == HIGH) {
//     // check channel A to see which way encoder is turning
//     if (digitalRead(encoder0PinA) == HIGH) {
//       encoder0Pos = encoder0Pos + 1;         // CW
//     }
//     else {
//       encoder0Pos = encoder0Pos - 1;         // CCW
//     }
//   }
//   // Look for a high-to-low on channel B
//   else {
//     // check channel B to see which way encoder is turning
//     if (digitalRead(encoder0PinA) == LOW) {
//       encoder0Pos = encoder0Pos + 1;          // CW
//     }
//     else {
//       encoder0Pos = encoder0Pos - 1;          // CCW
//     }
//   }


// }

// // ************** encoder 2 *********************

// void doEncoderC() {

//   // look for a low-to-high on channel A
//   if (digitalRead(encoder1PinA) == HIGH) {
//     // check channel B to see which way encoder is turning
//     if (digitalRead(encoder1PinB) == LOW) {
//       encoder1Pos = encoder1Pos - 1;         // CW
//     }
//     else {
//       encoder1Pos = encoder1Pos + 1;         // CCW
//     }
//   }
//   else   // must be a high-to-low edge on channel A
//   {
//     // check channel B to see which way encoder is turning
//     if (digitalRead(encoder1PinB) == HIGH) {
//       encoder1Pos = encoder1Pos - 1;          // CW
//     }
//     else {
//       encoder1Pos = encoder1Pos + 1;          // CCW
//     }
//   }

// }

// void doEncoderD() {

//   // look for a low-to-high on channel B
//   if (digitalRead(encoder1PinB) == HIGH) {
//     // check channel A to see which way encoder is turning
//     if (digitalRead(encoder1PinA) == HIGH) {
//       encoder1Pos = encoder1Pos - 1;         // CW
//     }
//     else {
//       encoder1Pos = encoder1Pos + 1;         // CCW
//     }
//   }
//   // Look for a high-to-low on channel B
//   else {
//     // check channel B to see which way encoder is turning
//     if (digitalRead(encoder1PinA) == LOW) {
//       encoder1Pos = encoder1Pos - 1;          // CW
//     }
//     else {
//       encoder1Pos = encoder1Pos + 1;          // CCW
//     }
//   }


// }