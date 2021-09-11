

void set_motors (int MAX_PWM, unsigned int forward_right, unsigned int backward_right, unsigned int forward_left, unsigned int backward_left);
void run_motors (void);
void stop_motors (void);
void run_forward (int right_PWM, int left_PWM);
void set_PWM_min (int PWM_Min_R, int PWM_Min_L, int PWM_Min_R_Rot, int PWM_Min_L_Rot);
void PWM_sign_change_counter(void);
