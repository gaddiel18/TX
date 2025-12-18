function [Mx,My,Mz]=rf_pulse(mx0,my0,mz0,B1,dB,tau_p,gyro)

dw=gyro*dB;
w1=gyro*B1;
w_eff=sqrt(w1*w1+dw*dw);
beta_=dw*mx0-w1*mz0;


Mx=dw/w_eff*my0*sin(w_eff*tau_p);
Mx=Mx+dw/w_eff^2*beta_*cos(w_eff*tau_p)+mx0-dw/w_eff^2*beta_;

My=my0*cos(w_eff*tau_p)-beta_/w_eff*sin(w_eff*tau_p);

Mz=mz0+w1/w_eff^2*beta_-my0*w1/w_eff*sin(w_eff*tau_p)-w1/w_eff^2*beta_*cos(w_eff*tau_p);