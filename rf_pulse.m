function [RF]=rf_pulse(B1,dB,tau_p,gyro)

dw=gyro*dB;
w1=gyro*B1;
w_eff=sqrt(w1*w1+dw*dw);
alpha_=w_eff*tau_p;


RF=zeros(3,3);

%% Mx
RF(1,1)=(1-(dw/w_eff)^2)+(dw/w_eff)^2*cos(alpha_);
RF(1,2)=dw/w_eff*sin(alpha_);
RF(1,3)=dw*w1/w_eff^2*(1-cos(alpha_));

%% My
RF(2,1)=-dw/w_eff*sin(alpha_);
RF(2,2)=cos(alpha_);
RF(2,3)=w1/w_eff*sin(alpha_);

%% Mz
RF(3,1)=w1*dw/w_eff^2*(1-cos(alpha_));
RF(3,2)=-w1/w_eff*sin(alpha_);
RF(3,3)=1+(w1/w_eff)^2*(cos(alpha_)-1);

%% Remove rotation to be added.

Rcorr=[cos(dw*tau_p) sin(dw*tau_p) 0
      -sin(dw*tau_p) cos(dw*tau_p) 0
      0                   0        1];

Rcorr_m1=[cos(dw*tau_p) -sin(dw*tau_p) 0
      sin(dw*tau_p) cos(dw*tau_p) 0
      0                   0        1];

RF=Rcorr_m1*RF*Rcorr;