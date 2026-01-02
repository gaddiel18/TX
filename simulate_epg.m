clear all
close all
set(0,'DefaultFigureWindowStyle','docked')

gyro=2.675*1e8;
t1=7*1e-3;
t2=23*1e-3;
k1=gyro*t1;
k2=gyro*t2;
k1_=round(1e12*k1);
k2_=round(1e12*k2);
k3_=k1_;
Lz=0.1;
G=0.00051;
B1=1e-5;
tau_p=1000*1e-6;
dz=0.001;
z=(-Lz/2:dz:Lz/2);
Nk=30;
%% Compute Signal
Nsamples=512;
tn=10*1e-6*(-Nsamples/2:1:Nsamples/2-1)';
Nz=length(z);

S=[ 1 1i 0
    1 -1i 0
    0 0 1];
Sm1=[0.5 0.5 0
    1/(2*1i) -1/(2*1i) 0
    0 0 1];

alpha_1=pi/2; alpha_2=pi; alpha_3=pi/3;


Rx1=[1 0 0
    0 cos(alpha_1) sin(alpha_1)
    0 -sin(alpha_1) cos(alpha_1)];



Rx2=[1 0 0
    0 cos(alpha_2) sin(alpha_2)
    0 -sin(alpha_2) cos(alpha_2)];

Rx3=[1 0 0
    0 cos(alpha_3) sin(alpha_3)
    0 -sin(alpha_3) cos(alpha_3)];

zi=-0.05;
[A,Na]=get_stimulated_echo_one_voxel(Nk,B1,zi,tau_p,gyro,G,S,Sm1,k1_,k2_,k3_);
St=get_echo_from_epg(Nsamples,G,gyro,dz,zi,A,Na,tn);
figure(1)
plot(real(St))
hold on
plot(imag(St),'r')
hold on
plot(abs(St),'k')
title('Stimulated Echo')
legend('real','imag','abs')


zi=0.05;
[A2,Na2]=get_stimulated_echo_one_voxel(Nk,B1,zi,tau_p,gyro,G,S,Sm1,k1_,k2_,k3_);
St2=get_echo_from_epg(Nsamples,G,gyro,dz,zi,A2,Na2,tn);
figure(2)
plot(real(St2))
hold on
plot(imag(St2),'r')
hold on
plot(abs(St2),'k')
title('Stimulated Echo')
legend('real','imag','abs')

figure(3)
plot(real(St+St2))
hold on
plot(imag(St+St2),'r')
hold on
plot(abs(St+St2),'k')
title('Stimulated Echo')
legend('real','imag','abs')

Signal=zeros(Nsamples,Nz);
for i=1:Nz
    [Ai,Nai]=get_stimulated_echo_one_voxel(Nk,B1,z(i),tau_p,gyro,G,S,Sm1,k1_,k2_,k3_);
    Signal(:,i)=get_echo_from_epg(Nsamples,G,gyro,dz,z(i),Ai,Nai,tn);
end

figure(4)
plot(real(Signal))
hold on
plot(imag(Signal),'r')
hold on
plot(abs(Signal),'k')
title('Stimulated Echo')
legend('real','imag','abs')


S_=sum(Signal,2)/Nz;
figure(5)
plot(real(S_))
hold on
plot(imag(S_),'r')
hold on
plot(abs(S_),'k')
title('Stimulated Echo')
legend('real','imag','abs')



function [A,Na]=shift(A,dk,Na,Nk)
A=shift_(A,dk,Na);
Na=Na*3;
%% Combine & Merge States
A=combine_and_merge(A,Na);
%% Tight to beginning
[A,Na]=tight(A,Na,Nk);
end

function A=shift_(A,dk,Na)

for k=1:Na

    %% Shift of positive state
    A(Na+k).active=1;
    A(Na+k).ki=A(k).ki-dk;
    A(Na+k).F(1)=A(k).F(1);
    A(k).F(1)=0;

    %% Shift of negative state
    A(2*Na+k).active=1;
    A(2*Na+k).ki=A(k).ki+dk;
    A(2*Na+k).F(2)=A(k).F(2);
    A(k).F(2)=0;

end
end



function A=combine_and_merge(A,Na)

for k=1:Na
    for k2=k+1:Na
        if(A(k).ki==A(k2).ki)
            A(k).F=A(k).F+A(k2).F;
            A(k2).F=0;
            A(k2).active=0;
        end
    end
end

end

function [A_,Na_]=tight(A,Na,Nk)

Na_=sum([A.active]);
A_=createArray(Nk,1,"epg");
k_=1;
for k=1:Na
    if(A(k).active==1)
        A_(k_).active=A(k).active;
        A_(k_).F=A(k).F;
        A_(k_).ki=A(k).ki;
        k_=k_+1;
    end
end

end

function [wk]=k_weight(k,dz,zc,gyro,tn,G)
k=k-gyro*tn*G;
if abs(k)>0
    wk=sin(k*dz/2)/(k*dz/2)*exp(1i*k*zc);
else
    wk=1;
end
end


function [St]=get_echo_from_epg(Nsamples,G,gyro,dz,zc,A,Na,tn)
St=zeros(Nsamples,1);
for n=1:Nsamples
    wk=zeros(Na,1);
    for k=1:Na
        wk(k)=k_weight(A(k).ki/1e12,dz,zc,gyro,tn(n),G)*A(k).F(1);
    end
    St(n)=sum(wk);
end
end



function [A,Na]=get_stimulated_echo_one_voxel(Nk,B1,zi,tau_p,gyro,G,S,Sm1,k1_,k2_,k3_)


A=createArray(Nk,1,"epg");
A(1).F=[0, 0,1]';
A(1).active=1;
Na=1;
dB=G*zi;
%% First Pulse Action
for k=1:Na
    RF=rf_pulse(B1,dB,tau_p/2,gyro);
    A(k).F=S*RF*Sm1*A(k).F;
end

%% First Interpulse Shift
dk=k1_;
[A,Na]=shift(A,dk,Na,Nk);

%% Second Pulse Action
for k=1:Na
    RF=rf_pulse(B1,dB,tau_p,gyro);
    A(k).F=S*RF*Sm1*A(k).F;
end



%% Second Interpulse Shift
dk=k2_;
[A,Na]=shift(A,dk,Na,Nk);

%% Third Pulse Action
for k=1:Na
    RF=rf_pulse(B1,dB,tau_p/2,gyro);
    A(k).F=S*RF*Sm1*A(k).F;
end

%% Go to Echo
dk=k3_;
[A,Na]=shift(A,dk,Na,Nk);
end