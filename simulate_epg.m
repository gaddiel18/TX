clear all
close all

gyro=2.675*1e8;
t1=7*1e-3;
t2=23*1e-3;
k1=gyro*t1;
k2=gyro*t2;
k1_=round(1e12*k1);
k2_=round(1e12*k2);
k3_=k1_;
S=[ 1 1i 0 
    1 -1i 0 
    0 0 1];
Sm1=[0.5 0.5 0
   1/(2*1i) -1/(2*1i) 0
    0 0 1];

alpha_1=pi/2;
Rx1=[1 0 0
    0 cos(alpha_1) sin(alpha_1)
    0 -sin(alpha_1) cos(alpha_1)];

alpha_3=pi/3;
Rx1=[1 0 0
    0 cos(alpha_1) sin(alpha_1)
    0 -sin(alpha_1) cos(alpha_1)];

alpha_2=pi;
Rx2=[1 0 0
    0 cos(alpha_2) sin(alpha_2)
    0 -sin(alpha_2) cos(alpha_2)];


Ns=30;

A=zeros(Ns,5);

A(1,1)=1;
A(1,5)=1;
Na=1;

%% First Pulse Action
for k=1:Na
    A(k,3:5)=S*Rx1*Sm1*A(k,3:5)';
end

%% First Interpulse Shift
dk=k1_;
A=shift(A,dk,Na);
Na=Na*3;

%% Combine & Merge States
A=combine_and_merge(A,Na);

%% Tight to beginning
[A,Na]=tight(A,Na,Ns);

%% Second Pulse Action
for k=1:Na
    A(k,3:5)=S*Rx2*Sm1*A(k,3:5)';
end

%% Second Interpulse Shift
dk=k2_;
A=shift(A,dk,Na);
Na=Na*3;

%% Combine & Merge States
A=combine_and_merge(A,Na);

%% Tight to beginning
[A,Na]=tight(A,Na,Ns);

%% Third Pulse Action
for k=1:Na
    A(k,3:5)=S*Rx1*Sm1*A(k,3:5)';
end

%% Go to Echo
dk=k3_;
A=shift(A,dk,Na);
Na=Na*3;

%% Combine & Merge States
A=combine_and_merge(A,Na);

%% Tight to beginning
[A,Na]=tight(A,Na,Ns);

Lz=0.1;
G=0.00051;
%% Compute Signal
Nsamples=128;
tn=10*1e-6*(-Nsamples/2:1:Nsamples/2-1)';
St=zeros(Nsamples,1);
for n=1:Nsamples
    wk=zeros(Na,1);
    for k=1:Na
        wk(k)=k_weight(A(k,2)/1e12,Lz,gyro,tn(n),G);
    end
    St(n)=sum(wk);
end

function A=shift(A,dk,Na)

for k=1:Na

    %% Shift of positive state
    A(Na+k,1)=1;
    A(Na+k,2)=A(k,2)-dk;
    A(Na+k,3)=A(k,3);
    A(k,3)=0;

    %% Shift of negative state
    A(2*Na+k,1)=1;
    A(2*Na+k,2)=A(k,2)+dk;
    A(2*Na+k,4)=A(k,4);
    A(k,4)=0;

end
end



function A=combine_and_merge(A,Na)

for k=1:Na
    for k2=k+1:Na
        if(A(k,2)==A(k2,2))
            A(k,3:5)=A(k,3:5)+A(k2,3:5);
            A(k2,:)=0;
        end
    end
end

end

function [A,Na]=tight(A,Na,Ns)

    Na_=sum(A(:,1));
    A_=zeros(Ns,5);
    k_=1;
    for k=1:Na
        if(A(k,1)==1)
        A_(k_,:)=A(k,:);
        k_=k_+1;
        end
    end

    Na=Na_;
    A=A_;
end

function [wk]=k_weight(k,Lz,gyro,tn,G)
        k=k-gyro*tn*G;         
        wk=sin(k*Lz/2)/(k*Lz/2);
end