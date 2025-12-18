clear all
close all
gyro=2.675*1e8;

addpath('C:\Users\Gaddiel.Ouaknin\Documents\BSSFP')
addpath('C:\Users\Gaddiel.Ouaknin\Documents\EddyCurrent\')

seq=sequence;


f = 'meas_MID00204_FID12351_AdjTra.dat';
p = 'C:\Users\Gaddiel.Ouaknin\Documents\TXADJUST\SiemensData';
twix = mapVBVD( [ p filesep f ] ,'removeOS',true);
data=twix{1}.image(:,:,:,:,:);
data=squeeze(data);
number_of_lines = twix{1}.image.NLin;
number_of_columns = min(twix{1}.image.NCol,size(data,1));
sampling_time=10*1e-6;
Nsamples=128;
t_sampled=(0:1:(Nsamples-1))*sampling_time;
figure(1)
plot(t_sampled, abs(data(:,1,1)),'.')
hold on
plot(t_sampled,abs(data(:,1,2)),'*')
xlabel('time')
ylabel('Signal (ABS)')
title('Siemens Data')
legend('SE','STE')

input_data_folder='\Tx 18 Dec 9 am\'; channels=9;
mrs=MRSparser();
S=mrs.get_signal(input_data_folder,'TxCalibration',channels);
figure(2)
plot(t_sampled,abs(S(1:Nsamples,1)),'.');
hold on
plot(t_sampled,abs(S(Nsamples+1:2*Nsamples,1)),'*');
xlabel('time')
ylabel('Signal (ABS)')
title('MRS Data')
legend('SE','STE')

figure(3)
plot(t_sampled, abs(data(:,1,1))/max(abs(data(:,1,1))),'.')
hold on
plot(t_sampled,abs(S(1:Nsamples,1))/max(abs(S(1:Nsamples,1))),'+');
xlabel('time')
ylabel('Signal (ABS)')
title('SE')
legend('Siemens','MRS')


figure(4)
plot(t_sampled, abs(data(:,1,2))/max(abs(data(:,1,2))),'.')
hold on
plot(t_sampled,abs(S(Nsamples+1:2*Nsamples,1))/max(abs(S(Nsamples+1:2*Nsamples,1))),'+');
xlabel('time')
ylabel('Signal (ABS)')
title('STE')
legend('Siemens','MRS')


peak_se=max(abs(S(1:Nsamples,1)));
peak_ste=max(abs(S(Nsamples+1:2*Nsamples,1)));
flip_angle=acos(peak_ste/peak_se);

Amp=pi/2/flip_angle;
flip_angle_deg=180*flip_angle/pi;

t1=6.92*1e-3;
t2=23.04*1e-3;
rf_1=pi/2;
rf_2=2*rf_1;
rf_3=rf_1;
Meqm=[0 0 1]';
seq.Gz=0.51*1e-3;
phi_1=gyro*t1*seq.Gz;
phi_2=gyro*t2*seq.Gz;


z=(-0.1:1e-5:0.1); z=z';
Nz=length(z);
ME1=zeros(Nz,3);
%% Compute First Echo %%
tau_p=500*1e-6;
seq=seq.generate_rx(rf_1);
B1=rf_1/gyro/tau_p;
for n=1:Nz
    [ME1(n,1),ME1(n,2),ME1(n,3)]= ...
        rf_pulse(Meqm(1),Meqm(2),Meqm(3),B1,seq.Gz*z(n),tau_p,gyro);
    seq=seq.generate_rz(phi_1*z(n));
    ME1(n,:)=seq.rotate_z(ME1(n,:)');
end
seq=seq.generate_rx(rf_2);
for n=1:Nz
    [ME1(n,1),ME1(n,2),ME1(n,3)]= ...
        rf_pulse(ME1(n,1),ME1(n,2),ME1(n,3),B1,seq.Gz*z(n),2*tau_p,gyro);
    seq=seq.generate_rz(phi_1*z(n));
    ME1(n,:)=seq.rotate_z(ME1(n,:)');
end

line_1_simulated=zeros(Nsamples,3);
for it=1:Nsamples
    Mit=ME1;
   for n=1:Nz
        seq=seq.generate_rz((it-1-Nsamples/2)*sampling_time*z(n)*gyro);
        Mit(n,:)=seq.rotate_z(Mit(n,:)');
   end
   % if it-1==Nsamples/2
   % disp('hi')
   % end
    line_1_simulated(it,:)=sum(Mit)/Nz;
end

figure(5)
plot(sqrt(line_1_simulated(:,1).^2+line_1_simulated(:,2).^2));

ME2=ME1;
seq=seq.generate_rx(rf_3);
for n=1:Nz
    seq=seq.generate_rz((phi_2-phi_1)*z(n));
    ME2(n,:)=seq.rotate_z(ME2(n,:)');
[ME2(n,1),ME2(n,2),ME2(n,3)]= ...
        rf_pulse(ME2(1),ME2(2),ME2(3),B1,seq.Gz*z(n),tau_p,gyro);
    seq=seq.generate_rz(phi_1*z(n));
    ME2(n,:)=seq.rotate_z(ME2(n,:)');
end


line_2_simulated=zeros(Nsamples,3);
for it=1:Nsamples
    Mit=ME2;
    for n=1:Nz
        seq=seq.generate_rz((it-1-Nsamples/2)*sampling_time*z(n)*gyro);
        Mit(n,:)=seq.rotate_z(Mit(n,:)');
    end
    line_2_simulated(it,:)=sum(Mit)/Nz;
end

figure(6)
plot(sqrt(line_2_simulated(:,1).^2+line_2_simulated(:,2).^2));


expected_ratio=0.5*sin(rf_2)*sin(rf_3)/sin(rf_2/2)/sin(rf_2/2);

Sx1=sum(ME1(:,1))/Nz;
Sy1=sum(ME1(:,2))/Nz;
Sx2=sum(ME2(:,1))/Nz;
Sy2=sum(ME2(:,2))/Nz;

simulated_ratio=sqrt(Sx2*Sx2+Sy2*Sy2)/sqrt(Sx1*Sx1+Sy1*Sy1);


