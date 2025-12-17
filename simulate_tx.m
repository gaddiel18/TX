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

figure(1)
plot(abs([data(:,1,1)' data(:,1,2)']'))


input_data_folder='\Tx 17 Dec\'; channels=9; 
mrs=MRSparser();
S=mrs.get_signal(input_data_folder,'Tx',channels); 
figure(2)
plot(abs(S(:,1)));



t1=6.92*1e-3;
t2=23.04*1e-3;
rf_1=pi/3;
rf_2=2*rf_1;
rf_3=rf_1;
Meqm=[0 0 1]';
seq.Gz=0.51*1e-3;
phi_1=gyro*t1*seq.Gz;
phi_2=gyro*t2*seq.Gz;


z=(-0.1:0.0001:0.1); z=z';
Nz=length(z);
ME1=zeros(Nz,3);
%% Compute First Echo %%
seq=seq.generate_rx(rf_1);
for n=1:Nz
    ME1(n,:)=seq.rotate_x(Meqm);
    seq=seq.generate_rz(phi_1*z(n));
    ME1(n,:)=seq.rotate_z(ME1(n,:)');
end
seq=seq.generate_rx(rf_2);
for n=1:Nz
    ME1(n,:)=seq.rotate_x(ME1(n,:)');
    seq=seq.generate_rz(phi_1*z(n));
    ME1(n,:)=seq.rotate_z(ME1(n,:)');
end

ME2=ME1;
seq=seq.generate_rx(rf_3);
for n=1:Nz
    seq=seq.generate_rz((phi_2-phi_1)*z(n));
    ME2(n,:)=seq.rotate_z(ME2(n,:)');
    ME2(n,:)=seq.rotate_x(ME2(n,:)');
    seq=seq.generate_rz(phi_1*z(n));
    ME2(n,:)=seq.rotate_z(ME2(n,:)');
end



expected_ratio=0.5*sin(rf_2)*sin(rf_3)/sin(rf_2/2)/sin(rf_2/2);

Sx1=sum(ME1(:,1))/Nz;
Sy1=sum(ME1(:,2))/Nz;
Sx2=sum(ME2(:,1))/Nz;
Sy2=sum(ME2(:,2))/Nz;

simulated_ratio=sqrt(Sx2*Sx2+Sy2*Sy2)/sqrt(Sx1*Sx1+Sy1*Sy1);


