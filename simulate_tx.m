clear all
close all
gyro=2.675*1e8;

Nsamples=512;
sampling_time=5120*1e-6/Nsamples;
t_sampled=(0:1:(Nsamples-1))*sampling_time;
T1=1;
T2=0.1;

addpath('C:\Users\Gaddiel.Ouaknin\Documents\BSSFP')
addpath('C:\Users\Gaddiel.Ouaknin\Documents\EddyCurrent\')

seq=sequence;
analyze_data=false;
simulate_data=true;
if analyze_data
    f = 'meas_MID00204_FID12351_AdjTra.dat';
    p = 'C:\Users\Gaddiel.Ouaknin\Documents\TXADJUST\SiemensData';
    twix = mapVBVD( [ p filesep f ] ,'removeOS',true);
    data=twix{1}.image(:,:,:,:,:);
    data=squeeze(data);
    number_of_lines = twix{1}.image.NLin;
    number_of_columns = min(twix{1}.image.NCol,size(data,1));
    figure(1)
    plot(t_sampled, abs(data(:,1,1)),'.')
    hold on
    plot(t_sampled,abs(data(:,1,2)),'*')
    xlabel('time')
    ylabel('Signal (ABS)')
    title('Siemens Data')
    legend('SE','STE')

    input_data_folder='\Tx 29 Dec 2pm\'; channels=5;
    mrs=MRSparser();
    S=mrs.get_signal(input_data_folder,'TxCal',channels);
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

    Nsamples=size(S,1);
    Nrf=size(S,2);

    ste_min=zeros(Nrf,1);
    ste_ix=zeros(Nrf,1);


    figure(7)
    for n=1:size(S,2)
        clf
        plot(abs(S(:,n)))
        [ste_min(n),ste_ix(n)]=min(abs(S(round(0.4*Nsamples):round(0.6*Nsamples),n)));
    end

    
    figure(777)
    plot(ste_min)

    [semm,ixx]=min(ste_min);


    figure(77)
    clf
    plot(abs(S(:,ixx)));
    xlabel('Signal with 180 pulse')
    ylabel('Signal')

  kMaxRfScale = 32767;
  noViews_=101;
  minDAC = 0;
  center = round(0.65 * kMaxRfScale);
  x = 4000; 
  rfCalibrationValue_=zeros(noViews_,1);
  calIdx=ixx;
  t = double(calIdx) / (noViews_ - 1);
  val = center + round(-x + 2.0 * x * t);
end

if simulate_data
    tau_p=500*1e-6;
    t1=7*1e-3;
    t2=23*1e-3;
    taup_1=2*tau_p;
    taup_2=2*tau_p;
    taup_3=2*tau_p;
    rf_1=0.5*pi;
    rf_2=rf_1;
    rf_3=rf_1;
    Meqm=[0 0 1]';
    seq.Gz=0.51*1e-3;
    phi_1=gyro*(t1-taup_1/2-taup_2/2)*seq.Gz;
    phi_2=gyro*(t2-taup_1/2-taup_2/2)*seq.Gz;


    R=1/gyro/seq.Gz/sampling_time/2;
    dz=2*R/Nsamples/1000;
    z=(-R:dz:R-dz); z=z';

    rho_z=(1-(z./R).^2);
    Nz=length(z);
    ME1=zeros(Nz,3);

    %% Compute First Echo %%


    %% From first RF to second RF
    %% alpha_1 -------*alpha_2
    seq=seq.generate_rx(rf_1);
    B1=rf_1/gyro/taup_1;
    for n=1:Nz
        [ME1(n,1),ME1(n,2),ME1(n,3)]= ...
            rf_pulse(Meqm(1),Meqm(2),Meqm(3),B1,seq.Gz*z(n),taup_1,gyro);
        seq=seq.generate_rz(phi_1*z(n));
        ME1(n,:)=seq.rotate_z(ME1(n,:)');
    end

    %% Decay of transverse magnetization
    ME1(:,1:2)=ME1(:,1:2)*exp(-t1/T2);

    %% Decay and rebirth of longtidunal magnetization
    ME1(:,3)=ME1(:,3)*exp(-t1/T1)+(1-exp(-t1/T1));


    %% From second RF to spin echo
    %% alpha_1 -------alpha_2-------*Spin Echo

    seq=seq.generate_rx(rf_2);
    for n=1:Nz
        [ME1(n,1),ME1(n,2),ME1(n,3)]= ...
            rf_pulse(ME1(n,1),ME1(n,2),ME1(n,3),B1,seq.Gz*z(n),taup_2,gyro);
        seq=seq.generate_rz(phi_1*z(n));
        ME1(n,:)=seq.rotate_z(ME1(n,:)');
    end

    %% Decay of transverse magnetization
    ME1(:,1:2)=ME1(:,1:2)*exp(-t1/T2);

    %% Decay and rebirth of longtidunal magnetization
    ME1(:,3)=ME1(:,3)*exp(-t1/T1)+(1-exp(-t1/T1));


    figure(101)
    plot(z,abs(complex(ME1(:,1),ME1(:,2))))
    xlabel('z')
    ylabel('Mxy at spin echo')

    line_1_simulated=zeros(Nsamples,3);
    for it=1:Nsamples
        Mit=ME1;
        for n=1:Nz
            seq=seq.generate_rz((it-1-Nsamples/2)*sampling_time*seq.Gz* z(n)*gyro);
            Mit(n,:)=seq.rotate_z(Mit(n,:)');
        end
        line_1_simulated(it,:)=sum(rho_z.*Mit);
    end

    figure(5)
    plot(t_sampled, sqrt(line_1_simulated(:,1).^2+line_1_simulated(:,2).^2));
    xlabel('time')
    ylabel('Signal (ABS)')
    title('Simulated Data: SE')

    %% From spin echo to third RF pulse
    %% alpha_1 -------alpha_2-------Spin Echo----------------*alpha_3



    ME2=ME1;

    %% Decay of transverse magnetization
    ME2(:,1:2)=ME2(:,1:2)*exp(-(t2-t1)/T2);

    %% Decay and rebirth of longtidunal magnetization
    ME2(:,3)=ME2(:,3)*exp(-(t2-t1)/T1)+(1-exp(-(t2-t1)/T1));
    seq=seq.generate_rx(rf_3);
    for n=1:Nz
        seq=seq.generate_rz((phi_2-phi_1)*z(n));
        ME2(n,:)=seq.rotate_z(ME2(n,:)');
        [ME2(n,1),ME2(n,2),ME2(n,3)]= ...
            rf_pulse(ME2(1),ME2(2),ME2(3),B1,seq.Gz*z(n),taup_3,gyro);
    end




    %% From third RF pulse to STE
    %% alpha_1 -------alpha_2-------Spin Echo----------------alpha_3-------*STE

    for n=1:Nz
        seq=seq.generate_rz(phi_1*z(n));
        ME2(n,:)=seq.rotate_z(ME2(n,:)');
    end

    %% Decay of transverse magnetization
    ME2(:,1:2)=ME2(:,1:2)*exp(-t1/T2);

    %% Decay and rebirth of longtidunal magnetization
    ME2(:,3)=ME2(:,3)*exp(-t1/T1)+(1-exp(-t1/T1));

    figure(102)
    plot(z,abs(complex(ME2(:,1),ME2(:,2))))
    xlabel('z')
    ylabel('Mxy at stimulated echo')


    line_2_simulated=zeros(Nsamples,3);
    for it=1:Nsamples
        Mit=ME2;
        for n=1:Nz
            seq=seq.generate_rz((it-1-Nsamples/2)*sampling_time*seq.Gz*z(n)*gyro);
            Mit(n,:)=seq.rotate_z(Mit(n,:)');
        end
        line_2_simulated(it,:)=sum(Mit.*rho_z);
    end

    figure(55)
    hold on
    plot(t_sampled, sqrt(line_2_simulated(:,1).^2+line_2_simulated(:,2).^2),'r');
    xlabel('time')
    ylabel('Signal (ABS)')
    title('Simulated Data: STE')


    figure(6)
    hold on
    plot(t_sampled, sqrt(line_2_simulated(:,1).^2+line_2_simulated(:,2).^2),'*');
    xlabel('time')
    ylabel('Signal (ABS)')
    title('Simulated Data STE')


    figure(7)
    hold on
    plot(t_sampled, sqrt(line_1_simulated(:,1).^2+line_1_simulated(:,2).^2),'*');
    xlabel('time')
    ylabel('Signal (ABS)')
    title('Simulated Data SE')



    expected_ratio=0.5*sin(rf_2)*sin(rf_3)/sin(rf_2/2)/sin(rf_2/2);

    Sx1=sum(ME1(:,1))/Nz;
    Sy1=sum(ME1(:,2))/Nz;
    Sx2=sum(ME2(:,1))/Nz;
    Sy2=sum(ME2(:,2))/Nz;

    simulated_ratio=sqrt(Sx2*Sx2+Sy2*Sy2)/sqrt(Sx1*Sx1+Sy1*Sy1);


    %% Analyze/Visualize Slice Profiles

    figure(8)
    clf
    Mxy1= ((sqrt(ME1(:,1).^2+ME1(:,2).^2)));
    plot(z,Mxy1)
    xlabel('z')
    ylabel('Mxy(z)')
    title('Transverse M at SE')

    figure(9)
    clf
    Mxy2= ((sqrt(ME2(:,1).^2+ME2(:,2).^2)));
    plot(z,Mxy2)
    xlabel('z')
    ylabel('Mxy(z)')
    title('Transverse M at SE')



    figure(10)
    clf
    Mxy1= rho_z.*((sqrt(ME1(:,1).^2+ME1(:,2).^2)));
    plot(z,Mxy1)
    xlabel('z')
    ylabel('rho(z)*Mxy(z)')
    title('Transverse M x rho(z) at SE')


    figure(11)
    clf
    Mxy2= rho_z.*((sqrt(ME2(:,1).^2+ME2(:,2).^2)));
    plot(z,Mxy2)
    xlabel('z')
    ylabel('rho(z)*Mxy(z)')
    title('Transverse M x rho(z) at STE')


    %% Fourrier Analysis
    line_1=zeros(Nsamples,2);
    line_2=zeros(Nsamples,2);
    for i=1:Nsamples
        line_1(i,:)= exp(1i*pi*(i-1))*line_1_simulated(i,1:2);
        line_2(i,:)= exp(1i*pi*(i-1))*line_2_simulated(i,1:2);
    end

    rho_1=ifft(line_1(:,1)+1i*line_1(:,2));
    rho_2=ifft(line_2(:,1)+1i*line_2(:,2));

    rho_1_recentered=zeros(Nsamples,1);
    rho_2_recentered=zeros(Nsamples,1);

    for i=1:Nsamples
        rho_1_recentered(i)= exp(1i*pi*(i-1))*rho_1(i);
        rho_2_recentered(i)= exp(1i*pi*(i-1))*rho_2(i);
    end

    figure(18)
    plot(abs(rho_1_recentered))
    figure(36)
    plot(abs(rho_2_recentered))
end