%analyze the gait data acquired
close all
 clearvars -except data timeStamps subject
fs=125000;Ts=1/fs;tF=7e-3;t=Ts:Ts:tF;N=tF*fs;
% Wnd=rectwin(N)';
f0=39e3;f1=41e3;

tx1=5*chirp(t,f0,tF,f1,'linear');
tx2=5*chirp(t,f1,tF,f0,'linear');
% output1=Wnd.*tx1;
% output2=Wnd.*tx2;
outputSignal1 =tx1;
outputSignal2 =tx2;
Temp=23;
v=331.5+0.6*Temp;
Range1=[];
Range2=[];
Range_ephp1=[];
Range_ephp2=[];
% l=3750;
l=5000;
m=l/2;
n=length(data)/m;
% Finding the correlation of the received signals with the stored template signals
for i=1:n
C11=xcorr(data(1+m*(i-1):i*m,1)-mean(data(1+m*(i-1):i*m,1)),outputSignal1);
C12=xcorr(data(1+m*(i-1):i*m,2)-mean(data(1+m*(i-1):i*m,2)),outputSignal1);
C13=xcorr(data(1+m*(i-1):i*m,3)-mean(data(1+m*(i-1):i*m,3)),outputSignal1);
C14=xcorr(data(1+m*(i-1):i*m,4)-mean(data(1+m*(i-1):i*m,4)),outputSignal1);
C21=xcorr(data(1+m*(i-1):i*m,1)-mean(data(1+m*(i-1):i*m,1)),outputSignal2);
C22=xcorr(data(1+m*(i-1):i*m,2)-mean(data(1+m*(i-1):i*m,2)),outputSignal2);
C23=xcorr(data(1+m*(i-1):i*m,3)-mean(data(1+m*(i-1):i*m,3)),outputSignal2);
C24=xcorr(data(1+m*(i-1):i*m,4)-mean(data(1+m*(i-1):i*m,4)),outputSignal2);

[~,p11(i)]=max(envelope(C11));
[~,p12(i)]=max(envelope(C12));
[~,p13(i)]=max(envelope(C13));
[~,p14(i)]=max(envelope(C14));
[~,p21(i)]=max(envelope(C21));
[~,p22(i)]=max(envelope(C22));
[~,p23(i)]=max(envelope(C23));
[~,p24(i)]=max(envelope(C24));
% p11=p11+2700;
% p12=p12+2700;
% p13=p13+2700;
% p14=p14+2700;
% p21=p21+2700;
% p22=p22+2700;
% p23=p23+2700;
% p24=p24+2700;
%multipath compensation
% [p_ephp11(i)]=multipath(C11,l);
% [p_ephp12(i)]=multipath(C12,l);
% [p_ephp13(i)]=multipath(C13,l);
% [p_ephp14(i)]=multipath(C14,l);
% [p_ephp21(i)]=multipath(C21,l);
% [p_ephp22(i)]=multipath(C22,l);
% [p_ephp23(i)]=multipath(C23,l);
% [p_ephp24(i)]=multipath(C24,l);
if ~mod(i,2)
TOF11(i/2)=(p21(i)+p11(i-1)-2*length(outputSignal1))/2;
TOF12(i/2)=(p22(i)+p12(i-1)-2*length(outputSignal1))/2;
TOF13(i/2)=(p23(i)+p13(i-1)-2*length(outputSignal1))/2;
TOF14(i/2)=(p24(i)+p14(i-1)-2*length(outputSignal1))/2;
TOF21(i/2)=(p11(i)+p21(i-1)-2*length(outputSignal2))/2;
TOF22(i/2)=(p12(i)+p22(i-1)-2*length(outputSignal2))/2;
TOF23(i/2)=(p13(i)+p23(i-1)-2*length(outputSignal2))/2;
TOF24(i/2)=(p14(i)+p24(i-1)-2*length(outputSignal2))/2;
%ephp
% TOF11_ephp(i/2)=(p_ephp21(i)+p_ephp11(i-1))/2;
% TOF12_ephp(i/2)=(p_ephp22(i)+p_ephp12(i-1))/2;
% TOF13_ephp(i/2)=(p_ephp23(i)+p_ephp13(i-1))/2;
% TOF14_ephp(i/2)=(p_ephp24(i)+p_ephp14(i-1))/2;
% TOF21(i/2)=(p_ephp11(i)+p_ephp21(i-1))/2;
% TOF22(i/2)=(p_ephp12(i)+p_ephp22(i-1))/2;
% TOF23(i/2)=(p_ephp13(i)+p_ephp23(i-1))/2;
% TOF24(i/2)=(p_ephp14(i)+p_ephp24(i-1))/2;

% p11old=p11+p11old-75;
% p12old=p12+p12old-75;
% p13old=p13+p13old-75;
% p14old=p14+p14old-75;
% p21old=p21+p21old-75;
% p22old=p22+p22old-75;
% p23old=p23+p23old-75;
% p24old=p24+p24old-75;

ran1(1,1)=TOF11(i/2)*v/fs;
ran1(2,1)=TOF12(i/2)*v/fs;
ran1(3,1)=TOF13(i/2)*v/fs;
ran1(4,1)=TOF14(i/2)*v/fs;
ran2(1,1)=TOF21(i/2)*v/fs;
ran2(2,1)=TOF22(i/2)*v/fs;
ran2(3,1)=TOF23(i/2)*v/fs;
ran2(4,1)=TOF24(i/2)*v/fs;

% ranephp1(1,1)=TOF_ephp11(i)*v/fs;
% ranephp1(2,1)=TOF_ephp12(i)*v/fs;
% ranephp1(3,1)=TOF_ephp13(i)*v/fs;
% ranephp1(4,1)=TOF_ephp14(i)*v/fs;
% ranephp2(1,1)=TOF_ephp21(i)*v/fs;
% ranephp2(2,1)=TOF_ephp22(i)*v/fs;
% ranephp2(3,1)=TOF_ephp23(i)*v/fs;
% ranephp2(4,1)=TOF_ephp24(i)*v/fs;

Range1=[Range1 ran1];
Range2=[Range2 ran2];
% Range_ephp1=[Range_ephp1 ranephp1];
% Range_ephp2=[Range_ephp2 ranephp2];
time(1,i/2)=timeStamps(l*(i/2-1)+1);

end
end
figure
subplot(2,4,1)
plot(Range1(1,:))
subplot(2,4,2)
plot(Range1(2,:))
subplot(2,4,3)
plot(Range1(3,:))
subplot(2,4,4)
plot(Range1(4,:))
subplot(2,4,5)
plot(Range2(1,:))
subplot(2,4,6)
plot(Range2(2,:))
subplot(2,4,7)
plot(Range2(3,:))
subplot(2,4,8)
plot(Range2(4,:))
% figure
% subplot(2,4,1)
% plot(Range_ephp1(1,:))
% subplot(2,4,2)
% plot(Range_ephp1(2,:))
% subplot(2,4,3)
% plot(Range_ephp1(3,:))
% subplot(2,4,4)
% plot(Range_ephp1(4,:))
% subplot(2,4,5)
% plot(Range_ephp2(1,:))
% subplot(2,4,6)
% plot(Range_ephp2(2,:))
% subplot(2,4,7)
% plot(Range_ephp2(3,:))
% subplot(2,4,8)
% plot(Range_ephp2(4,:))

 [cd11,cd12,cd13,cd14,cd21,cd22,cd23,cd24,ctime,~,~,~,~,S]=analyzecamdata1(subject,v);
%  Range11=(Range1-[4.4102;4.3330;4.3305;4.4256])*1000;
%  Range22=(Range2-[3.6329;3.7088;3.7055;3.6385])*1000;
%  for i=1:length(Range1)
%  Range11(:,i)=(Range1(:,i)-[4.4102;4.3330;4.3305;4.4256])*1000;
%  Range22(:,i)=(Range2(:,i)-[3.6329;3.7088;3.7055;3.6385])*1000;
%   end
 for i=1:length(Range1)
 Range11(:,i)=(Range1(:,i)-[4.6900;4.6900;4.6900;4.6900])*1000;
 Range22(:,i)=(Range2(:,i)-[4.6900;4.6900;4.6900;4.6900])*1000;
  end
% for i=1:length(Range1)
%  Range11(:,i)=(Range1(:,i)*1000-(mean(Range1'*1000)-mean([cd11;cd12;cd13;cd14]'))');
%  Range22(:,i)=(Range2(:,i)*1000-(mean(Range2'*1000)-mean([cd11;cd12;cd13;cd14]'))');
% end
%  Range_ephp11=(Range_ephp1-[0.5047;0.4144;0.4201;0.5057])*1000;
%  Range_ephp22=(Range_ephp2-[-0.2896;-0.2157;-0.2169;-.2850])*1000;
%  for i=2:length(Range11m)
%  delv(:,i)=((Range11m(:,i)-Range11m(:,i-1))/(time(i)-time(i-1)));
%  end
%  B = 1/100*ones(100,1);
% out = filter(B,1,delv);
%  for i=2:length(Range11m)
%     Range11(:,i)=Range11m(:,i)+(out(i)*40000)*0.009/2000;
%  end
figure
subplot(2,4,1)
plot(time,Range11(1,:))
 hold on
plot(ctime,cd11)
subplot(2,4,2)
plot(time,Range11(2,:))
 hold on
plot(ctime,cd12)
subplot(2,4,3)
plot(time,Range11(3,:))
 hold on
plot(ctime,cd13)
subplot(2,4,4)
plot(time,Range11(4,:))
 hold on
plot(ctime,cd14)
subplot(2,4,5)
plot(time,Range22(1,:))
 hold on
plot(ctime,cd11)
subplot(2,4,6)
plot(time,Range22(2,:))
 hold on
plot(ctime,cd12)
subplot(2,4,7)
plot(time,Range22(3,:))
 hold on
plot(ctime,cd13)
subplot(2,4,8)
plot(time,Range22(4,:))
 hold on
plot(ctime,cd14)


% figure
% subplot(2,4,1)
% plot(Range_ephp11(1,:))
% subplot(2,4,2)
% plot(Range_ephp11(2,:))
% subplot(2,4,3)
% plot(Range_ephp11(3,:))
% subplot(2,4,4)
% plot(Range_ephp11(4,:))
% subplot(2,4,5)
% plot(Range_ephp22(1,:))
% subplot(2,4,6)
% plot(Range_ephp22(2,:))
% subplot(2,4,7)
% plot(Range_ephp22(3,:))
% subplot(2,4,8)
% plot(Range_ephp22(4,:))
%calibrate
%  Range11=Range1-[4.5707;4.5736;4.5550;4.5389];
%  Range22=Range2-[3.8496;3.8485;3.8461;3.8321];
% Range11=lowpassfilt(Range11');
% Range22=lowpassfilt(Range22');
% Range11=Range11';
% Range22=Range22';
%prefiltering
ul1=1300;
ll1=300;
ul2=1300;
ll2=300;
[distance1,distance2,T1,T2,ind1,ind2] = Prefiltering(Range11,Range22,time,ul1,ll1,ul2,ll2);
% distance1=Range11;
% distance2=Range22;
% T1=time;
% T2=time;
% [distance1, T1]=velocity_limit(Range11,time);
% [distance2, T2]=velocity_limit(Range22,time);
% shift=-103;
% [distance11,CD11,timeshift,time1]=align_time(distance1(1,:),T1,cd11,ctime(1:end),shift);
% [distance12,CD12,timeshift,time1]=align_time(distance1(2,:),T1,cd12,ctime(1:end),shift);
% [distance13,CD13,timeshift,time1]=align_time(distance1(3,:),T1,cd13,ctime(1:end),shift);
% [distance14,CD14,timeshift,time1]=align_time(distance1(4,:),T1,cd14,ctime(1:end),shift);
% [distance21,CD21,timeshift,time2]=align_time(distance2(1,:),T2,cd21,ctime(1:end),shift);
% [distance22,CD22,timeshift,time2]=align_time(distance2(2,:),T2,cd22,ctime(1:end),shift);
% [distance23,CD23,timeshift,time2]=align_time(distance2(3,:),T2,cd23,ctime(1:end),shift);
% [distance24,CD24,timeshift,time2]=align_time(distance2(4,:),T2,cd24,ctime(1:end),shift);
% %filtering with low pass filter
% 
% [distance11]=lowpassfilt(dist1(1,:));
% [distance12]=lowpassfilt(dist1(2,:));
% [distance13]=lowpassfilt(dist1(3,:));
% [distance14]=lowpassfilt(dist1(4,:));
% [distance21]=lowpassfilt(dist2(1,:));
% [distance22]=lowpassfilt(dist2(2,:));
% [distance23]=lowpassfilt(dist2(3,:));
% [distance24]=lowpassfilt(dist2(4,:));
% % 
% d1=[CD11; CD12; CD13 ;CD14];
% d2=[CD21 ; CD22 ; CD23; CD24];
% distance1=[distance11; distance12; distance13 ;distance14];
% distance2=[distance21 ; distance22 ; distance23; distance24];
% T1=time1;
% T2=time2;

% % distance1=[distance11;distance12:distance13:distance14];
% % distance2=[distance21;distance22:distance23:distance24];
%unscekalman filtering and smoothing
% [r_SM1_UKF,r_MM_UKF]=Localization_tracking_TOA_3D_6state(d1,T1);
[r_SM1_UKF,r_MM_UKF]=Localization_tracking_TOA_3D_6state_testing_r(distance1,T1,S);
% [SM1_EKF1_r,MM_EKF1_r,SM1_EKF2_r,MM_EKF2_r]=Localization_EKF1_2_6state(d1,T1);
% [r_SM1_UKF,Pos_Newton_EKF_r]=Localization_newton_gaussian(distance1,T1);

% [l_SM1_UKF,l_MM_UKF]=Localization_tracking_TOA_3D_6state(d2,T2);
[l_SM1_UKF,l_MM_UKF]=Localization_tracking_TOA_3D_6state_testing_l(distance2,T2,S);
% [SM1_EKF1_l,MM_EKF1_l,SM1_EKF2_l,MM_EKF2_l]=Localization_EKF1_2_6state(d2,T2);
% [l_SM1_UKF,Pos_Newton_EKF_l]=Localization_newton_gaussian(distance1,T1);
Right=r_SM1_UKF;
Left=l_SM1_UKF;
% Right=r_SM1_UKF(:,50:end);
% Left=l_SM1_UKF(:,50:end);
% T1=T1(50:end);
% T2=T2(50:end);

figure;
subplot(211);plot(T1,distance1,'-o');
subplot(212);plot(T1,Right(1:3,:),'-o');
figure;
subplot(211);plot(T2,distance2,'-o');
subplot(212);plot(T2,Left(1:3,:),'-o');

% 
% %find time shift
timeshift1=time_shift_searching(Right,subject,T1,0);
timeshift2=time_shift_searching(Left,subject,T2,1);
% 
% %adjust the time shifts and find RMSE
[RMSE_DEV_r,X_DEV_r,Y_DEV_r,Z_DEV_r,xSM1_UKF_shift_r,ySM1_UKF_shift_r,zSM1_UKF_shift_r,xcamera_r,ycamera_r,zcamera_r,time_r]...
=time_shift_ultrasound_camera_simplified(Right,subject,T1,timeshift1,0);
[RMSE_DEV_l,X_DEV_l,Y_DEV_l,Z_DEV_l,xSM1_UKF_shift_l,ySM1_UKF_shift_l,zSM1_UKF_shift_l,xcamera_l,ycamera_l,zcamera_l,time_l]...
=time_shift_ultrasound_camera_simplified(Left,subject,T2,timeshift2,1);
crop
% RMSE after cropping the unwanted part
[RMSE_final_r,X_DEV_final_r,Y_DEV_final_r,Z_DEV_final_r,coefy_r,coefz_r]=find_rmse(xcamera_r,ycamera_r,zcamera_r,xSM1_UKF_shift_r,ySM1_UKF_shift_r,zSM1_UKF_shift_r);
[RMSE_final_l,X_DEV_final_l,Y_DEV_final_l,Z_DEV_final_l,coefy_l,coefz_l]=find_rmse(xcamera_l,ycamera_l,zcamera_l,xSM1_UKF_shift_l,ySM1_UKF_shift_l,zSM1_UKF_shift_l); 


%Find spatial and temporal parameters
% mark_peak(ycamera_r,ySM1_UKF_shift_r,time_r);
%1. step length and step time
[Stancetime_camera_r,Stancetime_ultrasound_r,Swingtime_camera_r,Swingtime_ultrasound_r,IndMax_cam_r,IndMax_us_r,IndMin_cam_r,IndMin_us_r,S1_camera_r,S1_ultrasound_r,Stridetime_camera_r,Stridetime_ultrasound_r,Stridevelocity_cam_r,Stridevelocity_us_r]= swing_time(ycamera_r(1:end),ySM1_UKF_shift_r(1:end),zcamera_r(1:end),zSM1_UKF_shift_r(1:end),0.0083333);
[Stancetime_camera_l,Stancetime_ultrasound_l,Swingtime_camera_l,Swingtime_ultrasound_l,IndMax_cam_l,IndMax_us_l,IndMin_cam_l,IndMin_us_l,S1_camera_l,S1_ultrasound_l,Stridetime_camera_l,Stridetime_ultrasound_l,Stridevelocity_cam_l,Stridevelocity_us_l]= swing_time(ycamera_l(1:end),ySM1_UKF_shift_l(1:end),zcamera_l(1:end),zSM1_UKF_shift_l(1:end),0.0083333);

%bland altman plots
% [means,diffs,meanDiff,CR,linFit] = BlandAltman(S1_camera, S1_ultrasound-mean(S1_ultrasound-S1_camera), 2);
% [means,diffs,meanDiff,CR,linFit] = BlandAltman(Swingtime_camera, Swingtime_ultrasound-mean(Swingtime_ultrasound-Swingtime_camera), 2);

%Mann-whitney U test for statistical significance.

%Wilcoxon Ranksum
[p,h,stats]=ranksum((S1_camera_l),(S1_ultrasound_l))
[p,h,stats]=ranksum((S1_camera_r),(S1_ultrasound_r))
[p,h,stats]=ranksum((Stancetime_camera_l),(Stancetime_ultrasound_l))
[p,h,stats]=ranksum((Stancetime_camera_r),(Stancetime_ultrasound_r))
[p,h,stats]=ranksum((Swingtime_camera_l),(Swingtime_ultrasound_l))
[p,h,stats]=ranksum((Swingtime_camera_r),(Swingtime_ultrasound_r))
[p,h,stats]=ranksum((MFC_cam_l),(MFC_ultrasound_l))
[p,h,stats]=ranksum((MFC_cam_r),(MFC_ultrasound_r))

mean(S1_ultrasound_l-S1_camera_l)
std(S1_ultrasound_l-S1_camera_l)
mean(S1_ultrasound_r-S1_camera_r)
std(S1_ultrasound_r-S1_camera_r)
mean(Stridetime_camera_l-Stridetime_ultrasound_l)
std(Stridetime_camera_l-Stridetime_ultrasound_l)
mean(Stridetime_camera_r-Stridetime_ultrasound_r)
std(Stridetime_camera_r-Stridetime_ultrasound_r) 
%mean and std of parameters
mean(S1_camera_l)
std(S1_camera_l)
mean(S1_ultrasound_l)
std((S1_ultrasound_l))

mean(S1_camera_r)
std(S1_camera_r)
mean(S1_ultrasound_r)
std((S1_ultrasound_r))

mean(Stridetime_camera_l)
std(Stridetime_camera_l)
mean(Stridetime_ultrasound_l)
std(Stridetime_ultrasound_l)
mean(Stridetime_camera_r)
std(Stridetime_camera_r)
mean(Stridetime_ultrasound_r)
std(Stridetime_ultrasound_r)
%Bland-Altman test for the trajectory along x, y and z
%x
[means,diffs,meanDiff,CR,linFit] = BlandAltman(xcamera_l, xSM1_UKF_shift_l, 2);
[means,diffs,meanDiff,CR,linFit] = BlandAltman(xcamera_r, xSM1_UKF_shift_r, 2);
%y
[means,diffs,meanDiff,CR,linFit] = BlandAltman(ycamera_l, ySM1_UKF_shift_l, 2);
[means,diffs,meanDiff,CR,linFit] = BlandAltman(ycamera_r, ySM1_UKF_shift_r, 2);
%z
[means,diffs,meanDiff,CR,linFit] = BlandAltman(zcamera_l, zSM1_UKF_shift_l, 2);
[means,diffs,meanDiff,CR,linFit] = BlandAltman(zcamera_r, zSM1_UKF_shift_r, 2);

%Bland-Altman test for the gait parameters
[means,diffs,meanDiff,CR,linFit] = BlandAltman(S1_camera_l, S1_ultrasound_l, 2);
[means,diffs,meanDiff,CR,linFit] = BlandAltman(S1_camera_r, S1_ultrasound_r, 2);

[means,diffs,meanDiff,CR,linFit] = BlandAltman(1000*Stridetime_camera_l, 1000*Stridetime_ultrasound_l, 2);
[means,diffs,meanDiff,CR,linFit] = BlandAltman(1000*Stridetime_camera_r, 1000*Stridetime_ultrasound_r, 2);
