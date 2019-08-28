function [SM1_UKF,MM_UKF]=Localization_tracking_TOA_3D_6state_testing_r(distance,T,S);
%   S1 = [0.01;0.01;0];
%   S2 = [200;0;0];
%   S3 = [200;0;250];
%   S4=[0;0;250];

%   S1 = [178;1061;614];
%   S2 =[378;1065;618];
%   S3 = [375;1061;866];
%   S4=[178;1058;867];
  
  
  sd =11;
%   sd =5;
%   S=[S1 S2 S3 S4];
  
  qx = 500;
  qy = 800;
  qz =700;
%   qx = 1000;
%   qy = 2000;
%   qz =1500;
R = sd^2;

  h_func = @bot_h3D;
  
  Y=distance;
 
  F = [0 0 0 1 0 0;
       0 0 0 0 1 0;
       0 0 0 0 0 1;
       0 0 0 0 0 0;
       0 0 0 0 0 0;
       0 0 0 0 0 0];
   L= [0 0 0;
       0 0 0;
       0 0 0;
       1 0 0;
       0 1 0;
       0 0 1];

    dt =0.04;
  [A,Q] = lti_disc(F,L,diag([qx^2 qy^2 qz^2]),dt);
 
  M_0 = [0;0;1000;0;0;0];
  P_0 = diag([1 1 1 10 10 10]);
  M = M_0;
  P = P_0;
  % Filter with UKF
  for k=1:size(Y,2)
    [M,P] = ukf_predict1(M,P,A,Q);
    [M,P] = ukf_update1(M,P,Y(:,k),h_func,R*eye(size(Y,1)),S);
    MM_UKF(:,k)   = M;
    PP_UKF(:,:,k) = P;
%     ME_UKF(k) = P(1,1) + P(2,2) + P(3,3);
  end
   
  % URTS Smoother   
  [SM1_UKF,SP1_UKF] = urts_smooth1(MM_UKF,PP_UKF,A,Q);
%   ME1_UKF = squeeze(SP1_UKF(1,1,:)+SP1_UKF(2,2,:)+SP1_UKF(3,3,:));
  

