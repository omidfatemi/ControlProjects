clc
clear all
close all
% kalman filter for rectilinear motion estimation, case of penalty kick,
% Rectilinear motion is characterized with the following
% dynamics: X_dot= A*X+G*v(n),y=C*X+w(n), v(n) is the disturbance due to wind 
% w(n) is the sensor noise both assumed as zero mean gaussian noises, and (x1,x2,x3) are the ball
% coordinates, since the motion is rectilinear, acceleration =0

A = [ 0   0   0   1   0   0
    0   0   0   0   1   0
    0   0   0   0   0   1
    0   0   0   0   0   0
    0   0   0   0   0   0
    0   0   0   0   0   0];

% B = [ 0   0   0   0   0   0]';

C = [ 1   0   0   0   0   0
      0   1   0   0   0   0
      0   0   1   0   0   0     ];

a=11;               %distance to goal in x1 direction
b=7.5/2;            %half of goal width
c=2.5;               %goal height
BallArea=0.3;        %Ball area estimate
n=ceil(2*b*c/BallArea)     %number of penalty kicks=goal area/Ball area        
YY=-b+2*rand(1,n)*b;       %randomly distributed penalty kicks within goal area
ZZ=rand(1,n)*c;            %randomly distributed penalty kicks within goal area
figure
scatter(YY,ZZ)
title(  + n + " random penalty kicks " )
axis([-3.75 3.75 0 2.5])
% DETERMINING initial velocity 
% T = a/v1 final trajectery time to reach the goal
alfa=b/a % alfa*v1 = maximum allowed initial speed of the ball in y axis direction
beta=c/a % beta*v1 = maximum allowed initial speed of the ball in z axis direction
v1min=20; %minimum assumed initial speed of the ball in x axis direction
v1max=40/sqrt(1+alfa^2+beta^2);
v1=v1min+rand(1,n)*(v1max-v1min); %initial random speed in x direction in 
% the range that guarantees the penalty kick reaches the goal and is inside the goal
T=a./v1;            %final time= a/v1
v2=YY./T;           %initial random speed in y direction
v3=ZZ./T;           %initial random speed in z direction
%%
success=0;
clear saved
% clear nPlots
for i=1:n
    
    V0=[v1(:,i) v2(:,i) v3(:,i)];
    tGoal=0.2; %required time for goal keeper to estimate ball position
    h = .5*tGoal/5000;      %time steps
    t = 0:h:tGoal;          %time vector
    x_ini=[0   0   0   V0]' ; % initial state
    x=zeros(size(A,1),size(t,2));%initialize state vectors
    x(:,1)=x_ini;
    Qb=.001; %state noise covariance matices
    Rb=.001; %output noise covariance matices
    time = length(t);
    v = sqrt(Qb)*randn(6,time);
    w = sqrt(Rb)*randn(size(C,1),time);
    L=zeros(size(A));
    P=zeros(6,6,n);
    P(:,:,1)=       [    0   0   0   0   0   0
        0   0   0   0   0   0
        0   0   0   0   0   0
        0   0   0   0.0015   0   0
        0   0   0   0   0.0015   0
        0   0   0   0   0   0.0015];
    % P is  posteriori estimate covariance matrix  (a measure of the estimated accuracy of the state estimate)
    % P is the initial  goalkeeper estimate  of the ball states
    % goalkeeper knows the ball position exactly, while has a initial
    % estimate of the ball speed
    G =eye(size(A,1));   %determines how the state disturbance is mapped into systwem states
    Q=1e-3*eye(size(A)); %E(ww')=Q, output noise covariance matrix?
    R=1e-5*eye(3) ;      %E(vv')=R,  state error covariance matrix
    xhat=zeros(size(x)); %initial estimate of system states
    
    %%
    for k=1:numel(t)-1
        y(:,k)  =   C*x(:,k)+ w(:,k); %compute system output with noise
        xdot=   A*x(:,k)+ v(:,k);      %compute system states with noise
        %  L  =   A*P*C'*(R+C*P*C')^-1;
        L  =   P(:,:,k)*C'*R^-1;       % compute kalman filter gain
        % P   =   (A-L*C)*P*(A-L*C)'+Q+R;
        Pdot  =   A*P(:,:,k)+P(:,:,k)*A'-P(:,:,k)*C'*R^-1*C*P(:,:,k)+Q;
        P(:,:,k+1)=P(:,:,k)+Pdot*h;
        xhatdot  = A*xhat(:,k) + L*(y(:,k) - C*xhat(:,k) ); 
        x(:,k+1)=x(:,k)+xdot*h;          %compute states 
        xhat(:,k+1)=xhat(:,k)+ xhatdot*h; %compute states estimate
    end
    if (norm(x(2:3,end)-xhat(2:3,end))<.3)
        success=success+1;
        saved(:,success)=[YY(i) ZZ(i)]';       
    end
    disp(i)
end
success/n
% less the sonsor covariance error R value(higher accuracy of the eyes of the goal keeper), the higher number of saves, which aligns with experience. 


%%
figure
scatter(YY,ZZ)
hold
scatter(saved(1,:),saved(2,:),'filled')
title(  +100*( success/n) + "%  penalty kicks were saved " )