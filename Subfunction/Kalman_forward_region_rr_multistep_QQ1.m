function [xr00,Pr0]=Kalman_forward_region_rr_multistep_QQ1(time,ym,t,sig2,QQg,q)

% Input data:
% time: the study months in year
% y: the observations involving SHCs (monthly SHC numbers * months)
% t: the correlation time of FOGM
% sig: the variance of FOGM
% q: the optimal order for Markov process

% Output data:
% xr00: Irregular parameters
% Pr00: The covariance matrices of xr


% Converts geoid coefficients (gc) to mass coefficients (mc)


%Initialize the first s months
M=size(ym,1);% The SHCs numbers
IM=eye(M);
Q=diag(ones(1,M),0)*0.000; % Process noise
I=diag(ones(1,M),0);
t0=time(1);
X01s=zeros(M,1); %预测值          
P01s=eye(M)*10^3; %预测值

%Design Matrix
H=IM;
EE=QQg;
% Update the first q months' values
for i=1:q
K=P01s*H'*inv(H*P01s*H'+EE(:,:,i)); % Gain matrix
% Update
vv(:,i)=ym(:,i)-H*X01s;
xxf=X01s+K*vv(:,i);
PPf=(I-K*H)*P01s;
Pr00(:,:,i)=PPf;
xr00(:,i)=xxf;
end

% The rest n-s months
for i=q+1:length(time)
   
   %Process noise
   Qr0=zeros(M,M);
   for jk=i-q:i-1
   dt=(time(i)-time(jk))*365.25;    
   Qr0=Qr0+eye(M)*sig2*exp(-2*dt/t);
   end
   Qr(:,:,i)=eye(M)*sig2-Qr0;
   
   %Prediction
    X01s=zeros(M,1);
    P01s0=zeros(M,M);
   for jk=i-q:i-1
   dt=(time(i)-time(jk))*365.25;    
   B=eye(M)*diag(exp(-dt/t));    
   X01s=X01s+B*xr00(:,jk);
   P01s0=P01s0+B*Pr00(:,:,jk)*B';
   end
   P01s=P01s0+Qr(:,:,i);
   
   K=P01s*H'*inv(H*P01s*H'+EE(:,:,i));
   %Update
   vv(:,i)=ym(:,i)-H*X01s;
   xxf=X01s+K*vv(:,i);
   PPf=(I-K*H)*P01s;
  
   Pr00(:,:,i)=PPf;
   
   xr00(:,i)=xxf;
end

Pr0(:,:,1)=Pr00(:,:,length(time));
Pr0(:,:,2)=Pr00(:,:,length(time)-1);
end
