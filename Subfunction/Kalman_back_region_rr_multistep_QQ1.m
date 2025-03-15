function [xr00,Pr0]=Kalman_back_region_rr_multistep_QQ1(time,ym,t,sig2,QQg,q,X00,P00)

% Input data:
% time: the study months in year
% y: the observations involving SHCs (monthly SHC numbers * months)
% t: the correlation time of FOGM
% sig: the variance of FOGM
% q: the optimal order for Markov process
% x00: initial parameter (usually obtain from least-square)
% P00: initial covariance matrix of parameter (usually obtain from least-square)

% Output data:
% xr00: Irregular parameters
% Pr00: The covariance matrices of xr


% Converts geoid coefficients (gc) to mass coefficients (mc)


%Initialize前s个
M=size(ym,1);% The SHCs numbers
IM=eye(M);
Q=diag(ones(1,M),0)*0.000; % Process noise
I=diag(ones(1,M),0);
t0=time(1);
% X01s=zeros(M,1); %预测值          
% P01s=eye(M)*10.^(14); %预测值

%Design Matrix
H=IM;
EE=QQg;
% Update
for i=1:q
Pr00(:,:,i)=P00(:,:,i);
xr00(:,i)=X00(:,i);
end

% The rest n-s months
for i=q+1:length(time)
   
   %过程噪声
   Qr0=zeros(M,M);
   for jk=i-q:i-1
   dt=-(time(i)-time(jk))*365.25;    
   Qr0=Qr0+eye(M)*sig2*exp(-2*dt/t);
   end
   Qr(:,:,i)=eye(M)*sig2-Qr0;
   
   %预测
    X01s=zeros(M,1);
    P01s0=zeros(M,M);
   for jk=i-q:i-1
   dt=-(time(i)-time(jk))*365.25;    
   B=eye(M)*diag(exp(-dt/t));    
   X01s=X01s+B*xr00(:,jk);
   P01s0=P01s0+B*Pr00(:,:,jk)*B';
   end
   P01s=P01s0+Qr(:,:,i);
   
   %更新
   K=P01s*H'*inv(H*P01s*H'+EE(:,:,i));
   %Update
   vv(:,i)=ym(:,i)-H*X01s;
   xxf=X01s+K*vv(:,i);
   PPf=(I-K*H)*P01s;
  
   Pr00(:,:,i)=PPf;
   xr00(:,i)=xxf;
end
Pr0(:,:,1)=Pr00(:,:,length(time));%n个月份
Pr0(:,:,2)=Pr00(:,:,length(time)-1);%n-1个月份
%滤波

end
