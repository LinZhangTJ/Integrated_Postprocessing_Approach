%%A more convenient toolbox will be released later; here we will first present the main code.
clear all;
%% Setp1: Preparing the input data:
%(1) GRACE Sperical Harmonic Coefficients (SHCs): each column represents
%SHCs of each month, 
load ITSG_vv1.mat %ITSG SHCs after removing mean field 2004.01-2009.12, replacing C20 and C30 terms. 
%(2) Converts geoid coefficients (gc) to mass coefficients (mc)
load ITSG_T.mat
%(3) Covariance matrices of regional mass changes propagated from those of
%SHCs after expanding 2 degrees
load QQ_yt_buffer_degree2.mat %here we take Yangzte basin as an example
%(4) Latitude (lat), longtitude (lon), row (ll), column (cc), region range (msk), converting matrix
%from GRACE SHCs to regional grids (MM).
load cs2gridM_yt_buffer_degree2.mat

%% Step2: Iteration process to derive deterministic parameters and irregular signals.
%% Step2.1: Intializing 
n=size(vv1,2);%The number of months
load GIA_new.mat; %GIA read from ICE6G-D model
%Deriving GIA changes within region 
giam=gia1.*msk;
for i=1:n
giam2=reshape(giam(:,:,i)',180*360,1);
giam2(giam2==0)=[];
giamm(:,i)=giam2;
end
%Deriving region grids
for i=1:n
vv2=T*vv1(:,i)*100;
gg1(:,i)=MM*vv2-giamm(:,i);
end
m=size(gg1,1);%The number of regional grids
IM=eye(m);
load time_new.mat;%2002.04-2023.12
tt=time;
t0=time(1);
invR=eye(7*m);% The inverse of regularization matrix

for ite=1:3 %Gnerally three iterations is ok according to the diffirence level of two adjacent iterations.

%% Step2.2: Deriving deterministic parameters
%Deriving normalized obervation vector, weight matrix and regularization matrix
[Uc, Sc, Vc] = svd(invR);
Vc = 1/2*(Uc+Vc);
invR1=Vc*diag(sqrt(diag(Sc)))*Vc'; %R^(-1/2)
R1=Vc*diag(sqrt(diag(1./Sc)))*Vc'; %R^1/2

 for i=1:n
[UQ, SQ, VQ] = svd(QQ(:,:,i));
VQ = 1/2*(UQ+VQ);
Pg = VQ * diag(sqrt(1./diag(SQ))) * VQ'; %1/2 

A0=[1*IM (tt(i)-t0)*IM ((tt(i)-t0).^2)*IM cos(2*pi*(tt(i)-t0))*IM sin(2*pi*(tt(i)-t0))*IM cos(4*pi*(tt(i)-t0))*IM sin(4*pi*(tt(i)-t0))*IM];
yy01(:,i)=Pg*gg1(:,i);
AA01((i-1)*m+1:i*m,:)=Pg*A0*invR1;
 end
yy=reshape(yy01,m*n,1); %observation vector 

%Initalizing true value and unit weight variance
[U, S, V] = svd(AA01);
mm=size(AA01,1);
nn=size(AA01,2);
S = S(1:min(mm,nn), :);
U = U(:, 1:min(mm,nn));
x_ls = V * diag(1./diag(S'*S)) * S * U' * yy;
x_prior=x_ls; %inital true value of parameter
v_ls = yy - AA01 * x_prior;
sigma2 = v_ls'*v_ls/(mm-nn);%inital unit weight variance
xd=zeros(7*m,1);

% Deriving regularized estimates by updating the true data and unit weight variance
ee=100;
while ee>0.01
[xd,Qd,sigma2] = MMSE_Tik(AA01,yy,x_prior,sigma2);
ee=rms(xd-x_prior)/rms(x_prior);
x_prior=xd;
end

xd=R1*xd;

%Update the inverse matrix of regularization matrix 
for i=1:7
    xx=xd(m*(i-1)+1:i*m);
    Pd=Qd(m*(i-1)+1:i*m,m*(i-1)+1:i*m);
    Cc(:,:,i)=Ecov(lat,lon,xx);
end
invR=blkdiag(Cc(:,:,1),Cc(:,:,2),Cc(:,:,3),Cc(:,:,4),Cc(:,:,5),Cc(:,:,6),Cc(:,:,7));


%% Step2.3: Deriving irregular signals

for i=1:n
    QQ1(:,:,i)=QQ(:,:,n+1-i);
end

for i=1:n
A1=[1*IM (tt(i)-t0)*IM ((tt(i)-t0).^2)*IM cos(2*pi*(tt(i)-t0))*IM sin(2*pi*(tt(i)-t0))*IM cos(4*pi*(tt(i)-t0))*IM sin(4*pi*(tt(i)-t0))*IM];
yy_res(:,i)=gg1(:,i)-A1*xd(1:7*m);
end

%Note: if the irregualar signals are over-estimated, please reduce the time
%variable and variance, otherwise, decrease the time variable and variance.
t=[30;60;90;120];%time variable
sig2=[10.^3;5*10.^3;10.^4;5*10.^4;10.^5]; %variance 
q=2; %The order for Markov process

num=0;
for i=1:size(t,1)
    for j=1:size(sig2,1)
        for k=1:size(q,1)
num=num+1
t1=t(i)/q;
%Foward KF
[xrf,Prf]=Kalman_forward_region_rr_multistep_QQ1(time(1:n),yy_res,t1,sig2(j),QQ,q(k));
%Backward KF
[xrb,~]=Kalman_back_region_rr_multistep_QQ1(flip(time(1:n)),fliplr(yy_res),t1,sig2(j),QQ1,q(k),[xrf(:,n) xrf(:,n-1)],Prf);
%Average value
xr(:,:,num)=(xrf+fliplr(xrb))/2;

end
    end
end


end
