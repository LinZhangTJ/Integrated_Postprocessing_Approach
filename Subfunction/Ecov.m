function C=Ecov(lat,lon,x)
%Input data:
%lat:latitude
%lon:longtitude
%x:parameter
%QQ noise covariance matrix of parameter

%Output data
%C:signal covariance matrix of parameter

lat1=lat*pi/180;
lon1=lon*pi/180;
R=6371;
n=size(lat,1);
for i=1:n
    for j=1:n
       %Spherical distance
       a = sin((lat1(j) - lat1(i))/2)^2 + cos(lat1(i)) * cos(lat1(j)) * sin((lon1(j) - lon1(i))/2)^2;
       c = 2 * atan2(sqrt(a), sqrt(1-a));
       S(i,j) = R * c;
      
       G(i,j)=x(i)*x(j);
    end
end

SG(:,1)=reshape(S,n*n,1);
SG(:,2)=reshape(G,n*n,1);

SG=sortrows(SG, 1);
SG=[SG;zeros(1,2)];
ss=zeros(1,2);
num=0;
sn=0;
for i=1:n*n
    if SG(i,1)==SG(i+1,1)
        ss=ss+SG(i,:);
        sn=sn+1;
    elseif SG(i,1)~=SG(i+1,1)
        num=num+1;
        sn=sn+1;
        ss=ss+SG(i,:);
        SG1(num,:)=ss;
        ssn(num)=sn;
        ss=zeros(1,2);
        sn=0;
    end
end

SG2=SG1./ssn';
SG3= SG2(SG2(:, 2) > 0, :);

%c*exp(-d/D)
y=log(SG3(:,2));
A=[ones(size(y,1),1) -SG3(:,1)];
x=inv(A'*A)*A'*y;
k1=exp(x(1));
k2=1/x(2);

C=k1*exp(-S/k2);

end



     