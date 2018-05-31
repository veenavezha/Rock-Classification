     function [P,mu,sigma]=maxlikl_JOINTpar_multi(y,pMarg,pJoint,d);

%
% data y
% pMarg p(x_i=k|data)
% pJoint p(x_i=k, x_{i+1}=l | data) 
% number of classes d

[n,m]=size(y);
% n is length of series (number of depths)
% m is number of multivariate measurements per depth

% Estimate P
C=sum(pJoint,3);
P=C./(sum(pMarg,1)'*ones(1,d));

% Estimate mu and sigma
lev=zeros(n,1);

differ=zeros(m,n,d);
for ii=1:d,
    mu(ii,1:m)=sum(pMarg(:,ii).*y)/sum(pMarg(:,ii)); 
    differ(1:m,1:n,ii)=(y-ones(n,1)*mu(ii,1:m))';
end
    
sigma=zeros(m,m);
for tt=1:n,
    for ii=1:d,
        sigma=sigma+pMarg(tt,ii)*(differ(:,tt,ii)*differ(:,tt,ii)');
    end
end
sigma=sigma/n;


