function hmm_em;

% EM algorithm to estimate model parameters in a hidden Markov model. 
% markov chain, x_i \in {0,...,d}, i=1,...,n.
% p(x_i=l|x_{i-1}=k,x_{i-2}=k2,....)=p(x_i=l|x_{i-1}=k)
% likelihood
% Data p(y_j|x_j=k)=mu_k+N(0,tau^2), j=1,...,n%% 
% Separate subfiles.
%
% This is done for the mine borehole data.
num_hard1 = csvread('C:\Users\veenas\Desktop\Latex-hmm\Codelibrary\HMMdata.csv', 1,1 );
%%%selecting the MWD data variables Penetration rotation and dampening
y=[num_hard1(:,2),num_hard1(:,6), num_hard1(:,7)];
n=size(y,1);
% Depth plots
figure(1); 
subplot(1,3,1)
plot(num_hard1(:,2),num_hard1(:,1),'k'); axis ij;
title('');
 xlabel('Penrate'); ylabel('Depth');
 
 figure(1); 
 subplot(1,3,2)
plot(num_hard1(:,6),num_hard1(:,1),'k'); axis ij;
title('');
 xlabel('Rot Press'); ylabel('Depth');
 
 figure(1); 
 subplot(1,3,3)
 plot(num_hard1(:,7),num_hard1(:,1),'k'); axis ij;
title('');
 xlabel('Damp Press'); ylabel('Depth');
%%%%%
%initializing parameters

 d=3;
 P=[0.7 0.2 0.1;0.15 0.7 0.15;0.1 0.1 0.8];
 mu=[2.5 60 71;2 57 74;1.6 56 75];
 sigma=[0.3 -0.17 0.17 ; -0.17 70 -17; 0.17 -17 120];

sigma=[cov(y)];
iter=100;
kk=1;
while (kk<iter),
    
    %% Expectation step    
    % Likelihood of data 
    phit=zeros(n,d);
    for t=1:n,
        for i=1:d,
            phit(t,i)=-log(det(sigma))-0.5*(y(t,:)-mu(i,:))*inv(sigma)*(y(t,:)-mu(i,:))';
        end;
    end;
    % Posterior computation by forward-backward routine
     [pMarg,pJoint]=forwbackJOINT(P,phit);  
    
    %% Maximization step
  [P,mu,sigma]=maxlikl_JOINTpar_multi(y,pMarg,pJoint,d);    

  %%%Save the estimated parameter value for each EM steps
    Pit(:,:,:,kk)=P;
    P;
    mit(:,:,kk)=mu;
    mu;
    sit(:,:,kk)=sigma;
    sigma;
    kk=kk+1;
end

mu
sigma
%%%% checking for convergence of EM estimated parameters
figure(4);
clf;
cc=1;
for ii=1:d,
    for jj=1:d,
        subplot(d,d,cc), plot(shiftdim(Pit(ii,jj,:),1));
        cc=cc+1;
    end
end
suptitle('Transition matrix parameters, EM steps');

figure(5);
clf;
cc=1;
for ii=1:d,
    for jj=1:d,
        subplot(d,d,cc), plot(shiftdim(mit(ii,jj,:),1));
        cc=cc+1;
    end
end
suptitle('Mean parameters, EM steps');

figure(6);
cc=1;
for ii=1:d,
    for jj=1:d,
        subplot(d,d,cc), plot(shiftdim(sit(ii,jj,:),1));        
        cc=cc+1;
    end
end
suptitle('covariance matrix in EM steps')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Calculate the MAP using Virterbi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawreal=25;
[pMarg,vitX,MMAP,sampV]=forwbackVIT(P,phit,drawreal);
MAP4=vitX;
% %1939-21 for uniform coloring in pic
% for ss=1:length(vitX)
%     if vitX(ss)==1
%         MAP4(ss)=2;
%     elseif vitX(ss)==2
%          MAP4(ss)=3;
%     else vitX(ss)==3
%         MAP4(ss)=1;
%     end
% end
figure(13);
subplot(1,8,1)
a=num_hard1(:,1);
x=[1:25];
imagesc(x,a,MAP4*ones(1,1)); colormap(flipud(gray));
ylabel('Depth');
xlabel('MAP-PRD');
set(gca,'xtick',[])
%  xlswrite('C:\Users\veenas\Desktop\Latex-hmm\Codelibrary\MAP7s.xls',vitX)
grid on,

%% ------------------

