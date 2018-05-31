function [pMarg,pJoint]=forwbackJOINT(P,phit);

%
% Forward - backward calculation.
% Marginal - pMarg p(x_i=k|data)
% Joint - pJoint p(x_i=k, x_{i+1}=l | data) 
%

% Input P is the Markov transition probabilities
% inputs phit are exponents in Normal probabilities of data

% For large time-series it is more stable to work with logarithms.
phit=exp(phit);
% See also subfile marglikl.

% FORWARD
d=size(P,1);
tot=size(phit,1);
probp=zeros(tot,d);
pii=P^100;
probp(1,1:d)=pii(1,:);
for t=1:tot,
    % Filtering probability, non-normalized
    probuun(t,1:d)=phit(t,1:d).*probp(t,1:d);
    % Filtering probability, normalized
    probu(t,1:d)=(1/sum(probuun(t,1:d)))*probuun(t,1:d);
    % Predictive probability, normalized
    probp(t+1,1:d)=sum(P(1:d,1:d).*(probu(t,1:d)'*ones(1,d)),1);
    
    % Joint prob x_{t-1}, x_t, and y_t, given y_1,...y_{t-1} 
    if (t>1)
        psmmat(1:d,1:d,t)=(ones(d,1)*phit(t,1:d)).*P(1:d,1:d).*...
	    (probu(t-1,1:d)'*ones(1,d));
        % Backward transition probability (re-used going back)
        Pback(1:d,1:d,t)=psmmat(1:d,1:d,t)./(ones(d,1)*sum(psmmat(1:d,1:d,t),1));
    else
        psmmat(1:d,1:d,t)=(ones(d,1)*phit(t,1:d)).*P(1:d,1:d).*...
	    (probp(t,1:d)'*ones(1,d));
        % Backward transition probability (re-used going back)
        Pback(1:d,1:d,t)=psmmat(1:d,1:d,t)./(ones(d,1)*sum(psmmat(1:d,1:d,t),1));
   end;
    
end;

% MARGINAL PROBS BY BACKWARD EVALUATION
pMarg=zeros(tot,d);
pJoint=zeros(d,d,tot);
vv=cumsum(probu(tot,1:d));
pMarg(tot,1:d)=vv-[0 vv(1:d-1)];
for s=1:(tot-1),
        pJoint(1:d,1:d,tot-s)=Pback(1:d,1:d,tot-s+1).*(ones(d,1)*pMarg(tot-s+1,1:d));
        pvec=sum(pJoint(1:d,1:d,tot-s),2);        
        vv=cumsum(pvec)';
        pMarg(tot-s,1:d)=vv-[0 vv(1:d-1)];
end;

