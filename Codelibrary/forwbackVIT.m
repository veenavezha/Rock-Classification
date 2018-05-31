function [pMarg,vitX,MMAP,samp]=forwbackVIT(P,phit,drawreal);

%
% Forward calculation and backward calculation giving
% Marginal probabilities pMarg
% Optimization, vitX 
% and samples samp if drawreal>0
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
deltaM=zeros(tot,d);
deltaJ=zeros(d,d,tot);
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
        psmmat=(ones(d,1)*phit(t,1:d)).*P(1:d,1:d).*...
	    (probu(t-1,1:d)'*ones(1,d));
        % Backward transition probability (re-used going back)
        Pback(1:d,1:d,t)=psmmat./(ones(d,1)*sum(psmmat,1));
        % Max trans (re-used going back)
        deltaJ(1:d,1:d,t)=(ones(d,1)*phit(t,1:d)).*P(1:d,1:d).*...
            (deltaM(t-1,1:d)'*ones(1,d));
        deltaMu=max(deltaJ(1:d,1:d,t));
        deltaM(t,1:d)=(1/sum(deltaMu))*deltaMu;
        gg=1;
    else
        deltaM(t,1:d)=probu(t,1:d);
   end;
    
end;

% MARGINAL PROBS BY BACKWARD EVALUATION
pMarg=zeros(tot,2);
vv=cumsum(probu(tot,1:d));
pMarg(tot,1:d)=vv-[0 vv(1:d-1)];
for s=1:(tot-1),
        pvec=sum(Pback(1:d,1:d,tot-s+1).*(ones(d,1)*pMarg(tot-s+1,1:d)),2);        
        vv=cumsum(pvec)';
        pMarg(tot-s,1:d)=vv-[0 vv(1:d-1)];
end;

% VITERBI OPTIMIZATION BACKWARD MAX
vitX=zeros(tot,1);
vitX(tot,1)=find(deltaM(tot,1:d)==max(deltaM(tot,1:d)));
for s=1:(tot-1),
    dvec=deltaJ(1:d,vitX(tot-s+1,1),tot-s+1);
    vitX(tot-s,1)=find(dvec==max(dvec));
end

%MMAP
MMAP=zeros(tot,1);
for s=1:(tot),
    MMAP(s,1)=find(pMarg(s,1:d)==max(pMarg(s,1:d)));    
end
% SAMPLE BACKWARD
if (drawreal>0) 
    real=drawreal;     % number of realizations
    samp=zeros(tot,real);
    vec=cumsum(probu(tot,1:d));
    for j=1:real,
        samp(tot,j)=min(find(vec>rand));
        for s=1:(tot-1),
            vec=cumsum(Pback(1:d,samp(tot-s+1,j),tot-s+1));
            samp(tot-s,j)=min(find(vec>rand));
        end;
    end;
end;
