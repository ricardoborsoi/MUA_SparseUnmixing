function dist = O_JSDiv(Ri,Rj)

% computes the Jensen-Shannon divergence between two spectra Ri and Rj

P_Ri = Ri/sum(Ri);
P_Rj = Rj/sum(Rj);
M = 0.5*P_Ri + 0.5*P_Rj;
d1 = Kullback_Leibler(P_Ri,M);
d2 = Kullback_Leibler(P_Rj,M);

dist = sqrt((d1+d2)/2);

function KL_PQ = Kullback_Leibler(P,Q)

if length(P)~=length(Q)
    error('probability distributions have not the same length. Please check')
end
Q_0 = Q==0;
L = log(P./Q); L(Q_0) = 0; % to avoid the case where P(i)~=0 and Q(i)=0 (set the log to 0)
P = P(P~=0); L = L(P~=0); % remove case where P(i) = 0 (=> log(P(i)/Q(i)) = 0)
                          % that case yields 0 normally (0*log(0) = 0)
KL_PQ = sum(P.*L);