function dist = O_KLDiv(Ri,Rj)

% computes the Kullback-Leibler divergence between the two spectra Ri and
% Rj

P_Ri = Ri/sum(Ri);
P_Rj = Rj/sum(Rj);

dist = 0.5*Kullback_Leibler(P_Ri,P_Rj) + 0.5*Kullback_Leibler(P_Rj,P_Ri);

function KL_PQ = Kullback_Leibler(P,Q)

if length(P)~=length(Q)
    error('probability distributions have not the same length. Please check')
end

Q_0 = Q<=0;
P(Q_0) = 0; % Q = 0 enforces P = 0 (absolute continuity)
P_0 = P<=0;
P(P_0) = []; Q(P_0) = []; % because x*log(x) => 0 when x => 0

L = log(P./Q);
KL_PQ = sum(P.*L);