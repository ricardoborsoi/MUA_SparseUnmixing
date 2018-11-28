function dist = O_Battacharyya(Hi,Hj)

% computes the Battacharyya distance bewteen the two sets of histograms Hi and
% Hj

BCcoeff = sum(sqrt(Hi).*sqrt(Hj),2);
if any(BCcoeff==0)
    BCcoeff(BCcoeff==0) = 1e-10;
end
dist = sum(-log10(BCcoeff));

