function dist = O_L1norm(Ri,Rj)

% computes the L1 norm (Manhattan norm) bewteen the two spectra Ri and Rj

dist = sum(abs(Ri-Rj));