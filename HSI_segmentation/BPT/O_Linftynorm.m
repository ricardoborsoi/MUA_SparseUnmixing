function dist = O_Linftynorm(Ri,Rj)

% computes the L-infinity norm between the two spectra Ri and Rj

dist = max(abs(Ri-Rj));