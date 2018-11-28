function dist = O_Euclidean(Ri,Rj)

% computes the Euclidean distance bewteen the two spectra Ri and Rj

dist = sqrt(sum((Ri-Rj).^2));