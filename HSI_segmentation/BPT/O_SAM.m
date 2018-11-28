function dist = O_SAM(Ri,Rj)

% Computes the spectral angle between the two spectra Ri and Rj

dist = acos(sum(Ri.*Rj)/(norm(Ri)*norm(Rj)));