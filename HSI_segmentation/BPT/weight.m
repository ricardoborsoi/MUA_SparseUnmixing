function w = weight(N1,N2)

% computes the similarity weight between neighborhood N1 and N2

h = 0.5; % filtering parameter - controls the decay of exponential function

N = (N1-N2).^2;
sim = squeeze(sum(sum(N)));
w = exp(-sim/h^2);