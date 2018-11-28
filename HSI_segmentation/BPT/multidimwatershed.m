function S = multidimwatershed(im,method)

% S = multidimwatershed(IM,method)
% Performs on the image IM one of the following multidimentional watershed
% segmentation by performing a watershed transformation of the image
% gradient, the gradient being defined according to the string 'method':
% - method = 'supremum' => performs morphological supremum gradient
% (structuring element is by default a 3x3 square).
% method = 'metric' => performs metric-based gradient (metric distance is
% by default the Euclidean distance).
% method = 'robust' => performs the robust color morphological gradient
% (the value retained by default is the second highest in the
% 8-neighborhood).
%--------------------------------------------------------------------------
% The output is the structure S with the following fields:
% - 'gradient' contains the multidimensional gradient of the input image.
% - 'nowhed' contains the watershed segmentation with border pixels.
% - 'whed' contains the watersged segmantion with border pixels removed

if isequal(method,'supremum')==1
    [grad,mapnowhed] = supgrad(im);
elseif isequal(method,'metric')==1
    [grad,mapnowhed] = mbgrad(im);
elseif isequal(method,'robust')==1
    [grad,mapnowhed] = rcmgrad(im);
else
    error('incorrect input method, choose between supremum, metric or robust')
end

im = reshape(im,size(im,1)*size(im,2),size(im,3));
mapwhed = whed(mapnowhed,im);

S = struct;
S.gradient = grad;
S.nowhed = mapnowhed;
S.whed = mapwhed;
