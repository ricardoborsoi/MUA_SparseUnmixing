function GofF = GofF_MS(pixels,~)

[nbpix,~] = size(pixels);
pixels_mean = repmat(mean(pixels),nbpix,1);
pixels_diff = pixels-pixels_mean;
pixels_norm = sum(pixels_diff.^2,2); % L2 norm of each pixel difference to mean
GofF = sum(pixels_norm);