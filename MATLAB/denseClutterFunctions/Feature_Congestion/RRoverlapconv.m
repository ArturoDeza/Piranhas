function out = RRoverlapconv(kernel, in)
% out = RRoverlapconv(kernel, in)
%   Filters the image in with filter kernel, where it only "counts" the
%   part of the filter that overlaps the image.  Rescales the filter so its
%   weights which overlap the image sum to the same as the full filter
%   kernel.

% Convolve with the original kernel
out = conv2(in, kernel, 'same');

% Convolve kernel with an image of 1's, of the same size as the input image
[m, n] = size(in);
rect = ones(m, n);
overlapsum = conv2(rect, kernel, 'same');

% Now scale the output image at each pixel by the relative overlap of the
% filter with the image
out = sum(kernel(:)) * out ./ overlapsum;

return;