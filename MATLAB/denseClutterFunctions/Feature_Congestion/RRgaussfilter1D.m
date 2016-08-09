function kernel = RRgaussfilter1D(halfsupport, sigma, center)
% kernel = RRgaussfilter1D(halfsupport, sigma)
%   Creates a 1D gaussian filter kernel, centered at center (default=0), with pixels from
%   -halfsupport:halfsupport, and standard deviation sigma.

t = -halfsupport:halfsupport;
if nargin<3,
    center = 0;
end
kernel = exp(-(t-center).^2/(2*sigma^2));
kernel = kernel/sum(kernel);

return;