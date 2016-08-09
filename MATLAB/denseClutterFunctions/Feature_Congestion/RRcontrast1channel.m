function contrast = RRcontrast1channel(pyr, DoG_sigma)
%
% contrast = RRcontrast1channel(pyr, DoG_sigma)
% Filters a Gaussian pyramid, pyr, with a 1-channel contrast feature detector.  
% 
% input:  
% pyr is a Gaussian pyramid. It can be computed from the following routine:
%   RRgaussianPyramid(im, numLevels)
% DoG_sigma (optional) is the size of the center-surround (Difference-of-Gaussian) 
% filter used for computing the contrast. Default = 2. Refer to DoG1filter.

% Code by Ruth Rosenholtz and Zhenlan Jin
% modified by Yuanzhen Li, Sep 2004

if nargin < 2
    DoG_sigma = 2;
end

levels = length(pyr);
% pyr should be a mx1 or 1xm cell structure. It can be computed using
% 'RRgaussianPyramid'
if (prod(size(pyr)) ~= levels)
    error('Gaussian pyramid passed to RRcontrast4channel is a cell matrix -- should be a vector!');
end

contrast = cell(levels, 1);

% Here we're using the difference-of-gaussian filters. Separable. 
% Refer to routine 'DoG1filter'.
[innerG1, outerG1] = DoG1filter(round(DoG_sigma*3), DoG_sigma);

% Do contrast feature computation with these filters:
for i=1:levels,
%    disp(i);
    % DOG1 channel:
    inner = filt2(innerG1, pyr{i});
    inner = filt2(innerG1', inner);
    outer = filt2(outerG1, pyr{i});
    outer = filt2(outerG1', outer);
    tmp = inner - outer;
    contrast{i}{1} = abs(tmp); %%.^2;
end

return;