function [color_clutter, contrast_clutter, orientation_clutter] = computeClutter(inputImage, numlevels, contrast_filt_sigma, contrast_pool_sigma, color_pool_sigma, contrast_pix, color_pix, orient_pix)

% [color_clutter, contrast_clutter, orientation_clutter] = computeClutter(inputImage, 
%       [numlevels], [contrast_filt_sigma], [contrast_pool_sigma], [color_pool_sigma], 
%       [contrast_pix], [color_pix], [orient_pix])
%
% Computes Feature Congestion clutter map(s) of an image. Returns:
%   "color_clutter", a cell structure containing the color clutter map(s)
%   "contrast_clutter", a cell structure containing the contrast clutter map(s)
%   "orientation_clutter", a cell structure containing the orientation clutter map(s)
% As for each of the three cell structure, *{1} is another cell structure containing 
%   the clutter maps at a number of scales specified by "numlevels", and *{2} is a 
%   single clutter map (same size as the input image) collapsed from all scales
%
% Inputs:
%   "inputImage" is the input image. It can be a string (file name of the image),
%       or an array (the image itself).
%   "numlevels" is the number of levels. Defaults to 3.
% For all the other input parameters, please see subroutines "colorClutter.m", 
% "contrastClutter.m", and "orientationClutter.m" for explanations.
%
% Reference: 
% Ruth Rosenholtz, Yuanzhen Li, and Lisa Nakano. "Measuring Visual Clutter". 
% Journal of Vision, 7(2), 2007. http://www.journalofvision.com/7/2/

% Ruth Rosenholtz, Yuanzhen Li, and Lisa Nakano, March 2007.

% set parameters:
if ~exist('numlevels')
    numlevels = 3;
end
% color 
if ~exist('color_pool_sigma')
    color_pool_sigma = 3;
end
if ~exist('color_pix')
    color_pix = 0;
end
% contrast 
if ~exist('contrast_filt_sigma')
    contrast_filt_sigma = 1;
end
if ~exist('contrast_pool_sigma')
    contrast_pool_sigma = 3*contrast_filt_sigma;
end
if ~exist('contrast_pix')
    contrast_pix = 0;
end
% orientation 
orient_pool_sigma = 7/2;
if ~exist('orient_pix')
    orient_pix = 0;
end

% load input image
if ischar(inputImage)
    try
        [tmp, map] = imread(inputImage);
    catch
        error(sprintf('Unable to open %s image file.', inputImage))
    end
    im = double(tmp);
else
    if isnumeric(inputImage)
        im = double(inputImage);
    end
end

% The input image, im, should be a MxNx3 matrix, in which case we assume it is an RGB image.  
%   If it's MxN, it's probably gray, and this program is not appropriate.
[m, n, d] = size(im);
if d==3,
    % we first convert it into the perceptually-based CIELab color space
    Lab = RGB2Lab(im);
else
    error('%s should be run on RGB color images.  Input image appears to be grayscale.\n', mfilename);
end
% get Gaussian pyramids (one for each of L,a,b)
RRLab(:,:,1) = RRgaussianPyramid(Lab(:,:,1), numlevels);
RRLab(:,:,2) = RRgaussianPyramid(Lab(:,:,2), numlevels);
RRLab(:,:,3) = RRgaussianPyramid(Lab(:,:,3), numlevels);

% compute the color clutter
[color_clutter_levels, color_clutter_map] = colorClutter(RRLab, numlevels, color_pool_sigma, color_pix);  
% compute the contrast clutter
[contrast_clutter_levels, contrast_clutter_map] = contrastClutter(RRLab, numlevels, contrast_filt_sigma, contrast_pool_sigma, contrast_pix);
% compute the orientation clutter
[orient_clutter_levels, orientation_clutter_map] = orientationClutter(RRLab, numlevels, orient_pix);

% output them in cell structures
color_clutter = {color_clutter_levels, color_clutter_map};
contrast_clutter = {contrast_clutter_levels, contrast_clutter_map};
orientation_clutter = {orient_clutter_levels, orientation_clutter_map};