function [clutter_levels, clutter_map] = orientationClutter(inputImage, numlevels, pix)
% [clutter_levels, clutter_map] = orientationClutter(inputImage, numlevels, [pix])
% 
% computes the orientation clutter map(s) of an image. Returns:
% "clutter_levels", a cell structure, containing the orientation clutter at a number 
%    of scales specified by numlevels;  -- cell(numlevels,1) --, the n'th level of 
%    which can be accessed using clutter_levels{n}{1}
% "clutter_map", an array of the same size as inputImage, is a single clutter map 
%    collapsed from clutter_levels, which is the clutter measure at multiple scales
%    now the "collapsing" is done by taking the maximal values across scales
% 
% Input:
% "inputImage" gives the input. It can be one of the following 3 things: 1. an RGB 
%    image; 2. a string, i.e., file name of an RGB image; 3, a cell structure 
%    containing Gaussian pyramids of the luminance(L) and the chrominance(a,b) 
%    channels. Refer to RRgaussianPyramid.m and computeClutter.m for how to build 
%    and use Gaussian pyramids.
% "numlevels" is the number of levels.
% Orientation clutter is computed as the "volume" of an orientation distribution
%    ellipsoid, which is the determinant of covariance matrix. Covariance 
%    matrix can be computed efficiently through linear filtering. More 
%    specifically, cov(X,Y) = E(XY)-E(X)E(Y), where E (expectation value) 
%    can be approximated by filtering with a Gaussian window. 
% "pool_sigma" is the sigma (standard deviation) of this Gaussian window, and
%    here is hard-wired to 7/2.
% 
% If "pix" = 1, displays the clutter maps.  If it's 0, does not display
% (useful for batch processing of many images).  Defaults to display.
%
% Reference:
% Ruth Rosenholtz, Yuanzhen Li, Jonathan Mansfield, and Zhenlan Jin. 
% Feature Congestion: A Measure of Display Clutter. CHI '05: Proc. of the SIGCHI conference 
% on Human factors in computing systems. May 2005. 761-770.  No orientation
% clutter in that paper, however.

% Ruth Rosenholtz, May 2006

pool_sigma = 7/2;

if ~exist('pix');
    pix = 1;
end

if iscell(inputImage)
    L_pyr = inputImage(:,:,1);
else
    if ischar(inputImage)
        try
            [tmp, map] = imread(inputImage);
        catch
            error(sprintf('Unable to open %s image file.', 'input image'))
        end
        [h, w, depth] = size(tmp);
        if prod(size(map))==0,
            im = double(tmp);
        else
            tmp = tmp(:);
            im = zeros(prod(size(tmp)), 3);
            im = map(double(tmp)+1, :);
            im = reshape(im, h, w, 3);
        end
    else
        if isnumeric(inputImage)
            im = inputImage;
        end
    end

    % The input image, im, had better be a MxNx3 matrix, in which case we assume it is an RGB image.
    %   If it's MxN, it's probably gray, and this program is not appropriate.
    [m, n, d] = size(im);
    if d==3,
        % we first convert it into the perceptually-based CIELab color space
        Lab = RGB2Lab(im);
    else
        error('%s should be run on RGB color images.  Input image appears to be grayscale.\n', mfilename);
    end
    % Get Gaussian pyramid for the luminance channel
    L_pyr = RRgaussianPyramid(Lab(:,:,1), numlevels);
end

% Compute clutter
[clutter_levels] = computeOrientationClutter(L_pyr);

% collapse over scales by taking the maximum 

% first get a Gaussian kernel to upsample the clutter maps on bigger scales
% so that the clutter maps would have the same sizes, and max can be taken
% across scales
kernel_1d = [0.05 0.25 0.4 0.25 0.05];
kernel_2d = conv2(kernel_1d, kernel_1d');
clutter_map = clutter_levels{1}{1};
for scale=2:numlevels
    clutter_here = clutter_levels{scale}{1};
    for kk = scale:-1:2
        clutter_here = upConv(clutter_here, kernel_2d, 'reflect1', [2 2]);
    end
    common_sz = min(size(clutter_map), size(clutter_here));
    clutter_map(1:common_sz(1), 1:common_sz(2)) = max(clutter_map(1:common_sz(1), 1:common_sz(2)), clutter_here(1:common_sz(1), 1:common_sz(2)));
end



% to display the clutter maps:
if (pix==1)
    min_min = min(clutter_map(:));
    max_max = max(clutter_map(:));
    figure; 
    if numlevels < 5
        nx = 2; ny = 2;
    elseif numlevels < 7
        nx = 3; ny = 2;
    elseif numlevels < 9
        nx = 4; ny = 2;
    else
        error('too many levels!!');
    end
    for scale = 1:numlevels
        subplot(ny, nx, scale); imshow(clutter_levels{scale}{1}, [min_min max_max]); title(strcat('orientation clutter at level', num2str(scale)));
    end
    figure;
    imshow(clutter_map, [min_min max_max]); title('collapsed orientation clutter map');
    drawnow;
end

return;