function [clutter_levels, clutter_map] = colorClutter(inputImage, numlevels, pool_sigma, pix)
% [clutter_levels, clutter_map] = colorClutter(inputImage, numlevels, [pool_sigma], [pix])

% computes the color clutter map(s) of an image. Returns:
% "clutter_levels", a cell structure, containing the color clutter at a number of scales 
%    specified by numlevels;  -- cell(numlevels,1) --, the n'th level of which can be 
%    accessed using clutter_levels{n}{1}
% "clutter_map", an array of the same size as the input image, is a single clutter 
%    map collapsed from clutter_levels, which is the clutter measure at multiple 
%    scales now the "collapsing" is done by taking the maximal values across scales
% 
% Inputs:
% "inputImage" gives the input. It can be one of the following 3 things: 1. an RGB 
%    image; 2. a string, i.e., file name of an RGB image; 3, a cell structure 
%    containing Gaussian pyramids of the luminance(L) and the chrominance(a,b) 
%    channels. Refer to RRgaussianPyramid.m and computeClutter.m for how to build 
%    and use Gaussian pyramids.
% "numlevels" is the number of levels.
% Color clutter is computed as the "volume" of a color distribution
%    ellipsoid, which is the determinant of covariance matrix. Covariance 
%    matrix can be computed efficiently through linear filtering. More 
%    specifically, cov(X,Y) = E(XY)-E(X)E(Y), where E (expectation value) 
%    can be approximated by filtering with a Gaussian window. 
% "pool_sigma" is the sigma (standard deviation) of this Gaussian window.
%    Defaults to 3.
% "pix" is for the option of displaying clutter maps. If it's 1, displays the clutter maps. 
%    If it's 0, does not display (useful for batch processing of many images). Defaults to 1.
%
% Reference: 
% Ruth Rosenholtz, Yuanzhen Li, and Lisa Nakano. "Measuring Visual Clutter". 
% Journal of Vision, 7(2), 2007. http://www.journalofvision.com/7/2/

% Ruth Rosenholtz and Yuanzhen Li, Sept 2004.

if ~exist('pool_sigma')
    pool_sigma = 3;
end
if ~exist('pix')
    pix = 1;
end

if iscell(inputImage)
    L_pyr = inputImage(:,:,1);
    a_pyr = inputImage(:,:,2);
    b_pyr = inputImage(:,:,3);
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

    % Get Gaussian pyramids (one for each of L,a,b)
    L_pyr = RRgaussianPyramid(Lab(:,:,1), numlevels);
    a_pyr = RRgaussianPyramid(Lab(:,:,2), numlevels);
    b_pyr = RRgaussianPyramid(Lab(:,:,3), numlevels);
end

% Compute clutter
[clutter_levels] = computeColorClutter(L_pyr, a_pyr, b_pyr, pool_sigma);

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
        subplot(ny, nx, scale); imshow(clutter_levels{scale}{1}, [min_min max_max]); title(strcat('color clutter at level', num2str(scale)));
    end
    figure;
    imshow(clutter_map, [min_min max_max]); title('collapsed color clutter map');
    drawnow;
end

return;