function [clutter_levels, clutter_map] = contrastClutter(inputImage, numlevels, filt_sigma, pool_sigma, pix)
% [clutter_levels, clutter_map] = contrastClutter(inputImage, numlevels, filt_sigma, [pool_sigma], [pix])
% 
% Computes the contrast clutter map(s) of an image. Returns:
% "clutter_levels", a cell structure, containing the contrast clutter at a number 
%    of scales specified by numlevels;  -- cell(numlevels,1) --, the n'th level of 
%    which can be accessed using clutter_levels{n}{1}
% "clutter_map", an array of the same size as inputImage, is a single clutter map 
%    collapsed from clutter_levels, which is the clutter measure at multiple scales
%    now the "collapsing" is done by taking the maximal values across scales
% 
% Inputs:
% "inputImage" gives the input. It can be one of the following 3 things: 1. an RGB 
%    image; 2. a string, i.e., file name of an RGB image; 3, a cell structure 
%    containing Gaussian pyramids of the luminance(L) and the chrominance(a,b) 
%    channels. Refer to RRgaussianPyramid.m and computeClutter.m for how to build 
%    and use Gaussian pyramids.
% "numlevels" is the number of levels.
% "filt_sigma" is the sigma of the center-surround DoG1 filter used for
%     computing the contrast 
%   contrast clutter is then computed as the sqrt(variance) of contrast over a local
%     window. The variance can be computed efficiently through linear filtering. 
%     More specifically, var(X) = E(X.^2)-E(X).^2, where E (expectation value) 
%     can be approximated by filtering with a Gaussian window. 
% "pool_sigma" (optional) is the sigma of this Gaussian window. Default =
%    3*filt_sigma. 
% "pix" is for the option of displaying clutter maps. If it's 1, displays the clutter maps. 
%    If it's 0, does not display (useful for batch processing of many images). Defaults to 1.
%
% Reference: 
% Ruth Rosenholtz, Yuanzhen Li, and Lisa Nakano. "Measuring Visual Clutter". 
% Journal of Vision, 7(2), 2007. http://www.journalofvision.com/7/2/

% Ruth Rosenholtz and Yuanzhen Li, Sept 2004.

if ~exist('pool_sigma')
    pool_sigma = 3*filt_sigma;
end
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
    % If im is an MxNx3 matrix, assume it's in RGB.  If it's MxN, assume it's
    % already a grayscale image.
    [m, n, d] = size(im);
    if d==3,
        % if it's RGB, we first convert it into the
        % perceptually-based CIELab color space
        Lab = RGB2Lab(im);
        % and we take the luminance L channel, on which contrast is computed
        L = Lab(:,:,1);
    else
        L = im;
    end

    % We then process the luminance L at multiple scales by creating a Gaussian
    % pyramid by alternately smoothing and subsampling the image
    L_pyr = RRgaussianPyramid(L, numlevels);
end

% We then compute a form of "contrast-energy" by filtering the luminance
% channel L by a center-surround filter and squaring (or taking the absolute 
% values of) the filter outputs. The center-surround filter is a DoG1 filter 
% with std 'filt_sigma'.
contrast = RRcontrast1channel(L_pyr, filt_sigma);
% contrast returned is 1 a cell structure -- cell(numlevels, 1) --.

% initiate clutter_map and clutter_levels:
[m, n] = size(contrast);
clutter_levels = cell(m,n);   
clutter_map = zeros(size(L_pyr(1,1)));

% Get a Gaussian filter for computing the variance of contrast
% Since we used a Gaussian pyramid to find contrast features, these filters 
% have the same size regardless of the scale of processing.
bigG = RRgaussfilter1D(round(pool_sigma*2), pool_sigma);

for scale=1:m,
    for channel=1:n
        % var(X) = E(X.^2) - E(X).^2
        % get E(X) by filtering X with a 1-D Gaussian window separably in x
        %   and y directions
        meanD = RRoverlapconv(bigG, contrast{scale}{channel});
        meanD = RRoverlapconv(bigG', meanD);
        % get E(X.^2) by filtering X.^2 with a 1-D Gaussian window separably in x
        %   and y directions
        meanD2 = RRoverlapconv(bigG, contrast{scale}{channel}.^2);
        meanD2 = RRoverlapconv(bigG', meanD2);
        
        % get variance by var(X) = E(X.^2) - E(X).^2
        stddevD = sqrt(abs(meanD2 - meanD.^2));
        clutter_levels{scale}{channel} = stddevD; 
    end
end

% collapse over scales by taking the maximum 

% first get a Gaussian kernel to upsample the clutter maps on bigger scales
% so that the clutter maps would have the same sizes, and max can be taken
% across scales
kernel_1d = [0.05 0.25 0.4 0.25 0.05];
kernel_2d = conv2(kernel_1d, kernel_1d');
clutter_map = clutter_levels{1}{1};
for scale=2:m
    clutter_here = clutter_levels{scale}{1};
    for kk = scale:-1:2
        clutter_here = upConv(clutter_here, kernel_2d, 'reflect1', [2 2]);
    end
    common_sz = min(size(clutter_map), size(clutter_here));
    clutter_map(1:common_sz(1), 1:common_sz(2)) = max(clutter_map(1:common_sz(1), 1:common_sz(2)), clutter_here(1:common_sz(1), 1:common_sz(2)));
end



% to display the clutter maps:
if (pix==1),
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
        subplot(ny, nx, scale); imshow(clutter_levels{scale}{1}, [min_min max_max]); title(strcat('contrast clutter at level', num2str(scale)));
    end
    figure;
    imshow(clutter_map, [min_min max_max]); title('collapsed contrast clutter map');
    drawnow;
end

return;