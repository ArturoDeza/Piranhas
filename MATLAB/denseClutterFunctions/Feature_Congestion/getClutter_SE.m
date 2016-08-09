function [clutter_se] = getClutter_SE(map_name, wlevels, wght_chrom)

% [clutter_se] = getClutter_SE(map_name, [wlevels], [wght_chrom])
% Subband Entropy measure of visual clutter.
% Outputs:
%     "clutter_se": the subband entropy clutter of the image.
% Inputs: 
%     "map_name": input image. It can be a string (file name of the image),
%       or an array (the image itself).
%     "wlevels":  the number of scales (optional, default 3)
%     "wght_chrom": the weight on chrominance (optional, default 0.0625)

% This measure (Subband Entropy) of visual clutter is based on the notion
% that clutter is related to the number of bits required for subband
% (wavelet) image coding.
%
% Reference: 
% Ruth Rosenholtz, Yuanzhen Li, and Lisa Nakano. "Measuring Visual Clutter". 
% Journal of Vision, 7(2), 2007. http://www.journalofvision.com/7/2/


% Ruth Rosenholtz, Yuanzhen Li, and Lisa Nakano, March 2007.

if ~exist('wlevels')
    wlevels = 3; 
end
if ~exist('wght_chrom')
    wght_chrom = 0.0625; 
end

% load input image
if ischar(map_name)
    try
        [tmp, map] = imread(map_name);
    catch
        error(sprintf('Unable to open %s image file.', map_name))
    end
    map = double(tmp);
else
    if isnumeric(map_name)
        map = double(map_name);
    end
end

[ht, wth, dchrom] = size(map);
% if the image is a color rgb image, convert it to luminance and
% chrominance (from RGB to CIE Lab):
if dchrom == 3
    map_lab = RGB2Lab(map);
else
    map_lab = map;
end

wor = 4;
% luminance channel:
map = map_lab(:,:,1);
en_band = band_entropy(map, wlevels, wor);
clutter_se = mean(en_band);

if dchrom == 1
    return;
end

% chrominance channels:
for jj = 2:3
    map = map_lab(:,:,jj);
    if max(map(:))-min(map(:)) < 0.008
        map = zeros(size(map));
    end
    en_band = band_entropy(map, wlevels, wor);
    clutter_se = clutter_se + wght_chrom*mean(en_band);
end

clutter_se = clutter_se/(1+2*wght_chrom);
