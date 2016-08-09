function [clutter_scalar_fc, clutter_map_fc] = getClutter_FC(filename, p);

% [clutter_scalar_fc, clutter_map_fc] = getClutter_FC(filename, [p]);
% computes Feature Congestion measure of visual clutter. 
% Outputs:
%   "clutter_scalar_fc" is a scalar, which gives the Feature Congestion 
%     clutter of the whole image.
%   "clutter_map_fc" is a clutter map (same size as the input image), 
%     which gives local clutter information.
% Inputs: 
%    "filename": the file name of an image
%    "p" (optional, default 1): a parameter when combining local clutter over 
%      space; the combination can be considered Minkowski distance of order p

% This measure (Feature Congestion) of visual clutter is related to the
% local variability in certain key features, e.g., color, contrast, and
% orientation.
%
% Reference: 
% Ruth Rosenholtz, Yuanzhen Li, and Lisa Nakano. "Measuring Visual Clutter". 
% Journal of Vision, 7(2), 2007. http://www.journalofvision.com/7/2/

% Ruth Rosenholtz, Yuanzhen Li, and Lisa Nakano, March 2007.

if ~exist('p')
    p = 1;
end

% compute local clutter in color, contrast, and orientation
% (please see "computeClutter.m" for info about the paramters and the outputs)
[color_clutter, contrast_clutter, orient_clutter] = computeClutter(filename, 3, 1, 3, 3, 0, 0, 0);

% combine color, contrast, and orientation:
clutter_map_fc = color_clutter{2}/0.2088 + contrast_clutter{2}/0.0660 + orient_clutter{2}/0.0269;

% combine over space using a Minkowski mean of order p, then take the average
clutter_scalar_fc = mean(clutter_map_fc(:).^p).^(1/p);