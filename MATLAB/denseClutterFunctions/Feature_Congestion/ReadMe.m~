% ReadMe.m
%
% Implements two measures of visual clutter (Feature Congestion and Subband
% Entropy) proposed in:
% Ruth Rosenholtz, Yuanzhen Li, and Lisa Nakano. "Measuring Visual Clutter". 
% Journal of Vision, 7(2), 2007. http://www.journalofvision.com/7/2/
%
% Ruth Rosenholtz, Yuanzhen Li, and Lisa Nakano, May 2007.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contents:
%   
% getClutter_FC: computes Feature Congestion clutter, outputs both a scalar 
%    (clutter of the whole image) and a map (local clutter).
% getClutter_SE: computes Subband Entropy clutter, outputs only a scalar.
% colorClutter: computes clutter maps indicating local variability in color
% contrastClutter: computes clutter maps indicating local variability in contrast
% orientationClutter: computes clutter maps indicating local variability in orientation
% 
% (Please see individual routines for more info about parameters and outputs.)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples:
%
% get Feature Congestion clutter of a test map:
%[clutter_scalar_fc, clutter_map_fc] = getClutter_FC('test.jpg');
[clutter_scalar_fc, clutter_map_fc] = getClutter_FC('../Z1C1/1.jpg');

% display the clutter map and output the scalar 
figure, imshow((clutter_map_fc-min(clutter_map_fc(:)))/(max(clutter_map_fc(:))-min(clutter_map_fc(:))));
title('Feature Congestion clutter map');
clutter_scalar_fc
%
% get Subband Entropy clutter of the test map
clutter_se = getClutter_SE('../Z1C1/1.jpg')
%
% compute and display color clutter map(s)
[clutter_levels, clutter_map] = colorClutter('../Z1C1/1.jpg', 3);
% compute and display contrast clutter map(s)
[clutter_levels, clutter_map] = contrastClutter('../Z1C1/1.jpg', 3, 1);
% compute and display orientation clutter map(s)
[clutter_levels, clutter_map] = orientationClutter('../Z1C1/1.jpg', 3);

