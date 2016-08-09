function [out] = RROrientationOppEnergy(L_pyr, numlevels)
% OPP_ENERGY    This runs the oriented opponent energy calculation that
% 		serves as the first stages in Bergen & Landy's (1990)
% 		texture segmentor, except it uses DOOG filters (which actually
% 		don't work as well, but at least we can more easily control the
% 		scale).

hvdd = cell(numlevels, 1);
hv = cell(numlevels, 1);
dd = cell(numlevels, 1);
out = cell(numlevels, 1);
total = cell(numlevels, 1);
noise = 1.0;    % Was 1.5
filterScale = 16/14*1.75;
poolScale = 1.75;
% These probably seem like arbitrary numbers, but it's just trying to get
% three very different feature extraction methods to operate at basically
% the same scales.

for scale=1:numlevels
    % Check this is the right order for Landy/Bergen. RRR
    hvdd{scale} = orient_filtnew(L_pyr{scale},filterScale); %filt with 4 oriented filters 0, 45, 90, 135.  Was sigma = 16/14, orient_filtnew, then 16/14*1.75 to match contrast and other scales.
    % Eventually make this sigma a variable that's passed to this routine.
    %hvdd{scale} is the 4 output images concatenated together, 
    %in the order horizontal, vertical, up-left, and down-right.
    hvdd{scale} = hvdd{scale}.^2;    %local energy
    hvdd{scale} = poolnew(hvdd{scale}, poolScale); %Pools with a gaussian filter.  Was effectively sigma=1, then 1.75 to match 1.75 above.
    % RRR Should look at these results and see if this is the right amount of
    % pooling for the new filters.  It was right for the Landy-Bergen
    % filters.
    hv{scale} = HV(hvdd{scale}); %get the difference image between horizontal and vertical: H-V (0-90)
    dd{scale} = DD(hvdd{scale}); %get the difference image between right and left: R-L (45-135)
    % Normalize by the total response at this scale, assuming the total
    % response is high enough.  If it's too low, we'll never see this
    % orientation.  I'm not sure what to do here -- set it to zeros and
    % it's like that's the orientation.  Maybe output the total response
    % and decide what to do later.  RRR
    total{scale} = sumorients(hvdd{scale})+noise; %add noise based upon sumorients at visibility threshold
    hv{scale} = hv{scale}./total{scale};%normalize the hv and dd image
    dd{scale} = dd{scale}./total{scale};
    out{scale} = [hv{scale} dd{scale}];%out is the 2 output images concatenated 
                                       %together, in the order of hv, dd
end