function clutter_levels = computeOrientationClutter(L_pyr)
%
% [clutter_levels] = computeOrientationClutter(L_pyr)
% computes the orientation clutter maps. Returns:
% clutter_levels, a cell structure, containing the orientation clutter at a 
%   number of scales specified by numlevels;  -- cell(numlevels,1) --, the n'th 
%   level of which can be accessed using clutter_levels{n}{1}
% 
% input:
% L_pyr
%   the Gaussian pyramid of L (from CIELab color space)
%   the Gaussian pyramid is computed by alternately blurring and
%   subsampling the L channels
% Orientation clutter is computed as the "volume" of an orientation distribution
%   ellipsoid, which is the determinant of covariance matrix. Treats cos(2 theta)
%   and sin(2 theta) (computed from OrientedOppEnergy) as a two-vector, and gets
%   The covariance of this two-vector.  The covariance 
%   matrix can be computed efficiently through linear filtering. More 
%   specifically, cov(X,Y) = E(XY)-E(X)E(Y), where E (expectation value) 
%   can be approximated by filtering with a Gaussian window. 
% poolScale is set to 7/2.
%
% Reference (though there is no orientation clutter in this reference):
% Ruth Rosenholtz, Yuanzhen Li, Jonathan Mansfield, and Zhenlan Jin. 
% Feature Congestion: A Measure of Display Clutter. CHI '05: Proc. of the SIGCHI conference 
% on Human factors in computing systems. May 2005. 761-770.  
%
% Based upon RRcomputeOrientationSaliency
% Ruth Rosenholtz, May 2006
%
% This currently seems far too dependent on luminance contrast.  Check into
% why this is so -- I thought we were normalizing by local contrast.

noise = 0.001;  % Was eps, but that gave too much orientation noise in the saliency maps.  Then changed to 0.000001
poolScale = 7/2;

numlevels = length(L_pyr);
Dc = cell(numlevels, 1);  % mean "cos 2 theta" at distractor scale
Ds = cell(numlevels, 1);  % mean "sin 2 theta" at distractor scale

% Get approximations to cos(2theta) and sin(2theta) from oriented opponent
% energy, at each of the numlevels of the pyramid
[angles] = RROrientationOppEnergy(L_pyr, numlevels);

% Compute the two-vector [meancos, meansin] at each scale, as well as the
% things we need to compute the mean and covariance of this two-vector at
% the larger, distractor scale.

bigG = RRgaussfilter1D(round(8*poolScale), 4*poolScale);  
maxbigG = max(bigG)^2;

covMx = cell(numlevels, 2, 2);  

for i=1:numlevels
    [m,n]=size(angles{i});
    cmx=angles{i}(:,1:n/2);
    smx=angles{i}(:,n/2+1:n);
          
    % Pool to get means at distractor scale. In pooling, don't pool over the target
    % region (implement this by pooling with a big Gaussian, then
    % subtracting the pooling over the target region computed above.  Note,
    % however, that we first need to scale the target region pooling so
    % that its peak is the same height as this much broader Gaussian used
    % to pool over the distractor region.
    Dc{i} = RRoverlapconv(bigG, cmx);
    Dc{i} = RRoverlapconv(bigG', Dc{i});
    Ds{i} = RRoverlapconv(bigG, smx);
    Ds{i} = RRoverlapconv(bigG', Ds{i});

    % Covariance matrix elements.  Compare with computations in
    % RRStatisticalSaliency.  I tried to match computeColorClutter, but I
    % don't remember the meaning of some of the terms I removed.  XXX
    covMx{i}{1}{1} = RRoverlapconv(bigG, cmx.^2);
    covMx{i}{1}{1} = RRoverlapconv(bigG', covMx{i}{1}{1}) - Dc{i}.^2 + noise;
    covMx{i}{1}{2} = RRoverlapconv(bigG, cmx.*smx);
    covMx{i}{1}{2} = RRoverlapconv(bigG', covMx{i}{1}{2}) - Dc{i}.*Ds{i};
    covMx{i}{2}{2} = RRoverlapconv(bigG, smx.^2);
    covMx{i}{2}{2} = RRoverlapconv(bigG', covMx{i}{2}{2}) - Ds{i}.^2 + noise;

    % Get determinant of covariance matrix, which is the volume of the
    % covariance ellipse
    detIm = covMx{i}{1}{1}.*covMx{i}{2}{2} - covMx{i}{1}{2}.^2;
    % Take the square root considering variance is squared, and the square
    % root again, since this is the area and the contrast measure is a "length"
    clutter_levels{i}{1} = detIm.^(1/4);
end

return;