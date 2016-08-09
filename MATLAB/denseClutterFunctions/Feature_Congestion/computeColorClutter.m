function [clutter_levels] = computeColorClutter(L_pyr, a_pyr, b_pyr, sigD)
%
% [clutter_levels] = computeColorClutter(L_pyr, a_pyr, b_pyr, sigD)
% computes the color clutter maps. Returns:
% clutter_levels, a cell structure, containing the color clutter at a 
%   number of scales specified by numlevels;  -- cell(numlevels,1) --, the n'th 
%   level of which can be accessed using clutter_levels{n}{1}
% 
% input:
% L_pyr, a_pyr, b_pyr
%   the Gaussian pyramids of L, a, b (CIELab color space)
%   the Gaussian pyramids are computed by alternately blurring and
%   subsampling the L,a,b channels
% Color clutter is computed as the "volume" of a color distribution
%   ellipsoid, which is the determinant of covariance matrix. Covariance 
%   matrix can be computed efficiently through linear filtering. More 
%   specifically, cov(X,Y) = E(XY)-E(X)E(Y), where E (expectation value) 
%   can be approximated by filtering with a Gaussian window. 
% sigD is the sigma (standard deviation) of this Gaussian window. 
%
% Reference:
% Ruth Rosenholtz, Yuanzhen Li, Jonathan Mansfield, and Zhenlan Jin. 
% Feature Congestion: A Measure of Display Clutter. CHI '05: Proc. of the SIGCHI conference 
% on Human factors in computing systems. May 2005. 761-770.  

% Ruth Rosenholtz and Yuanzhen Li, Sept 2004


% the pyramids have this number of levels:
numlevels = length(L_pyr);
if length(a_pyr)~=numlevels | length(b_pyr)~=numlevels,
    error('L, a, and b channels must have the same number of levels in the Gaussian pyramid\n');
end

% initiatialization
covMx = cell(numlevels, 3, 3);  
clutter_levels = cell(numlevels, 1);
DL = cell(numlevels, 1);
Da = cell(numlevels, 1);
Db = cell(numlevels, 1);

% sensitivitis to the L,a,and b channels are different, therefore we use
% deltaL2, deltaa2, and deltab2 to "scale" the L,a,b axes when computing
% the covariance matrix. Eventually these numbers should be vary according
% to the spatial scales, mimicing our visual system's sensitivity function
deltaL2 = (0.0007)^2;
deltaa2 = 0.1^2;
deltab2 = 0.05^2;

% Get a Gaussian filter for computing the covariance
bigG = RRgaussfilter1D(round(2*sigD), sigD);

for i=1:numlevels,
    % get E(X) by filtering X with a 1-D Gaussian window separably in x
    %   and y directions:
    DL{i} = RRoverlapconv(bigG, L_pyr{i});
    DL{i} = RRoverlapconv(bigG', DL{i});    % E(L)
    Da{i} = RRoverlapconv(bigG, a_pyr{i});
    Da{i} = RRoverlapconv(bigG', Da{i});    % E(a)
    Db{i} = RRoverlapconv(bigG, b_pyr{i});
    Db{i} = RRoverlapconv(bigG', Db{i});    % E(b)
    
    % Covariance matrix 
    % covMx(L,a,b) = | cov(L,L)  cov(L,a)  cov(L,b) |
    %                | cov(a,L)  cov(a,a)  cov(a,b) |
    %                | cov(b,L)  cov(b,a)  cov(b,b) |
    % where cov(X,Y) = E(XY) - E(X)E(Y)
    %   and if X is the same as Y, then it's the variance var(X) =
    %   E(X.^2)-E(X).^2
    % and as cov(X,Y) = cov(Y,X), covMx is symmetric
    % covariance matrix elements:
    covMx{i}{1}{1} = RRoverlapconv(bigG, L_pyr{i}.^2);
    covMx{i}{1}{1} = RRoverlapconv(bigG', covMx{i}{1}{1}) - DL{i}.^2 + deltaL2;  % cov(L,L) + deltaL2
    covMx{i}{1}{2} = RRoverlapconv(bigG, L_pyr{i}.*a_pyr{i});
    covMx{i}{1}{2} = RRoverlapconv(bigG', covMx{i}{1}{2}) - DL{i}.*Da{i};        % cov(L,a)
    covMx{i}{1}{3} = RRoverlapconv(bigG, L_pyr{i}.*b_pyr{i});
    covMx{i}{1}{3} = RRoverlapconv(bigG', covMx{i}{1}{3}) - DL{i}.*Db{i};        % cov(L,b)
    covMx{i}{2}{2} = RRoverlapconv(bigG, a_pyr{i}.^2);
    covMx{i}{2}{2} = RRoverlapconv(bigG', covMx{i}{2}{2}) - Da{i}.^2 + deltaa2;  % cov(a,a) + deltaa2
    covMx{i}{2}{3} = RRoverlapconv(bigG, a_pyr{i}.*b_pyr{i});
    covMx{i}{2}{3} = RRoverlapconv(bigG', covMx{i}{2}{3}) - Da{i}.*Db{i};        % cov(a,b)
    covMx{i}{3}{3} = RRoverlapconv(bigG, b_pyr{i}.^2);    
    covMx{i}{3}{3} = RRoverlapconv(bigG', covMx{i}{3}{3}) - Db{i}.^2 + deltab2;  % cov(b,b) + deltab2

    % Get the determinant of covariance matrix
    % which is the "volume" of the covariance ellipsoid
    detIm = covMx{i}{1}{1}.*(covMx{i}{2}{2}.*covMx{i}{3}{3}-covMx{i}{2}{3}.*covMx{i}{2}{3})...
        - covMx{i}{1}{2}.*(covMx{i}{1}{2}.*covMx{i}{3}{3}-covMx{i}{2}{3}.*covMx{i}{1}{3})...
        + covMx{i}{1}{3}.*(covMx{i}{1}{2}.*covMx{i}{2}{3}-covMx{i}{2}{2}.*covMx{i}{1}{3});
    % take the square root considering variance is squared, and the cube
    % root, since this is the volume and the contrast measure is a "length"
    clutter_levels{i}{1} = sqrt(detIm).^(1/3); 
end

return;