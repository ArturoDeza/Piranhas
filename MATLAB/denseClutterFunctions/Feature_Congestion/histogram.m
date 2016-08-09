function result=histogram(x, nbins)

% to calculate the histogram of signal "x", given the number of bins
% "nbins"; used uniform binning

edges = min(x(:)):(max(x(:))-min(x(:)))/(nbins-1):max(x(:));
result = histc(x(:), edges);