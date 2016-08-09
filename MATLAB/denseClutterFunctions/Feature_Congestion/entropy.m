function estimate=entropy(x, nbins)

% to calculate the entropy of signal "x", given the number of bins "nbins"
% used uniform binning in the calculation

nsamples = length(x(:)); 

if ~exist('nbins')
    nbins = ceil(sqrt(nsamples));
end
h=histogram(x,nbins);

sum_h = sum(h);
if sum_h == 0
    sum_h = 1;
end
h = h/sum_h;

h_ = h;
h_(find(h==0)) = 0.01;
logf = log(h_);
estimate = sum(-logf.*h);
