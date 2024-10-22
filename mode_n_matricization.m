function X = mode_n_matricization(X,n)
% MODE_K_MATRICIZATION takes a tensor X as input and a mode n such
% that n<ndims(X) and returns the mode-n matricization of X.
% INPUTS tensor X, mode n.
% OUPUT mode-n matricization of X.
tsr_dim = size(X);
X = permute(X,[n, 1:n-1, n+1:ndims(X)]);
X = reshape(X,[tsr_dim(n), numel(X)/tsr_dim(n)]);
end