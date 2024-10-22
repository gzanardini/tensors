function Z = mode_n_product(X,Y,n)
    % MODE_N_PRODUCT takes tensor X and compatible matrix Y and performs mode-n product between X and Y.
    % INPUT tensor X, matrix Y.
    % OUTPUT tensor Z.
    tsr_dim = size(X);
    X = mode_n_matricization(X,n);
    Z = Y * X;
    Z = reshape(Z,[size(Y, 1),tsr_dim([1:n-1, n+1:end])]);
    Z = permute(Z, [2:n, 1, n+1:numel(tsr_dim)]);
end