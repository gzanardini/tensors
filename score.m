function ranks = score(X, rho)
ranks = zeros(ndims(X), 1);
sz = size(X);
V = cell(ndims(X), 1);

% Initialize eig vectors

for n = 1:ndims(X)
    X_n = mode_n_matricization(X, n);
    [eig_vec, D] = eig((sz(n)/prod(sz)).*X_n*X_n');
    [~, ind] = sort(diag(D));
    V{n} = eig_vec(:,ind);
end

H = mode_n_product(X, V{1}', 1);
for n = 2:ndims(X)
H = mode_n_product(H, V{n}', n);
end

for i = 1:ndims(X)
    H_i = mode_n_matricization(H, i);
    u_j = zeros(1, size(H_i, 2));
    for h = 1:size(H_i, 2)
        u_j(h) = H_i(:, h)' * H_i(:, h); 
    end
    [~, u_j_sorted] = sort(u_j, 'descend');  % Sort u_j in descending order

    threshold = ceil(rho * prod(sz) / size(X, i));
    P_G = H_i(:, u_j_sorted(1:threshold));

    tmp = (1 / threshold) * (P_G * P_G');

    lambda = diag(tmp);
    lambda = sort(lambda, 'descend');

    min_rank = size(X, i);  % Initialize minimum rank
    min_rank_val = 10^6;  % Initialize large minimum rank value

    for r = 1:(size(X, i) - 1)
        denominator = 0;
        nominator = 1;
        
        for m = r:size(X, i)
            nominator = nominator * lambda(m)^(1 / (size(X, i) - r));
            denominator = denominator + (1 / (size(X, i) - r)) * lambda(m);
        end
        
        min_rank_val_new = -2 * log((nominator / denominator)^(threshold * (size(X, i) - r))) + r * (2 * size(X, i) - r) * log(threshold);
        
        if min_rank_val_new < min_rank_val
            min_rank = r;
            min_rank_val = min_rank_val_new;
        end
    end
    
    ranks(i) = min_rank;
end

end