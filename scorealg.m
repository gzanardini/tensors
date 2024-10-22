function rank_SCORE = scorealg(data4d, s)
    % Estimates the rank according to an information theoretic criteria
    % for tensors as described at Yokota, Tatsuya, Namgil Lee, and Andrzej
    % Cichocki. "Robust multilinear tensor rank estimation using higher
    % order singular value decomposition and information criteria."
    % IEEE Transactions on Signal Processing 65.5 (2016): 1196 - 1206.
    %
    % Returns multilinear rank "rank_SCORE" that holds(R_1, R_2, ..., R_N)

    disp('>> Perform Rank Estimation via SCORE Algorithm...');
    
    % rho_hat is suggested to be between 0.0001 and 0.01
    rho_hat = 10^(-5);
    n = ndims(data4d);
    rank_SCORE = zeros(n, 1);  % Preallocate rank_SCORE array
    
    for i = 1:n
        mode_i_H = tens2mat(s, i);  % Convert tensor s to matrix for mode i
        u_j = zeros(1, size(mode_i_H, 2));
        
        for h = 1:size(mode_i_H, 2)
            u_j(h) = mode_i_H(:, h)' * mode_i_H(:, h);  % Compute u_j
        end
        
        [~, high_u_j_idx] = sort(u_j, 'descend');  % Sort u_j in descending order
        
        % Calculate the threshold
        threshold = ceil(rho_hat * prod(size(data4d)) / size(data4d, i));
        
        % Truncate the mode matrix based on threshold
        mode_i_H_P = mode_i_H(:, high_u_j_idx(1:threshold));
        tmp = (1 / threshold) * (mode_i_H_P * mode_i_H_P');
        lambda_ = diag(tmp);
        lambda_ = sort(lambda_, 'descend');  % Sort eigenvalues in descending order
        
        min_rank = size(data4d, i);  % Initialize minimum rank
        min_rank_val = 10^6;  % Initialize large minimum rank value
        
        for r = 1:(size(data4d, i) - 1)
            denominator = 0;
            nominator = 1;
            
            for m = r:size(data4d, i)
                nominator = nominator * lambda_(m)^(1 / (size(data4d, i) - r));
                denominator = denominator + (1 / (size(data4d, i) - r)) * lambda_(m);
            end
            
            min_rank_val_new = -2 * log((nominator / denominator)^(threshold * (size(data4d, i) - r))) + ...
                r * (2 * size(data4d, i) - r) * log(threshold);
            
            if min_rank_val_new < min_rank_val
                min_rank = r;
                min_rank_val = min_rank_val_new;
            end
        end
        
        rank_SCORE(i) = min_rank;
    end    
    disp(['The estimated rank is [', num2str(rank_SCORE'), ']']);
end