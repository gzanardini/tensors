%% Generating imaging domain

t_samples = (0:4:120)+eps;

TIC_1 = zeros(10,10,10, length(t_samples));
k_1 = 0.5.*ones(10,10,10); % 0.2.*abs(randn(10,10,10));          % local diffusion-related parameter (0.35<)
mu_1 = 2.5/k_1 ;    % mean transit time (lambda/k)
alpha_1 = 1000;        % scale parameter ()

for idx = 1:length(t_samples) %TIC1
    t = t_samples(idx);  
    TIC_1(:,:,:, idx) = alpha_1 .* sqrt(k_1 ./ (2 * pi * t)) .*exp(-((k_1 .* (t - mu_1).^2) ./ (2 * t)));
end

TIC_2 = zeros(4,4,4, length(t_samples));
TIC_3 = zeros(4,4,4, length(t_samples));

k_2 = 3.*ones(4,4,4); % 0.55 + 0.10 .*randn(4,4,4);  
mu_2 = 4.5/k_2 ;  
alpha_2 = 1600;   

k_3 = 2.*ones(4,4,4); %0.8 + 0.10.*randn(4,4,4);  
mu_3 = 6/k_3;    
alpha_3 = 1200;  

for idx = 1:length(t_samples)
    t = t_samples(idx);  
    TIC_2(:,:,:, idx) = alpha_2 .* sqrt(k_2./(2*pi*t)).*exp(-((k_2.*(t - mu_2).^2)./(2 * t))); %TC_2
    TIC_3(:,:,:, idx) = alpha_3 .* sqrt(k_3./(2*pi*t)).*exp(-((k_3.*(t - mu_3).^2)./(2 * t))); %TC_3
end

G = TIC_1; 
G(2:5,2:5,2:5,:) = TIC_2;
G(2:5,6:9,6:9,:) = TIC_3;

%%  TIC curves

figure;
hold('on')
plot(squeeze(G(1,1,1,:)),DisplayName='TIC 1')
plot(squeeze(G(3,3,3,:)),DisplayName='TIC 2')
plot(squeeze(G(3,7,7,:)),DisplayName='TIC 3')
legend()
hold('off')

%% Slicing volume visual

G_test=zeros(10,10,10);
G_test(2:5,2:5,2:5,:) = 2;
G_test(2:5,6:9,6:9,:) = 5;
vol_viz(G_test,'Slices of the imaged volume - No time index')

%% Signal visualization
G_test=G(:,:,:,6);      % for noisy use Y
vol_viz(G_test, 'Original signal at 16sec')

%% Noise
noise_param = 10^1;
N = raylrnd(noise_param, [10,10,10,length(t_samples)]);
SNR = snr(G, N)
Y = (G).*N; 

figure;
histogram(Y,'Normalization','pdf')
xlabel('bins')
ylabel('frequency')

% Taking logarithm
log_Y=log(Y+ eps);

figure;
hist=histogram(log_Y,'Normalization','pdf');
[Ns, Edges]=histcounts(log_Y,'Normalization','probability');

%% SVD -- Eqn 3
Y_4 = mode_n_matricization(log_Y,4);
[U_4, S_4, V_4] = svd(Y_4, 'econ');
semilogy(diag(S_4))

%% MLSVD -- Eqn 4
[C,U1,U2,U3,U4]=mlsvd_4d(log_Y);
rx = 5; ry = 5; rz = 5; rt = 5;

% Truncation
U1t = U1(:,1:rx); 
U2t = U2(:,1:ry); 
U3t = U3(:,1:rz); 
U4t = U4(:,1:rt);
Ct = C(1:rx,1:ry,1:rz,1:rt);

%% Score Alg
rhos = 0.0001:0.0001:0.01;
n_trials = length(rhos);
ranks=cell(n_trials,1);

for i=1:n_trials
    ranks{i} = score(log_Y, rhos(i))';
    display(num2str(ranks{i}));

end
%%
string_ranks = cellfun(@(v) mat2str(v), ranks, 'UniformOutput', false);

[uniqueVectors, ~, idx] = unique(string_ranks);
counts = histcounts(idx, unique(idx));

figure;
bar(counts);
set(gca, 'XTickLabel', uniqueVectors); 
xlabel('Unique Vectors');
ylabel('Frequency');
title('Histogram of Unique Vectors');
xtickangle(45); 

%% Reconstruction

[C,U1,U2,U3,U4]=mlsvd_4d(log_Y);
rx = 3; ry = 4; rz = 4; rt = 4;
%rx = 10; ry = 10; rz = 10; rt = 30;

% Truncation
U1t = U1(:,1:rx); 
U2t = U2(:,1:ry); 
U3t = U3(:,1:rz); 
U4t = U4(:,1:rt);
Ct = C(1:rx,1:ry,1:rz,1:rt);

rec_logY=mode_n_product(log_Y,(U1t*U1t'),1);
rec_logY=mode_n_product(rec_logY,(U2t*U2t'),2);
rec_logY=mode_n_product(rec_logY,(U3t*U3t'),3);
rec_logY=mode_n_product(rec_logY,(U4t*U4t'),4);

G_hat=exp(rec_logY);

MSE=norm(G_hat-G,'fro')/norm(G,'fro');
disp(MSE)

%% Variation of reconstruction

rec_logY=mode_n_product(Ct,U1t,1);
rec_logY=mode_n_product(rec_logY,U2t,2);
rec_logY=mode_n_product(rec_logY,U3t,3);
rec_logY=mode_n_product(rec_logY,U4t,4);

G_hat=exp(rec_logY);

MSE=norm(G_hat-G,'fro')/norm(G,'fro');
disp(MSE)

%errors are scaling with alpha

%%
vol_viz(rec_logY(:,:,:,6), 'Reconstructed signal')