close all;
clear variables;
clc;
format long

%% Generating imaging domain

t_samples = (0:4:120)+eps;

TIC_1 = zeros(10,10,10, length(t_samples));
k_1 = 0.5.*ones(10,10,10);   % local diffusion-related parameter (0.35<)
mu_1 = 30.*ones(10,10,10);         % 2.5/k_1 ;    % mean transit time (lambda/k)
alpha_1 = 1000;        % scale parameter ()

for idx = 1:length(t_samples) %TIC1
    t = t_samples(idx);  
    TIC_1(:,:,:, idx) = alpha_1 .* sqrt(k_1 ./ (2 * pi * t)) .*exp(-((k_1 .* (t - mu_1).^2) ./ (2 * t)));
end

TIC_2 = zeros(4,4,4, length(t_samples));
TIC_3 = zeros(4,4,4, length(t_samples));

k_2 = 3.*ones(4,4,4); 
mu_2 = 25.*ones(4,4,4);           % 4.5/k_2 ;  
alpha_2 = 1600;   

k_3 = 2.*ones(4,4,4);  
mu_3 = 15.*ones(4,4,4);           % 6/k_3;    
alpha_3 = 1200;  

for idx = 1:length(t_samples)
    t = t_samples(idx);  
    TIC_2(:,:,:, idx) = alpha_2 .* sqrt(k_2./(2*pi*t)).*exp(-((k_2.*(t - mu_2).^2)./(2 * t))); %TC_2
    TIC_3(:,:,:, idx) = alpha_3 .* sqrt(k_3./(2*pi*t)).*exp(-((k_3.*(t - mu_3).^2)./(2 * t))); %TC_3
end

G = TIC_1; 
G(2:5,2:5,2:5,:) = TIC_2;
G(2:5,6:9,6:9,:) = TIC_3;

%%% visual
G_test=zeros(10,10,10);
G_test(2:5,2:5,2:5,:) = 2;
G_test(2:5,6:9,6:9,:) = 5;
vol_viz(G_test,'Imaged volume')
colorbar();
saveas(gcf,'figures/regions.eps','epsc');

%% Signal Domain visualization
G_test=G(:,:,:,4);      % for noisy use Y
vol_viz(G_test, 'Original signal at 16sec')
colorbar();

%%  TIC curves

figure;
hold('on')
title('TIC Curves in each region')
plot(squeeze(G(1,1,1,:)),DisplayName='TIC 1')
plot(squeeze(G(3,3,3,:)),DisplayName='TIC 2')
plot(squeeze(G(3,7,7,:)),DisplayName='TIC 3')
legend()
hold('off')
saveas(gcf,'figures/TICs.eps','epsc')

%% Noise

noise_param = 10^.18;

%%% 10^4.2  SNR -60
%%% 10^3.18 SNR -40
%%% 10^2.2  SNR -20
%%% 10^1.18 SNR 0
%%% 10^.18  SNR +20

N = raylrnd(noise_param, [10,10,10,length(t_samples)]);
SNR = snr(G, N);
Y = (G).*N; 
%%
figure;
histogram(Y,'Normalization','pdf')
title('Y - pdf')
xlabel('bins')
ylabel('frequency')

% Taking logarithm
log_Y=log(Y+ eps);

figure;
hist=histogram(log_Y,'Normalization','pdf');
[Ns, Edges]=histcounts(log_Y,'Normalization','pdf');
title('log(Y) - pdf');

%% SVD -- Eqn 3

Y_4 = mode_n_matricization(log_Y, 4);
[U_4, S_4, V_4] = svd(Y_4, 'econ');

figure;
semilogy(diag(S_4), 'o-')
ylabel('\sigma')
title(['SVD of Y_{(4)} - SNR: ' num2str(SNR)])

svd_rank = 7;
U_4t = U_4(:,1:svd_rank);
V_4t = V_4(:,1:svd_rank);
S_4t = S_4(1:svd_rank,1:svd_rank);

rec_Y = U_4t* S_4t * V_4t';
rec_Y = reshape(rec_Y, size(G));
G_hat = exp(rec_Y);

MSE=norm(G_hat-G,'fro')/numel(G);
disp(MSE)

%% Score Alg
rhos = logspace(-5, -3, 50);
n_trials = length(rhos);
ranks=cell(n_trials,1);

for i=1:n_trials
    ranks{i} = score(log_Y, rhos(i))';
    %display(num2str(ranks{i}));

end

string_ranks = cellfun(@(v) mat2str(v), reshape(ranks, [], 1), 'UniformOutput', false);

figure;
histogram(categorical(string_ranks)); 
xlabel('Unique Vectors');
ylabel('Frequency');
title('Histogram of Unique Vectors - \rho sweep from 10^{-4} to 10^{-3}');
xtickangle(45); 

%% For different noise realizations
SNRs=[-60 -50 -40 -30 -20 -10 0 10 20];
noise_factors=[];

for k=1:length(SNRs)
    fac=sqrt(var_G/(0.5*(4-pi)*10^(0.1*SNRs(k))));
    noise_factors=[noise_factors fac];
    
end

rho=0.01;
n_trials = 100;

ranks=cell(length(noise_factors),n_trials);

for i = 1:length(noise_factors)
    for k=1:n_trials
        N = raylrnd(noise_factors(i), [10,10,10,length(t_samples)]);
        Y = (G).*N; 
        log_Y=log(Y+ eps);
        ranks{i,n} = score(log_Y, rho)';
    end
end

string_ranks = cellfun(@(v) mat2str(v), reshape(ranks, [], 1), 'UniformOutput', false);
figure;
histogram(categorical(string_ranks)); 
xlabel('Unique Vectors');
ylabel('Frequency');
title('Histogram of Unique Vectors - \rho sweep from 10^{-4} to 10^{-3}');
xtickangle(45);


%% Reconstruction

[C,U1,U2,U3,U4]=mlsvd_4d(log_Y);
rx = 2; ry = 3; rz = 3; rt = 3;
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

MSE=norm(G_hat-G,'fro')/numel(G);
disp(MSE)

vol_viz(G_hat(:,:,:,6),'ReconAAAAA')
%% Variation of reconstruction
% 
% rec_logY=mode_n_product(Ct,U1t,1);
% rec_logY=mode_n_product(rec_logY,U2t,2);
% rec_logY=mode_n_product(rec_logY,U3t,3);
% rec_logY=mode_n_product(rec_logY,U4t,4);
% 
% G_hat=exp(rec_logY);
% 
% MSE=norm(G_hat-G,'fro')/norm(G,'fro');
% disp(MSE)
% 
% %errors are scaling with alpha

%%

vol_viz(G_hat(:,:,:,6), 'Reconstructed signal')
%% calculate scale factor for rayleygh noise given SNR

var_G=var(G(:));

%%% Heuristic stuff
%%% 10^4.2  SNR -60
%%% 10^3.18 SNR -40
%%% 10^2.2  SNR -20
%%% 10^1.18 SNR 0
%%% 10^.18  SNR +20
%%% noise_factors=[10^4.2 10^3.18 10^2.2 10^1.18 10^.18];

SNRs=[-60 -50 -40 -30 -20 -10 0 10 20];
noise_factors=[];

for k=1:length(SNRs)
    fac=sqrt(var_G/(0.5*(4-pi)*10^(0.1*SNRs(k))));
    noise_factors=[noise_factors fac];
    
end

%%%% Results for report


rho_star=0.01;
num_trials=10;

numel_G=numel(G);

MSEs=zeros(length(SNRs),num_trials);
ranks=cell(length(SNRs),num_trials);
rell_errs=zeros(length(SNRs),num_trials);

norm_G=norm(G,'fro');

for k=1:length(noise_factors) %% iterate SNR from -60 to +20
    for n = 1:num_trials    
            %%% sample the noise
            N = raylrnd(noise_factors(k), [10,10,10,length(t_samples)]);
            Y = (G).*N; 
            log_Y=log(Y+eps);
    
            ranks{k,n} = score(log_Y, rho_star)'; %% rank estimate

            rx=ranks{k,n}(1);
            ry=ranks{k,n}(2);
            rz=ranks{k,n}(3);
            rt=ranks{k,n}(4);            
            
            [C,U1,U2,U3,U4]=mlsvd_4d(log_Y); %% MLSVD of logY
    
            %%% truncation
            U1t = U1(:,1:rx); 
            U2t = U2(:,1:ry); 
            U3t = U3(:,1:rz); 
            U4t = U4(:,1:rt);
            Ct = C(1:rx,1:ry,1:rz,1:rt);

            %%% reconstruction
            rec_logY=mode_n_product(log_Y,(U1t*U1t'),1);
            rec_logY=mode_n_product(rec_logY,(U2t*U2t'),2);
            rec_logY=mode_n_product(rec_logY,(U3t*U3t'),3);
            rec_logY=mode_n_product(rec_logY,(U4t*U4t'),4);
            
            G_hat=exp(rec_logY);

            rell_errs(k,n)=norm(G_hat-G,'fro')/norm_G;
            MSEs(k,n)=norm(G_hat-G,'fro')/numel_G;

    end
end

avgMSE=mean(MSEs,2);
avgRel_errs=mean(rell_errs,2);

%%
figure;
hold off;
plot(SNRs,avgMSE+eps,'Marker','O');
set(gca, 'YScale', 'log') 
xlabel('SNR')
ylabel('avg MSE')
title('avgMSE vs SNR - 1000 runs')
saveas(gcf,'figures/mse vs snr.eps','epsc')

figure;
plot(SNRs,avgRel_errs, 'Marker', 'O');
set(gca, 'YScale','log')
xlabel('SNR')
ylabel('Relative Reconstruction error');
title('avg Relative error vs SNR - 1000 runs')
saveas(gcf,'figures/relerrs.eps','epsc')

