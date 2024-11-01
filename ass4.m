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
% G_test=zeros(10,10,10);
% G_test(2:5,2:5,2:5,:) = 2;
% G_test(2:5,6:9,6:9,:) = 5;
% vol_viz(G_test,'Y-Z slices along the X dimension')
% colorbar();
% saveas(gcf,'figures/regions.eps','epsc');

%% Signal Domain visualization
G_test=G(:,:,:,6);      % for noisy use Y
vol_viz(G_test, 'Original signal at 20sec')
colorbar();

%%  TIC curves
figure;
hold('on')
title('TIC in each region of the clean signal')
plot(t_samples, squeeze(G(1,1,1,:)),DisplayName='TIC 1')
plot(t_samples, squeeze(G(3,3,3,:)),DisplayName='TIC 2')
plot(t_samples, squeeze(G(3,7,7,:)),DisplayName='TIC 3')
legend()
hold('off')
%saveas(gcf,'figures/TICs.eps','epsc')

%% Generating Noise parameters

var_G=var(G(:));
SNRs=[-60 -50 -40 -30 -20 -10 0 10 20];
noise_factors=[];
for k=1:length(SNRs)
    fac=sqrt(var_G/(0.5*(4-pi)*10^(0.1*SNRs(k))));
    noise_factors=[noise_factors fac];
    
end

%% Generating Noisy signal

N = raylrnd(noise_factors(3), [10,10,10,length(t_samples)]);
SNR = snr(G, N);
Y = (G).*N; 

G_test=Y(:,:,:,6);   
vol_viz(G_test, 'Noisy signal at 20sec')
colorbar();

figure;
hold('on')
title(['TIC in each region of the noisy signal - SNR: ' num2str(SNRs(3)) 'dB'])
plot(t_samples, squeeze(Y(1,1,1,:)),DisplayName='TIC 1')
plot(t_samples, squeeze(Y(3,3,3,:)),DisplayName='TIC 2')
plot(t_samples, squeeze(Y(3,7,7,:)),DisplayName='TIC 3')
legend()
hold('off')
%% Log Y + Noise Histograms 
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
title(['SVD of Y_{(4)} - SNR: ' num2str(SNRs(3))])

svd_rank = 3;
U_4t = U_4(:,1:svd_rank);
V_4t = V_4(:,1:svd_rank);
S_4t = S_4(1:svd_rank,1:svd_rank);

rec_Y = U_4t* S_4t * V_4t';
rec_Y = reshape(rec_Y, size(G));
G_hat_SVD = exp(rec_Y);

MSE=norm(G_hat_SVD-G,'fro')/numel(G);
disp(MSE)

figure;
hold('on')
title(['TIC in each region of the reconstructed signal using SVD - SNR: ' num2str(SNRs(3))])
plot(t_samples, squeeze(G_hat_SVD(1,1,1,:)),DisplayName='TIC 1')
plot(t_samples, squeeze(G_hat_SVD(3,3,3,:)),DisplayName='TIC 2')
plot(t_samples, squeeze(G_hat_SVD(3,7,7,:)),DisplayName='TIC 3')
legend()
hold('off')

%% SVD For different SNRs
n_trials=1000;
svd_rank = 3;
MSE_SVD = zeros(length(SNRs), n_trials);

for k=1:length(noise_factors) %% iterate SNR from -60 to +20
    for n = 1:n_trials    
            %%% sample the noise
            N = raylrnd(noise_factors(k), [10,10,10,length(t_samples)]);
            Y = (G).*N; 
            log_Y=log(Y+eps);
            Y_4 = mode_n_matricization(log_Y, 4);
            [U_4, S_4, V_4] = svd(Y_4, 'econ');
            U_4t = U_4(:,1:svd_rank);
            V_4t = V_4(:,1:svd_rank);
            S_4t = S_4(1:svd_rank,1:svd_rank);
            
            rec_Y = U_4t* S_4t * V_4t';
            rec_Y = reshape(rec_Y, size(G));
            G_hat = exp(rec_Y);
            
            MSE_SVD(k,n)=norm(G_hat-G,'fro')/numel(G);
            
    end
end
avgMSE_SVD=mean(MSE_SVD,2);

%% MLSVD -- Eqn 4
% Singular values of each mode
X = log_Y;      % X = G for singular vals of clean signal 
[~, S_1, ~] = svd(mode_n_matricization(X, 1), 'econ');
[~, S_2, ~] = svd(mode_n_matricization(X, 2), 'econ');
[~, S_3, ~] = svd(mode_n_matricization(X, 3), 'econ');
[~, S_4, ~] = svd(mode_n_matricization(X, 4), 'econ');

hold on;
subplot(411)
plot(diag(S_1), 'o-', DisplayName='Mode 1')
ylabel('\sigma')
title(['\Sigma_{(1)} - SNR: ' num2str(SNRs(3))])
set(gca, 'YScale','log')

subplot(412)
plot(diag(S_2), 'o-', DisplayName='Mode 2')
ylabel('\sigma')
title(['\Sigma_{(2)} - SNR: ' num2str(SNRs(3))])
set(gca, 'YScale','log')

subplot(413)
plot(diag(S_3), 'o-', DisplayName='Mode 3')
ylabel('\sigma')
title(['\Sigma_{(3)} - SNR: ' num2str(SNRs(3))])
set(gca, 'YScale','log')

subplot(414)
plot(diag(S_4), 'o-', DisplayName='Mode 4')
set(gca, 'YScale','log')
xlim([1 31])
ylabel('\sigma')
title(['\Sigma_{(4)} - SNR: ' num2str(SNRs(3))])
hold off;
%% Score Alg - To find opt rho
rhos = logspace(-4, -1, 100);
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
title('Histogram of Unique Vectors - \rho sweep from 10^{-4} to 10^{-1}');
xtickangle(45); 

%% For different noise realizations at rho = 0.01

rho=0.01;
n_trials = 1000;
ranks=cell(length(noise_factors),n_trials);

for i = 1:length(noise_factors)
    for k=1:n_trials
        N = raylrnd(noise_factors(i), [10,10,10,length(t_samples)]);
        Y = (G).*N; 
        log_Y=log(Y+ eps);
        ranks{i,k} = score(log_Y, rho)';
    end
end

string_ranks = cellfun(@(v) mat2str(v), reshape(ranks, [], 1), 'UniformOutput', false);
figure;
histogram(categorical(string_ranks)); 
xlabel('Unique Vectors');
ylabel('Frequency');
title('Estimated Ranks \rho = 0.01');
xtickangle(45);

%% Reconstruction using MLSVD

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

G_hat_MLSVD=exp(rec_logY);
%G_hat_MLSVD = G_hat_MLSVD./norm(log_Y, 'fro');         %Uncomment to rescale
MSE=norm(G_hat_MLSVD-G,'fro')/numel(G);
disp(MSE)

vol_viz(G_hat(:,:,:,6),'Reonstructed signal at 16 sec')

figure;
hold('on')
title(['TIC in each region of the reconstructed signal - SNR: ' num2str(SNRs(3))])
plot(t_samples, squeeze(G_hat_MLSVD(1,1,1,:)),DisplayName='TIC 1')
plot(t_samples, squeeze(G_hat_MLSVD(3,3,3,:)),DisplayName='TIC 2')
plot(t_samples, squeeze(G_hat_MLSVD(3,7,7,:)),DisplayName='TIC 3')
legend()
hold('off')

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

%% Calculate scale factor for rayleigh noise given SNR

%%% Heuristic stuff
%%% 10^4.2  SNR -60
%%% 10^3.18 SNR -40
%%% 10^2.2  SNR -20
%%% 10^1.18 SNR 0
%%% 10^.18  SNR +20
%%% noise_factors=[10^4.2 10^3.18 10^2.2 10^1.18 10^.18];


%% Results for report -- MLSVD 

rho_star=0.01;
num_trials=1000;
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
            %G_hat = G_hat./norm(log_Y, 'fro');         %Uncomment to rescale
            rell_errs(k,n)=norm(G_hat-G,'fro')/norm_G;
            MSEs(k,n)=norm(G_hat-G,'fro')/numel_G;

    end
end

avgMSE=mean(MSEs,2);
avgRel_errs=mean(rell_errs,2);

figure;
hold('on')
title(['TIC in each region of the reconstructed signal - SNR: ' num2str(SNRs(3))])
plot(t_samples, squeeze(G_hat(1,1,1,:)),DisplayName='TIC 1')
plot(t_samples, squeeze(G_hat(3,3,3,:)),DisplayName='TIC 2')
plot(t_samples, squeeze(G_hat(3,7,7,:)),DisplayName='TIC 3')
legend()
hold('off')

%% Plotting the AvgMSE vs SNR

figure;
hold on;
plot(SNRs,avgMSE+eps,'Marker','O',DisplayName='MLSVD + SCORE');
set(gca, 'YScale', 'log') 
plot(SNRs,avgMSE_SVD,'Marker','O',DisplayName='SVD');
set(gca, 'YScale', 'log') 
legend()
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

%% Correlation

G_4 = mode_n_matricization(G, 4);
G_SVD_4 = mode_n_matricization(G_hat_SVD, 4);
G_MLSVD_4 = mode_n_matricization(G_hat_MLSVD, 4);
corr1 = round(corr(G(:), G_hat_SVD(:)), 3);
corr2 = round(corr(G(:), G_hat_MLSVD(:)), 3);

%% Visualizing the rank 2 nature of mode 1 unfolding

G_test=zeros(10,10,10);
G_test(2:5,2:5,2:5,:) = 2;
G_test(2:5,6:9,6:9,:) = 5;
imagesc(mode_n_matricization(G_test,1))
title('Mode 1 unfolded matrix')
