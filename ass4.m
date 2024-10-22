%% Generating imaging domain

t_samples = (4:4:120)+eps;

TIC_1 = zeros(10,10,10, length(t_samples));
k_1 = 0.2.*abs(randn(10,10,10));          % local diffusion-related parameter (0.35<)
mu_1 = 2.5/k_1 ;    % mean transit time (lambda/k)
alpha_1 = 1;        % scale parameter ()

for idx = 1:length(t_samples) %TIC1
    t = t_samples(idx);  
    TIC_1(:,:,:, idx) = alpha_1 .* sqrt(k_1 ./ (2 * pi * t)) .*exp(-((k_1 .* (t - mu_1).^2) ./ (2 * t)));
end

TIC_2 = zeros(4,4,4, length(t_samples));
TIC_3 = zeros(4,4,4, length(t_samples));

k_2 = 0.55 + 0.10 .*randn(4,4,4);  
mu_2 = 4.5/k_2 ;  
alpha_2 = 0.52;   

k_3 = 0.8 + 0.10.*randn(4,4,4);  
mu_3 = 6/k_3;    
alpha_3 = 0.6;  

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

plot(squeeze(G(3,3,3,:)),DisplayName='TIC 2')
plot(squeeze(G(3,7,7,:)),DisplayName='TIC 3')
legend()
hold('off')

%% Slicing volume visual

G_test=zeros(10,10,10);
G_test(2:5,2:5,2:5,:) = 2;
G_test(2:5,6:9,6:9,:) = 5;

x=0:1:9;
y=0:1:9;
z=0:1:9;

[x,y,z]=meshgrid(x,y,z);

xslice = [3,3,3];    % location of y-z planes
yslice = 3;          % location of x-z plane
zslice = [2,0];         % location of x-y planes

figure();
slice(x,y,z,G_test,xslice,yslice,zslice)
xlabel('x')
ylabel('y')
zlabel('z')
title('Slices of the imaged volume - No time index')

%% Signal visualization
G_test=G(:,:,:,6);      % for noisy use Y
x=0:1:9;
y=0:1:9;
z=0:1:9;

[x,y,z]=meshgrid(x,y,z);

xslice = [3,3,3];    % location of y-z planes
yslice = 3;          % location of x-z plane
zslice = [2,0];         % location of x-y planes

figure();
slice(x,y,z,G_test,xslice,yslice,zslice)
xlabel('x')
ylabel('y')
zlabel('z')
title('Original signal at 16sec')


%% Noise
noise_param = 10^-4;
N = raylrnd(noise_param, [10,10,10,length(t_samples)]);
Y = (G+eps).*N; 

figure;
histogram(Y*1e6,'Normalization','pdf')
xlabel('bins')
ylabel('frequency')

% Taking logarithm
log_Y=log(Y);

figure;
hist=histogram(log_Y,'Normalization','pdf');
[Ns, Edges]=histcounts(log_Y,'Normalization','probability');

%% SVD -- Eqn 3
Y_4 = mode_n_matricization(log_Y,4);
[U_4, S_4, V_4] = svd(Y_4*Y_4', 'econ');
plot(sqrt(diag(S_4)))

%% MLSVD -- Eqn 4
[C,U1,U2,U3,U4]=mlsvd_4d(Y_4*Y_4');
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
histogram(ranks)


