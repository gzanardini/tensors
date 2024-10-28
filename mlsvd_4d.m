function [C,U1,U2,U3,U4]=mlsvd_4d(T)
% mlsvd_4d computes MLSVD for 4th order tensors.

% Computing the factor matrices
[U1,~,~] = svd(mode_n_matricization(T, 1));
[U2,~,~] = svd(mode_n_matricization(T, 2));
[U3,~,~] = svd(mode_n_matricization(T, 3));
[U4,~,~] = svd(mode_n_matricization(T, 4));

% Computing the core tensor
C = mode_n_product(T, U1', 1);
C = mode_n_product(C, U2', 2);
C = mode_n_product(C, U3', 3); 
C = mode_n_product(C, U4', 4); 

end