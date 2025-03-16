% Load your matrix A from 'Test.mat'
load('Test.mat'); 
%whos('-file', 'test.mat'); % A is 10 by 10 matrix

b = rand(10, 1);  
k = 5; 

%------TSVD-----%
[U, S, V] = csvd(A); 
%[x_tsvd, rho_tsvd, eta_tsvd] = tsvd(U, S, V, b, k);  % TSVD solution using regtools
k = min(k, min(size(S)));
singular_values_svd=diag(S);
%formula for TSVD
TSVD_solution = V(:, 1:k) * diag(1 ./ diag(S(1:k, 1:k))) * U(:, 1:k)' * b;

%-------TGSVD-----%
L = eye(size(A));  
[U_gsvd, sm, X, V_gsvd, W] = cgsvd(A, L);
%[x_tgsvd, rho_tgsvd, eta_tgsvd] = tgsvd(U_gsvd, S_gsvd, X_gsvd, b, k);  % TGSVD solution using regtools
singular_values_cgsvd= 1 ./ diag(sm(size(A,1):end, size(A,2):end));

%formula for TGSVD with p=n
n = size(U_gsvd, 1);
%Try to find the same order%
Sigma_gsvd_flipped = diag(sm);
Sigma_gsvd_flipped = [Sigma_gsvd_flipped(1:end-k); flipud(Sigma_gsvd_flipped(end-k+1:end))];
U_gsvd_flipped = U_gsvd;
U_gsvd_flipped(:, end-k+1:end) = fliplr(U_gsvd(:, end-k+1:end));

X_gsvd_flipped = X;
X_gsvd_flipped(:, end-k+1:end) = fliplr(X(:, end-k+1:end));
Sigma_gsvd_inv = zeros(n, n); 
Sigma_gsvd_inv(end-k+1:end, end-k+1:end) = diag(1./diag(Sigma_gsvd_flipped(end-k+1:end, end-k+1:end))); 


TGSVD_solution = X_gsvd_flipped * Sigma_gsvd_inv * U_gsvd_flipped' * b;
%-----Coincide----%;
difference = norm(TSVD_solution - TGSVD_solution);

fprintf('Difference between TSVD and permuted TGSVD solutions: %f\n', difference);
figure;
plot(1:length(TSVD_solution), TSVD_solution, 'b-o', 'DisplayName', 'TSVD Solution'); 
hold on;
plot(1:length(TGSVD_solution), TGSVD_solution, 'r-*', 'DisplayName', 'TGSVD Solution');
hold off;
xlabel('Index');
ylabel('Solution Value');
title('Comparison of TSVD and TGSVD Solutions');
legend('show');
% 1. Singular values
[singular_values_cgsvd_sorted, idx_cgsvd] = sort(singular_values_cgsvd, 'descend');
figure;
semilogy(singular_values_svd, 'b-o');
hold on;
semilogy(singular_values_cgsvd, 'r-x');
legend('SVD', 'GSVD');
title('Comparison of Singular Values');
xlabel('Index');
ylabel('Singular Value');
hold off;

