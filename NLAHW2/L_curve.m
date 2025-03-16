%-----------Problem------%
%I asked Chatgpt to create me an ill-conditionned problem in order to have
%a big curvature%
n = 100; 
[U, ~, V] = svd(randn(n, n)); 
sing_values = logspace(0, -10, n); 
A = U * diag(sing_values) * V'; 
x_true = randn(n, 1); 
b = A * x_true; 
noise_level = 0.2; 
noise = noise_level * norm(b) * randn(n, 1) / norm(randn(n, 1)); 
b_noisy = b + noise; 
%end Chatgpt%

%-----Norms-----%
[U, S, V] = svd(A); 
lambda_values = logspace(-10, 1, 50); 
num_lambda = length(lambda_values);
norm_Ax_minus_b = zeros(num_lambda, 1);
norm_Lx = zeros(num_lambda, 1);
for i = 1:num_lambda
    lambda = lambda_values(i);
    x_lambda = V * diag(diag(S) ./ (diag(S).^2 + lambda^2)) * U' * b_noisy; % Regularized solution
    norm_Ax_minus_b(i) = norm(A * x_lambda - b_noisy); % Residual norm
    norm_Lx(i) = norm(lambda * x_lambda); % Solution norm
end
% ------Lcurve-----%
loglog(norm_Lx, norm_Ax_minus_b, 'b-o');
xlabel('log(||Lx_{\lambda}||_2)');
ylabel('log(||Ax_{\lambda} - b||_2)');
title('L-curve for Tikhonov Regularization');

%-----Maximum curvature---%

[~, corner_index] = min(abs(log(norm_Ax_minus_b) - log(norm_Lx))); 
%maybe not the best curve, if I have time do dy/dx
optimal_lambda = lambda_values(corner_index+1);
hold on;
loglog(norm_Lx(corner_index), norm_Ax_minus_b(corner_index), 'ro', 'MarkerFaceColor', 'r');
hold off;
%----Result-----%
disp(['Optimal lambda: ', num2str(optimal_lambda)]);
