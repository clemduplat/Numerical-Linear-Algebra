%----------Import matrices---------%
matrix_filename = 'Test.mtx';  
fileID = fopen(matrix_filename, 'r');
data = textscan(fileID, '%f %f %f', 'HeaderLines', 1);
fclose(fileID);

rows = data{1};
columns = data{2};
values = data{3};
n = max([rows; columns]);
A = zeros(n, n);
for i = 1:length(values)
    A(rows(i), columns(i)) = values(i);
end
kmax = 100;  
r = randn(n, 1);  
%r=r/norm(r);
nrm_A = norm(A);

%----------------Testing my functions----------------%
W_j = zeros(kmax,kmax);
W_j_1 = zeros(kmax,kmax);
W_j_2 = zeros(kmax,kmax);
%---Normal Lanczos---%
tic;
[Q_k,T_k,r,err_ind] = Lanczos_HW1(A,kmax,r,nrm_A);
elapsed_time = toc;
fprintf('Elapsed time for first algorithm: %.4f seconds\n', elapsed_time);

%---Lanczos 1---%
tic;
[Q_k1, T_k1, r1, err_ind1]=Lanczos_1(A,kmax,r,nrm_A);
elapsed_time = toc;
fprintf('Elapsed time for second algorithm: %.4f seconds\n', elapsed_time);

%---Lanczos 2---%
tic;
[Q_k2, T_k2, r2, err_ind2]=Lanczos_2(A,kmax,r,nrm_A);
elapsed_time = toc;
fprintf('Elapsed time for third algorithm: %.4f seconds\n', elapsed_time);

%-------Computing-----%
ritz_values = eig(T_k);
ritz_values1 = eig(T_k1);
ritz_values2 = eig(T_k2);
orthogonality_loss=zeros(kmax, 1);
orthogonality_loss1=zeros(kmax, 1); 
orthogonality_loss2=zeros(kmax, 1);
true_values=eigs(A,kmax);
for j=1:kmax
    W_j(1:j, j) = Q_k(:, 1:j)' * Q_k(:, j);
    for i = 1:j-1
        W_j(j, i) = W_j(i, j); % Copy the upper triangle to the lower triangle
    end
    W_j_1(1:j, j) = Q_k1(:, 1:j)' * Q_k1(:, j);
    for i = 1:j-1
        W_j_1(j, i) = W_j_1(i, j); % Copy the upper triangle to the lower triangle
    end
    W_j_2(1:j, j) = Q_k2(:, 1:j)' * Q_k2(:, j);
    for i = 1:j-1
        W_j_2(j, i) = W_j_2(i, j); % Copy the upper triangle to the lower triangle
    end
    orthogonality_loss1(j)=norm(Q_k1(:,1:j)' * Q_k1(:,1:j) - eye(j), 'inf');
    orthogonality_loss2(j)=norm(Q_k2(:,1:j)' * Q_k2(:,1:j) - eye(j), 'inf');
    orthogonality_loss(j)=norm(Q_k(:,1:j)' * Q_k(:,1:j) - eye(j), 'inf');
    %ortho_loss = norm(W_j(1:j, 1:j) - eye(j), 'fro');
    %fprintf('Iteration %d: Orthogonality Loss = %.2e\n', j, ortho_loss);
    
end
%--------------Ritz values-----------%
figure;
subplot(2, 1, 1);
plot(1:kmax,true_values,'-m',1:kmax, ritz_values1, 'b.-', 1:kmax, ritz_values2, 'r.-',1:kmax, ritz_values, 'g.-');
title('Convergence of Ritz Values');
xlabel('Iteration');
ylabel('Ritz Values');
legend('True values','Lanczos1', 'Lanczos2','LanczosH');

%--------------Orthogonality loss-----------%
subplot(2, 1, 2);
plot(1:kmax, orthogonality_loss1, 'b.-', 1:kmax, orthogonality_loss2, 'r.-',1:kmax, orthogonality_loss, 'g.-');
title('Orthogonality Loss');
xlabel('Iteration');
ylabel('Loss');
legend('Lanczos1', 'Lanczos2','LanczosH');

%--------------Heatmap-----------%
figure;
imagesc(W_j);
colorbar;
title('Matrix W_j LanczosHw');
xlabel('j');
ylabel('j');
figure;
imagesc(W_j_1);
colorbar;
title('Matrix W_j Lanczos1');
xlabel('j');
ylabel('j');
figure;
imagesc(W_j_2);
colorbar;
title('Matrix W_j Lanczos2');
xlabel('j');
ylabel('j');
