%----------------Load data-----------%
load('HW1.mat');

% ---------Plot the full AFM mesh-------------%
figure(1)
plot_afm_mesh(nodes, elements);


%-----Convert sparse matrix to full matrix----%
fullK = full(K); 
fullM = full(M);


%------Perform Lanczos Algorithm for tridiagonalization----%
% Simplified eigenvalue problem with matrix A=K^-1*M
kmax=1000;
A=fullM\fullK;


[n,m]=size(A);
r = randn(n, 1);  
%r=r/norm(r);
nrm_A = norm(A);
[Q_k,T_k,r,err_ind] = Lanczos_2(A,kmax,r,nrm_A);

%-------Find first k eigenavalues----%
%5 dominant eigenvalues 
k=5;
tol=1e-2;
k_eigenvalues = findFirstKEigenvalues(T_k, k, tol);



%----- k eigenmodes-----%
T_k=full(T_k);
[V_Tk, e] = eig(T_k);
eigenmodes = Q_k*V_Tk;
k_eigenmodes=eigenmodes(:,1:k);

for i = 1:k
    % -----displacement vector for the i-th mode-------%
    mode = k_eigenmodes(:, i);
    dU = zeros(size(nodes, 1) * 3, 1); 
    dU(actualDofs) = mode; 
    dU = reshape(dU, [], 3); 
    
    %--------- Amplify the deformations-------%
    deformationScale = 200; 
    displaced_nodes = nodes + dU * deformationScale;
    
    %-----Plot of the deformed cantilever for i-th mode----%
    figure;
    plot_afm_mesh(nodes, elements); 
    hold on
    plot_afm_frame(displaced_nodes, elements); 
end

% Report your findings
% ...
%-----Compare------%
%--------Calculated eigenvalues----%
disp('Calculated eigenvalues:');
disp(k_eigenvalues);
%---------True eigenvalues-------%
d = eigs(T_k, k, 'smallestabs');
disp('True eigenvalues using eigs:');
disp(d);