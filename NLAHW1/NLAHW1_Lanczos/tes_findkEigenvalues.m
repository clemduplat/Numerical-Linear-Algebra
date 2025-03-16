%-------test with simple matrices-----%

%[Q_k1, T_k1, r1, err_ind1]=Lanczos_1(A,kmax,r,nrm_A);
tol=1e-10;
k=4;

B = [2, 1, 0, 0;
     1, 4, 1, 0;
     0, 1, 5, 1;
     0, 0, 1, 3];


kth_eigenvalue = findKthEigenvalue(B, k, tol);
disp(kth_eigenvalue);
firstk_eigenvalues = findFirstKEigenvalues(B, k, tol);
disp(firstk_eigenvalues);
%------Working well------%