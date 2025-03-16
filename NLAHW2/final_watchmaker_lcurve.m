%-------Download informations-----%
load('Temp.mat'); 

%whos('-file', 'Temp.mat');
%K         500x500            2000000  double optained by mid-point rule discretization
%f         500x1                 4000  double exact solution
%g         500x1                 4000  double perturbed right hand side solution
%-------------------------------------------------------------------------%
%-------Inspecting the problem-----%
%condition number
cond_K=cond(K);
disp(['Condition Number of K (if large then the problem is illposed): ', num2str(cond_K)]);
%using mldivide from matlab -> will provide a base-line solution
solution_mldivide=mldivide(K,g);
%-------------------------------------------------------------------------%
%------Applying regularization methods----%
%--------Tikhonov regularization----------%
%-------------Classical: L=I------------%
%I have a problem in this part
L_id = eye(size(K));
[U, S, V] = csvd(K);
corner_classic = l_curve(U, S, g, 'Tikh');
x_lambda_std = tikhonov(U, S, V, g, corner_classic);

%---------Advanced: L different I--------%
n = size(K, 2);


[L_second_derivative, W] = get_l(n, 2);
%L Sobolv%
D = sparse(1:n-1, 2:n, 1, n, n) - sparse(1:n, 1:n, 1, n, n);
L_sobolev = D;
%----Adding some constraint for the problem----%

%scaling_factor = 1e-3;
%Imposing L=id%
%L_second_derivative=L_id;
%L_second_derivative= scaling_factor * L_second_derivative;%imposing it is small
[U_gen, sm, X_gen, V_gen, W_gen] = cgsvd(K, L_second_derivative); % Using the provided cgsvd function
corner_svd_gen=l_curve(U_gen,sm,g,'Tikh');
x_lambda_gen= tikhonov(U_gen, sm, X_gen, g, corner_svd_gen);

[U_class, sm_class, X_class, V_class, W_class] = cgsvd(K, L_id);
corner_svd_class=l_curve(U_class,sm_class,g,'Tikh');
x_lambda_class= tikhonov(U_class, sm_class, X_class, g, corner_svd_class);

%-------------------------------------------------------------------------%
%-----------SVD-based methods-------------%
%-----TSVD------%

%tsvd(U,s,V,b,k)
%find optimal !!!k!!! with a discretized l_curve for the trucnated SVD
[corner_tsvd,rho_tsvd,eta_tsvd,reg_param_tsvd]=l_curve(U,S,g,'tsvd');
%k=corner(rho_tsvd,eta_tsvd,reg_param_tsvd);
[x_tsvd,rho,eta] = tsvd(U,S,V,g,corner_tsvd);

%------DSVD----%
%not a k, a lambda
[corner_dsvd,rho_dsvd,eta_dsvd,reg_param_dsvd]=l_curve(U,S,g,'dsvd');
%[corner_dsvd,rho_dsvd,eta_dsvd,reg_param_dsvd]=l_curve(U_gen,sm,g,'dsvd');
%lambda_dsvd=corner(rho_dsvd,eta_dsvd);
%disp(corner_dsvd);
[x_dsvd,rho,eta] = dsvd(U,S,V,g,corner_dsvd);
%[x_dsvd,rho,eta] = dsvd(U_gen,sm,V_gen,g,corner_dsvd);

%------TGSVD-----%
%!!!also a k here!!!!%
[corner_tgsvd,rho_tgsvd,eta_tgsvd,reg_param_tgsvd]=l_curve(U_gen,sm,g,'tsvd');
%k=corner(rho_tgsvd,eta_tgsvd,reg_param_tgsvd);
[x_tgsvd,rho,eta] = tgsvd(U_gen,sm,X_gen,g,corner_tgsvd);

%-----DGSVD----%
%not a k, a lambda
[corner_dgsvd,rho_dgsvd,eta_dgsvd,reg_param_dgsvd]=l_curve(U_gen,sm,g,'dsvd');
%disp(corner_dgsvd);
[x_dgsvd,rho,eta] = dsvd(U_gen,sm,V_gen,g,corner_dgsvd);
%-------------------------------------------------------------------------%
%---------Conjugate Graident-Based Method--------%
%Does not include an L matrix + convergence can be delayed due to roundoff
%errors
reorth = 0; 
k=20;
[X_classic,rho,eta,F] = cgls(K,g,k,reorth,S);


%-------------------------------------------------------------------------%
%------Graphical representation------%
% Plotting the solutions with the exact solution f

figure;
plot(f, 'b-'); hold on;
plot(x_lambda_std, 'r--');
title('Solution using Tikihon classic without L');
legend('Exact Solution', 'Tikihon');
hold off;

figure;
plot(f, 'b-'); hold on;
plot(x_lambda_class, 'r--');
title('Solution using Tikihon classic (using L=I)');
legend('Exact Solution','Tikhonov');
hold off;

figure;
plot(f, 'b-'); hold on;
plot(x_lambda_gen, 'r--');
title('Solution using Tikihon advanced');
legend('Exact Solution','Tikhonov advanced');
hold off;

figure;
plot(f, 'b-'); hold on;
plot(x_tsvd, 'r--');
title('Solution using TSVD');
legend('Exact Solution', 'TSVD Solution');
hold off;

figure;
plot(f, 'b-'); hold on;
plot(x_dsvd, 'r--');
title('Solution using DSVD');
legend('Exact Solution', 'DSVD solution');
hold off;

figure;
plot(f, 'b-'); hold on;
plot(x_tgsvd, 'r--');
title('Solution using TGSVD');
legend('Exact Solution', 'TGSVDsolution');
hold off;

figure;
plot(f, 'b-'); hold on;
plot(x_dgsvd, 'r--');
title('Solution using DGSVD');
legend('Exact Solution', 'DGSVD solution');
hold off;

figure;
plot(f, 'b-'); hold on;
plot(X_classic(:,k), 'r--');
title('Solution using CGLS with k=30');
legend('Exact Solution', 'CGLS');
hold off;


