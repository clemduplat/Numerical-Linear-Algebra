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
lambda_svd = gcv(U,S,g,'Tikh');
%lambda_svd=gcv(U,S,g,'tikhv');
x_lambda_std = tikhonov(U, S, V, g, lambda_svd);

%---------Advanced: L different I--------%
%Second order derivative to avoid oscillations
n = size(K, 2);
[L_second_derivative, W] = get_l(n, 2);


[U_gen, sm, X_gen, V_gen, W_gen] = cgsvd(K, L_second_derivative); % Using the provided cgsvd function
corner_svd_gen=gcv(U_gen,sm,g,'Tikh');
x_lambda_gen= tikhonov(U_gen, sm, X_gen, g, corner_svd_gen);


%-------------------------------------------------------------------------%
%-----------SVD-based methods-------------%
%-----TSVD------%

%tsvd(U,s,V,b,k)
%find optimal !!!k!!! with a discretized l_curve for the trucnated SVD
[corner_tsvd,G,reg_param]=gcv(U,S,g,'tsvd');
%disp(corner_tsvd); %equal to 37
%Why does it not work ??? OK!
[x_tsvd,rho,eta] = tsvd(U,S,V,g,corner_tsvd);

%------DSVD----%
%not a k, a lambda
corner_dsvd=gcv(U,S,g,'dsvd');
%lambda_dsvd=corner(rho_dsvd,eta_dsvd);
[x_dsvd,rho,eta] = dsvd(U,S,V,g,corner_dsvd);



%------TGSVD-----%
%!!!also a k here!!!!%
corner_tgsvd=gcv(U_gen,sm,g,'tsvd');
%k=corner(rho_tgsvd,eta_tgsvd,reg_param_tgsvd);
[x_tgsvd,rho,eta] = tgsvd(U_gen,sm,X_gen,g,corner_tgsvd);


%-----DGSVD----%
%not a k, a lambda
corner_dgsvd=gcv(U_gen,sm,g,'dsvd');
[x_dgsvd,rho,eta] = dsvd(U_gen,sm,V_gen,g,corner_dgsvd);
%-------------------------------------------------------------------------%
%------Graphical representation------%
% Plotting the solutions with the exact solution f
figure;
plot(f, 'b-'); hold on;
plot(x_lambda_std, 'r--');
title('Solution using Tikihon classic (using L=I)');
legend('Exact Solution', 'Tikihon classic');
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


