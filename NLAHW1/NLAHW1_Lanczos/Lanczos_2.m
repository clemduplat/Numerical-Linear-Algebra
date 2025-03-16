function [Q_k,T_k,r,err_ind] = Lanczos_2(A,kmax,r,nrm_A)
 

% INIT
n=size(A,1);
eta = ((eps)^(3/4))/sqrt(kmax);   % intermed. orth level
delta = sqrt(eps/kmax);         % threshold \delta for semi-orthogonality
reorth_prev = 1;                % RECOMMENDED: reorth. against prev 2 vectors
                                % you can turn this off, but results can be
                                % quite bad
err_ind = 0;                    % output error flag

%roundoff err. estim.
eps1 = sqrt(n)*eps/2;           % eps_mach = eps/2 !!! (in eta and delta 2*eps_mach is used)
gamma = 1/sqrt(2);		 % factor to check sufficient reduction
% Alloc. 
alpha = zeros(kmax+1,1);  beta = zeros(kmax+1,1);
Q_k = zeros(n,kmax);
q = zeros(n,1); beta(1)=norm(r);
full_reorth=1;			%toggle when needed
w = zeros(n,1); % Initialize loss of orthogonality vector
w_old = zeros(n,1); % Initialize old loss of orthogonality vector
for j=1:kmax
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%    PART ONE: REGULAR LANCZOS     %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  q_old = q;
  if beta(j)==0
    q = r;
  else
    q = r / beta(j);
  end
  Q_k(:,j) = q;
  u = A*q;
  r = u - beta(j)*q_old;
  alpha(j) = q'*r;
  r = r - alpha(j)*q;
  beta(j+1) = norm(r);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% PART TWO: REORTH w.r.t. j & j-1  %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  if beta(j+1)<gamma*beta(j) && reorth_prev 
    if  j==1
      h=0;
      for i=1:2
          dh = q'*r;    
          r = r-q*dh;
          h = h+dh;
      end
      alpha(j) = alpha(j) + h;
    elseif j>1
      h1 = q_old'*r;
      h2 = q'*r;
      r = r  - (q_old*h1 + q*h2); 
      if beta(j)~=0               
        beta(j) = beta(j) + h1;
      end
      alpha(j) = alpha(j) + h2;
    end        
    beta(j+1) = norm(r);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%  PART THREE: COMPUTE ORTH. LOSS  %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % In this part you should calculate loss of orth.
  
  % Lanczos1.m                  : compute w_{j,\infty}
  % Lanczos2.m & Lanczos3.m     : compute w vec. through recurrence (see 'update_w' below)
  if j > 1 
     [w, w_old] = update_w(w, w_old, alpha, beta,j);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%        PART FOUR: REORTH.        %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  % In this part you should reorthogonalize whenever loss. of orth. is
  % detected
  
  % Lanczos1.m & Lanczos2.m     : full reorth (see boolean full_reorth)
  % Lanczos3.m                  : select L and orth w.r.t. vectors in L
  %                               DON'T FORGET: update w afterwards!
  % DON'T FORGET: in next run of main loop, reorth against L (or full) again!!!
  % Write your reorth in a function of the form 
  % [r,nrmnew] = reorth(Q_k,r,nrm,L), with nrm the norm of r on input, nrmnew norm on output;
  if full_reorth || (norm(w)>eta)
    [r,nrmnew]=reorth(Q_k,r,norm(r),1:j); 
    beta(j+1)=nrmnew; % Update beta 
  end
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%   PART FIVE: DEFLATE & ESCAPE    %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % don't worry about this part, can be omitted, but sometimes needed if
  % eigenvalues to be computed are too small or if you are unlucky
  
  if  j<kmax && beta(j+1) < n*nrm_A*eps  
    % beta too small: attempt to deflate and escape to new invariant
    % subspace
    beta(j+1) = 0;
    FAIL = 1;
    for try_ind=1:3    
      r = rand(n,1)-0.5;  
      r = A*r;
      nrm=norm(r); % not necessary to compute the norm accurately here.
      L = 1:j;
      [r,nrmnew] = reorth(Q_k,r,nrm,L);
      w(L) = eps1;
      if nrmnew > 1e-14
        % SUCCES!
        % Continue iteration.
        FAIL=0;
        break;
      end
    end
    if FAIL
      err_ind = -j; % negative error index: algorithm terminated because unsolvable failure at j
      break;
    else
      r=r/nrmnew; % Continue with new r
      orth_jm1_again = 1;
    end    
  elseif j<kmax && ~full_reorth && beta(j+1)*delta < nrm_A*eps1
    % recurrence no longer reliable
    % maybe eigenvalues to be approx.'d too small?
    % switch to full reorth
    reorth_prev = 1;
    err_ind = j;                % positive error index: means from that index
                                % going forth, full reorth is used
  end
end

T_k = spdiags([[beta(2:j);0] alpha(1:j) beta(1:j)],-1:1,j,j);


end
function [w, w_old] = update_w(w, w_old, alpha, beta, j)
    %Then, with $w_{j0}:=0$:
    %\begin{alignat*}{3}
    %&w_{kk}&&=1&& \text{ for }k=1,\ldots,j\\
    %&w_{kk-1}&&=\underbrace{\mathbf{q}_k^{\ast}\mathbf{q}_{k-1}}_{:=\psi_k} && \text{ for } k=2,\ldots,j\\
    %&w_{jk+1}&&=w_{k+1j}&& \text{ for } k=1,\ldots,j-1\\
    %&\beta_{j}w_{j+1k}&&=\beta_{k}w_{jk+1}+&&(\alpha_k-\alpha_j)w_{jk}+\beta_{k-1}w_{jk-1}-\beta_{j-1}w_{j-1k}+\underbrace{\mathbf{q}_{j}^{\ast}\mathbf{f}_k-\mathbf{q}_{k}^{\ast}\mathbf{f}_j}_{:=\theta_{j,k}}
    %\end{alignat*}
    %where the last equation holds for $k=1,\ldots,j-1$.
    %You can take $\theta_{jk}=\epsilon(\beta_{j}+\beta_{k})z_1$ with $z_1\sim \mathcal{N}(0,.3)$ 
    % and $\psi_{k}=\sqrt{n}\epsilon(\beta_{1}/\beta_{j}) z_2$ with $z_2\sim \mathcal{N}(0,.6)$. 
    % You can update the selected elements in $w$ to $\epsilon z_3$, with $z_3\sim\mathcal{N}(0,3/2)$. 
    % You do not need to worry where these numbers come from.
    length_w = length(w); 
    epsilon = eps;  
    theta = zeros(length_w, 1);
    psi = zeros(length_w, 1);

    for k = 1:min(j, length_w)
        theta(k)= epsilon*(beta(k)+beta(j+1))*normrnd(0, 0.3); 
    end

    for k = 2:min(j, length_w)
        psi(k)= sqrt(length_w)*epsilon*(beta(1)/beta(k))*normrnd(0, 0.6); 
    end
    w_new= w;
    w_new(2:end)= w_old(1:end-1);
    w_new(1)= 0;
    z_3= epsilon*normrnd(0,1.5);
    %-------- Update w using recurrence---------%
    for k = 1:length_w
        if k == 1
            w(k) = 1;  % w_kk = 1 for k=1
        elseif k==j
            w(k)=psi(k);
        else 
            if k < j
                % ----recurrence formula---- %
                w(k)= (beta(k)* w_new(k+1) + (alpha(k) - alpha(j)) * w_old(k) + ...
                        beta(k-1)* w_old(max(k-1,1)) - beta(j-1)* w_old(max(j-1,1)) + theta(k))*z_3 / beta(j);
            end
        end
    end

    % Update w_old for the next iteration
    w_old = w_new;
end
function [r,nrmnew] = reorth(Q_k,r,nrm,L)
    for i=L
        coef=Q_k(:,i)' * r;  % Compute the coefficient for the i-th basis vector
        r=r - coef*Q_k(:,i); % Subtract the projection of r onto the i-th basis vector
    end
    nrmnew=norm(r); % Compute the new norm after reorthogonalization
end



