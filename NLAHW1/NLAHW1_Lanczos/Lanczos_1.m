function [Q_k,T_k,r,err_ind] = Lanczos_1(A,kmax,r,nrm_A)
 

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
  
  Wj = Q_k(:,1:j)' * Q_k(:,1:j); 
  Ij = eye(size(Wj));            
  w_inf=norm(Wj-Ij,'Inf');
  
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
  
  %full_reorth || 
  if (w_inf>eta)
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

function [r,nrmnew] = reorth(Q_k,r,nrm,L)
%-------Full orthogonalization----%
    for i=L
        coef= Q_k(:,i)' *r;  % Compute the coefficient for the i-th basis vector
        r=r - coef*Q_k(:,i); % Subtract the projection of r onto the i-th basis vector
    end
    nrmnew=norm(r); % Compute the new norm after reorthogonalization
end


