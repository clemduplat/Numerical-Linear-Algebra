function  [U, H, R1] = barnoldi(A, B, m, verbose)
%barnoldi [U, H, R1] = barnoldi(A, B, M),
%	Block Arnoldi algorithm
%            A  matrix n x n
%            B  block of rhs  n x s
%            m  Krylov subspace dimension wanted
%            verbose  Optional argument for extra fprintfs (1 yes , 0 no)
%
%            U = [U1 U2 U3 ...] matrix of orthogonal vectors
%            H - same as arnoldi H  but block Hessenberg
%            R1 - initial R from QR factorization such that  B = U1*R1;
%	QR orthogonalization has been used (matlab function qr)

if (nargin < 4)
  verbose = 0;
end

[n, s] = size(B);
H = zeros((m+1)*s,m*s);
U = zeros(n,m*s);

%% alocate space
sm1 = s * (m+1);
sm  = s * m;
%
% matlab 4.2c seems to want this to be full.  matlab 5 takes forever if it
%  isn't.  Although this was tested only for a matrix of one column
%[U(:, 1:s), R1] = qr(B,0);
[U(:, 1:s), R1] = qr(full(B),0);
  
for (j = 1:m)
  %% set the indexes
  jms = (j-1) * s + 1;
  js  = j * s;
  j1s = (j+1) * s;
  js1 = js + 1;
  
  Up = A * U(:, jms:js);
  
  if (verbose == 1 & rem(j,5) == 1)
    fprintf(1,'Done with %d block-matvec\n',j);
  end
  
  %% new basis block
  for (i = 1:j)
    ims = (i-1)*s + 1;
    is = i*s;
    
    H(ims:is, jms:js) = U(:, ims:is)' * Up;
    Up = Up - U(:, ims:is) * H(ims:is, jms:js);
  end;

  [U(:, js1:j1s), w] = qr(Up,0);
  H(js1:j1s, jms:js) = w;
  
  %% updating R from previous rotations
%  R(1:j1s, jms:js) = R(1:j1s, jms:js) + Q(1:j1s, 1:j1s) * H(1:j1s, 1:s);
  
  %% rotations to anihilate new elements
%  for (jj = 1:s)
%    k1 = jms - 1 + jj;
%    for (ii = k1+s:-1:k1+1)
%      GG = givens(R(ii-1, k1), R(ii, k1));
%      Q(ii-1:ii, 1:j1s) = GG * Q(ii-1:ii, 1:j1s);
%      R(ii-1:ii, k1:js) = GG * R(ii-1:ii, k1:js);
%    end;
%  end;
end;

fprintf(1,'%f mat-vecs performed to produce %f x %f H matrix\n',[s*m j1s js]);




