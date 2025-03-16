function [Cr, Gr, X] = prima(C, G, B, m, verbose)
%  [Cr, Gr, X] = prima(C, G, B, m)
% Compute the congruence transformed matrices from the arnoldi vectors
% as in Altan's paper
%
%  C dx/dt = G x + Bu
%
%  m*(columns of B) is the order of the reduced model
%
%  verbose:  Optional argument for extra fprintfs (1 yes , 0 no)

if (nargin < 5)
  verbose = 0;
end

block_size = size(B,2);

% G is assumed sparse and invertible

fprintf(1,'doing invG_CB\n');
%invG_CB = G\[C B];
%invG_C = invG_CB(:, 1:size(C,2));
%invG_B = invG_CB(:, size(C,2) + (1:size(B,2)));

% let's load it instead
invG_C = G\C;
invG_B = G\B;

%load /tmp/invPR

%invG_C = invPR_L;
%invG_B = invPR_B;

fprintf(1,'calling barnoldi\n');
[X, H, R1] = barnoldi(invG_C, invG_B, m, verbose);

fprintf(1,'computing reduced system\n');
Cr = X'*C*X;
Gr = X'*G*X;




