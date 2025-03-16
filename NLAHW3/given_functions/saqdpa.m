function [poles, residues, rightev, leftev, nr_solves] = ...
           saqdpa(A, E, M, b, c, d, s0, options)
% [poles, residues, leftev, rightev, nr_solves] = ...
%     saqdpa(A, E, M, b, c, d, s0, options)
% computes the dominant poles of the quadratic transfer function
%               H(s) = c'*(s^2M + sE + A)^{-1}b + d
% of the second-order system
%               M*x'' + E * x' + A*x = b*u
%                                 y  = c' * x + d * u
% using the Subspace Accelerated Quadratic Dominant Pole Algorithm of 
% Joost Rommes and Nelson Martins. A second-order system yields a 
% Quadratic Eigenvalue Problem (QEP) of the form
%               (s^2 * M + s * E + A)*x = 0, 
% M, E, and K being matrices of size n x n.
% This algorithm is specifically designed for large sparse 
% matrices A, E and M, that allow *cheap* solves  using lu or \.
%
%
% Input:
%      A,M,E: (nxn) system matrices
%      b: (nx1) input vector
%      c: (nx1) output vector
%      d: (1x1) direct term of system
%      s0: (kx1) initial estimates for pole
%          if size(s0,1) > 1, s0(i) will be injected after the
%              the i-th pole has been found (ie, one solve with
%              (s0(i)*E-A) in the next iteration)
%      options: structure with
%         nwanted : number of wanted poles (5)
%         tol  : tolerance for eigenpair residual (1e-10)
%         displ: show convergence information (1)
%         selection strategy: select pole p and residue with maximal ('LM'):
%                          one of
%                           'LM' : largest |res| magnitude
%                           'LS' : largest |res|/|p| magnitude
%                           'LR' : largest |res|/|re(p)| magnitude
%         kmin: minimal searchspace dimension (min(1,n))       
%         kmax: maximal searchspace dimension (min(10,n))
%         maxrestart: maximal nr of restarts (100)
%
%         Function handles to user functions for MVs and solves
%         can be used instead of A and E:
%         f_ax    :   Y = f_ax( tr, X) should return
%                         A*X if tr = 'N'
%                         A'*X if tr = 'T'
%         f_ex    :   Y = f_ex( tr, X) should return
%                         E*X if tr = 'N'
%                         E'*X if tr = 'T'
%         f_mx    :   Y = f_mx( tr, X) should return
%                         M*X if tr = 'N'
%                         M'*X if tr = 'T'
%         f_semax :   Y = f_semax( tr, X, s) should return
%                         (s^2M + sE + A)*X  if tr = 'N'
%                         (s^2M + sE + A)'*X if tr = 'T'
%         f_semax_s : Y = f_semax_s( tr, X, s) should return
%                         inv(s^2M + sE + A) * X  if tr = 'N'
%                         inv((s^2M + sE + A)') * X if tr = 'T'
%         Note that when using function handles, still dummy A, E, M need
%         to be supplied as arguments to sadpa_adi, typically 
%         just A=E=M=sparse(n,n)
%
%         Advanced options:
%         use_lu: if use_lu==1, use LU for solves (default), otherwise \ (0)
%         use_lu_w_amd : if 1 use complete UMFPACK reordering
%         dpa_bordered: if dpa_bordered==1, use bordered matrices in solves
%                           ie [s^2M + sE + A, -b ; c', d] (default with use_lu=0)
%                       else use just s^2M + sE + A (default with use_lu=1) 
%         yEx_scaling: if yEx_scaling == 1, scale approximate eigenvectors
%                         during selection such that y'Ex = 1
%                      else scale ||x|| = ||y|| = 1 (default)
%         rqitol: tolerance on ||r|| = ||Ax + lambda Ex + lambda^2Mx||. If ||r||<rqitol,
%                 refine using twosided Rayleigh quotient iteration (1e-4). Note
%                 rqitol > tol
%         turbo_deflation: if turbo_deflation == 1, use b and c to efficiently deflate
%                             (default, see Chapter 3 of PhD Thesis Joost Rommes)
%                          else use classical deflation via explicit
%                               projections with found eigenvectors
%.         
% Output:
%      poles: converged poles
%      residues: residues needed to construct bode plot
%      leftev: left eigenvectors
%      rightev: right eigenvectors
%      nr_solves: number of LU factorizations used
%
% When using this software, please cite
%
%   Joost Rommes and Nelson Martins
%   Efficient computation of transfer function dominant poles of 
%   large second-order dynamical systems
%   SIAM Journal on Scientific Computing 30 (4), 2008, pp. 2137-2157
%
%
% Joost Rommes (C) 2006 -- 2010
% rommes@gmail.com
% http://sites.google.com/site/rommes

if nargin < 8
    options = struct([]) ;
end

[n, nwanted, tol, strategy, use_lu, use_lu_w_amd, kmin, krestart, maxrestarts, ...
    displ, dpa_bordered, yEx_scaling, rqitol, tdefl, ...
    f_ax, f_ex, f_mx, f_semax, f_semax_s] = ...
     init( A, options ) ;

if isa( f_semax_s, 'function_handle' ) && dpa_bordered
    error( 'Bordered DPA does not support function handles for solves' ) ;
end

b_f_ax = isa( f_ax, 'function_handle' ) ;
b_f_ex = isa( f_ex, 'function_handle' ) ;
b_f_mx = isa( f_mx, 'function_handle' ) ;
b_f_semax = isa( f_semax, 'function_handle' ) ;
b_f_semax_s = isa( f_semax_s, 'function_handle' ) ;


conjtol = tol * 10^2 ;
imagtol = 1e-8 ;

nr_shifts = length(s0) ;
shifts = s0 ; %pole approximations
shift_nr = 2 ;
just_found = 0 ;

k = 0 ;
rhs = sparse([ zeros(n, 1); 1 ] );

X = [] ;
V = [] ;

nr_converged = 0 ;
poles = [] ;
leftev = [] ;
rightev = [] ;
nr_solves = 0 ;

Q = zeros(n,0) ;
Qt= zeros(n,0) ;
Qgen = zeros(2*n,0) ;
Qgent = zeros(2*n,0) ;
Qsparse = zeros(2*n,0) ;
Qtsparse = zeros(2*n,0) ;
Bg_K = blkdiag( -A, M ) ;

R = [] ;
Rt= [] ;
G = [] ;
H = [] ;
AX = [] ;
MX = [] ;
K = [] ; %reduced problem lambda^2 K + lambda G + H

st = shifts(1) ; %s(1,1) ;
st_newton = 0 ;
nrestart = 0 ;

bt_d = [ zeros(n,1) ; b ] ;
ct_d = [c ; zeros(n,1) ] ;

if b_f_semax_s == 0
    bK_d =  (A\ bt_d(1:n,1)) ; 
    cK_d =  (A'\c) ;
else
    bK_d =  f_semax( 'N', bt_d(1:n,1), 0 ) ; 
    cK_d =  f_semax( 'T', c, 0) ;    
end

while nrestart < maxrestarts
    k = k+1 ;

    if just_found && shift_nr < nr_shifts
        just_found = 0 ;
        st = shifts( shift_nr ) ;
        shift_nr = shift_nr + 1 ;
    end


    if tdefl == 0 %dpa standard
        if dpa_bordered %original dpa form
	        sEmA = [st^2*M + st * E + A, -b; c', d ] ;    

            if use_lu
                [xu,vu] = solve_with_lu( sEmA, b, c, use_lu_w_amd ) ;
                %[L,U] = lu( sEmA ) ;
                %xu = U \ ( L \ rhs ) ;    
                %vu = L' \ ( U' \ rhs ) ;
            else
                xu = sEmA \ rhs ;    
                vu = sEmA' \ rhs ;
            end
            xu = xu(1:n) ;
            vu = vu(1:n) ;
        else %dpa variant without border (exact Newton)
            if b_f_semax_s == 0
                sEmA = st^2*M + st * E + A ;    
    
    	        if use_lu
                    [xu,vu] = solve_with_lu( sEmA, b, c, use_lu_w_amd ) ;
                    %[L,U] = lu( sEmA ) ;
                    %xu = U \ ( L \ b ) ;
                    %vu = L' \ ( U' \ c ) ;
                else
                    xu = sEmA \ b ;
                    vu = sEmA' \ c ;
                end
            else
                xu = f_semax_s( 'N', b, st ) ;
                vu = f_semax_s( 'T', c, st ) ;
            end
        end

	    if b_f_mx == 0
            Mxu = M*xu ;
        else
            Mxu = f_mx( 'N', xu ) ;
        end
	    if b_f_ex == 0
            Exu = E*xu ;
        else
            Exu = f_ex( 'N', xu ) ;
        end
        st_newton = st - ( (c'*xu) / (vu'*(2*st*(Mxu) + Exu))) ;
	
    else %turbo deflation
        if b_f_semax_s == 0
            sEmA = st^2*M + st * E + A ;    
   	        if use_lu
                [xu2,vu2] = solve_with_lu( sEmA, st*bt_d(n+1:end,1) + bt_d(1:n,1),...
                              st'*ct_d(n+1:end,1) + ct_d(1:n,1), use_lu_w_amd ) ;
	    	    %[L,U] = lu( sEmA ) ;
		        %xu2 = U \ (L \ (st*bt_d(n+1:end,1) + bt_d(1:n,1)) ) ;
        		%vu2 = L' \ (U' \ (st'*ct_d(n+1:end,1) + ct_d(1:n,1))) ;
	        else
                xu2 = sEmA \ (st*bt_d(n+1:end,1) + bt_d(1:n,1)) ;
		        vu2 = sEmA' \ (st'*ct_d(n+1:end,1) + ct_d(1:n,1)) ;
	        end
        else
            xu2 = f_semax_s( 'N', st*bt_d(n+1:end,1) + bt_d(1:n,1), st ) ;
            vu2 = f_semax_s( 'T', st'*ct_d(n+1:end,1) + ct_d(1:n,1), st ) ;
        end
		xu = ( -bK_d + xu2 ) / st ;
	    vu = ( -cK_d + vu2 ) / (st') ;        
    end %end if tdefl == 0

    nr_solves = nr_solves + 1 ;
    
    x = xu(1:n) ;

    if tdefl == 0
		tmp = [ x ; st_newton*x ] ;
		tmp = mgs_skew( Qgen, Qtsparse, tmp, 1, 1 ) ;
	else
		tmp = x ;
	end
	[x,h] = mgs(X, tmp(1:n,1),1) ;
	X(:, k) = x / h(end) ;

    v = vu(1:n) ;	

	if tdefl == 0
		tmp = [v ; st_newton'*v ] ;
		tmp = mgs_skew( Qgent, Qsparse, tmp, 1, 1) ;
	else
		tmp = v ;
    end
	[v,h] = mgs(V, tmp(1:n,1),1) ;
    V(:, k) = v / h(end) ;
    
    if b_f_ax == 0
        AX = [AX, A*X(:, k)] ;
        H = [H, V(:,1:k-1)'*AX(:,k) ; V(:,k)'*AX] ;
    else
        AX = [AX, f_ax('N', X(:, k))] ;
        H = [H, V(:,1:k-1)'*AX(:,k) ; V(:,k)'*AX] ;
    end

    if b_f_ex == 0
        G = [ G, V(:,1:k-1)'*E*X(:,k) ; V(:,k)'*E*X ] ;
    else
        vE = f_ex( 'T', V(:,k) ) ;
        G = [ G, V(:,1:k-1)'*f_ex('N',X(:,k)) ; vE'*X ] ;
    end

    if b_f_mx == 0
        K = [ K, V(:,1:k-1)'*M*X(:,k) ; V(:,k)'*M*X ] ;
    else
        vM = f_mx( 'T', V(:,k) ) ;
        K = [ K, V(:,1:k-1)'*f_mx('N',X(:,k)) ; vM'*X ] ;
    end
    
    [SH, SHt, UR, URt] = select_max_res( H, G, K, X, V, b, c, ...
                                 strategy, bt_d, ct_d, tdefl, ...
                                 A, E, M, f_ex, f_mx, yEx_scaling) ;
    
    found = 1 ;
    do_rqi = 1 ;
    while found
        schurvec = X * UR(:,1) ;
        theta = SH(1,1) ;
        schurvec = schurvec / norm(schurvec) ;
        lschurvec = V * URt(:,1) ;
        lschurvec = lschurvec / norm(lschurvec) ;
        st = theta ;

        if b_f_semax == 0
            r = A*schurvec + theta * (E* schurvec) + theta^2*(M*schurvec) ;
        else
            r = f_semax( 'N', schurvec, theta ) ;
        end
        nr = norm(r) ;
        
        if displ
            fprintf('Iteration (%d): theta = %s, norm(r) = %s\n', k, num2str(theta),num2str(nr) ) ;
        end
        
        %if the direct rqi residual is smaller, then use that approx.
        %this avoids unnecessary stagnation due to rounding errors        
        if nr < rqitol && do_rqi 
			rqopt.x0 = schurvec ;
			rqopt.v0 = lschurvec ;
			rqopt.sigma = theta ;
			rqopt.tol = tol ;
            [x_rq, x_rqi, v_rqi, bla, blabla, nr_solvesp] = ...
                        rqi_quad(A, E, M, use_lu, use_lu_w_amd,...
                        f_ax, f_ex, f_mx, f_semax, f_semax_s, rqopt) ;
			nr_solves = nr_solves + nr_solvesp ;
            v_rq = x_rq' ;
			theta = x_rq ;
			if b_f_semax == 0
                rqi_res = A*x_rqi + theta * (E* x_rqi) + theta^2*(M*x_rqi) ;
            else
                rqi_res = f_semax( 'N', x_rqi, theta ) ;
            end

			nrq = norm(rqi_res) ;
			if nrq < nr %& nrq < tol
                schurvec = x_rqi ;
                lschurvec = v_rqi ;
                theta = x_rq ;
                nr = nrq ;
            end

            do_rqi = 0 ; %to avoid duplicate convergence if more than one found
            
            if abs(imag( theta ) )/abs(theta) < imagtol
            %possible real pole
			   if b_f_semax == 0
                   rres = A*real(schurvec) + real(theta) * (E*real(schurvec)) + ...
                          real(theta)^2*(M*real(schurvec)) ;
               else
                   rres = f_semax( 'N', real(schurvec), real(theta) ) ;
               end
               nrr = norm( rres) ;
               if nrr < nr
                   schurvec = real(schurvec) ;
                   lschurvec = real(lschurvec) ;
                   theta = real(theta ) ;
                   nr = nrr ;
                   %st = theta ;
               end
            end
        end %end rqi
        
        found = nr < tol && nr_converged < nwanted ;
        
        if found
            Q = [Q, schurvec] ;
            Qgen = [Qgen, [schurvec ; theta*schurvec] ] ;
			Qsparse = [Qsparse, Bg_K*Qgen(:,end) ] ;
            R = [R, theta] ;
            

            Qt = [Qt, lschurvec] ;
            Qgent = [Qgent,[lschurvec; conj(theta)*lschurvec] ] ;
			Qtsparse = [Qtsparse, Bg_K'*Qgent(:,end)] ;
            Rt = [ Rt, conj(theta)] ;
            
			nqqt = Qtsparse(:,end)'* Qgen(:,end) ;
            Q(:,end) = Q(:,end) / nqqt ;
			Qgen(:,end) = Qgen(:,end) / nqqt ;
            Qsparse(:, end) = Qsparse(:, end) / nqqt ;

            poles = [ poles, theta ] ;

            if displ
                fprintf( 'Dominant pole %s found in iteration %d\n', num2str(theta), k ) ;
            end
            nr_converged = nr_converged + 1 ;

            J = 2 : k ;
            X = X * UR(:, J) ; 
            V = V * URt(:, J) ;
			
			if tdefl
				bt_d = bt_d - Qsparse * (Qgent'*bt_d) ;
				ct_d = ct_d - Qtsparse * (Qgen'*ct_d) ;

                if b_f_semax_s == 0
                    bK_d =  (A\ bt_d(1:n,1)) ; 
                    cK_d =  (A'\ct_d(1:n,1)) ;
                else
                    bK_d =  f_semax( 'N', bt_d(1:n,1), 0 ) ; 
                    cK_d =  f_semax( 'T', ct_d(1:n,1), 0) ;    
                end
			end

            if abs(imag( theta ) )/abs(theta) < imagtol %relative
            %in case of real pole, do reortho. If a complex pair
            %is detected, do ortho after defl of complex conj
				for kk = 1 : size(V,2)
					tmp = [ X(:,kk) ; SH(kk,kk)*X(:,kk) ] ;
					tmp = mgs_skew( Qgen, Qtsparse, tmp, 1, 1 ) ;
					[X(:, kk), vv] =  mgs(X(:,1:kk-1), tmp(1:n,1)) ;
					X(:, kk) = X(:, kk) / vv(end) ;
					tmp = [ V(:,kk) ; SHt(kk,kk)*V(:,kk) ] ;
					tmp = mgs_skew( Qgent, Qsparse, tmp, 1, 1 ) ;
					[V(:, kk), vv] =  mgs(V(:,1:kk-1), tmp(1:n,1)) ;
					V(:, kk) = V(:, kk) / vv(end) ;
				end
            end
            k = k - 1 ;
            
            %conjugate pair check
            cce = conj( theta ) ;
            if abs(imag( cce )) / abs(cce) >= imagtol %relative
                pairdifs = abs( poles - cce ) ;
                if min( pairdifs ) > tol
                    %possible conjugate pair
                    %directly deflate conjugate
                    ccv = conj( schurvec ) ;
                    ccv = ccv / norm(ccv) ;

                    if b_f_semax == 0
                        r = A*ccv + cce*E*ccv +cce^2*M*ccv;
                    else
                        r = f_semax( 'N', ccv, cce ) ;
                    end
                    if displ
                        fprintf('(conjugate) Iteration (%d): directly deflate conj...\n', k ) ;
                    end

                    if norm(r) < conjtol
                        Q = [Q, ccv] ;
						Qgen = [Qgen, [ccv ; cce*ccv] ] ;
						Qsparse = [Qsparse, Bg_K*Qgen(:,end) ] ;
                        R = [ R, cce ] ;
                           
                        ccet = conj(cce) ;
                        ccvt = conj( lschurvec ) ;

                        Qt = [Qt, ccvt];
						Qgent = [Qgent,[ccvt; ccet*ccvt] ] ;
						Qtsparse = [Qtsparse, Bg_K'*Qgent(:,end)] ;
                        Rt = [ Rt, ccet ] ;
                        
						nqqt = Qtsparse(:,end)'* Qgen(:,end) ;
                        Q(:,end) = Q(:,end) / nqqt ;
                        Qgen(:,end) = Qgen(:,end) / nqqt ;
						Qsparse(:, end) = Qsparse(:, end) / nqqt ;

                        poles = [poles, cce] ;
                        if displ
                            fprintf( '(conjugate) ...dominant pole %s found in iteration %d\n', num2str(cce), k ) ;
					    end
                        if tdefl
							bt_d = bt_d - Qsparse * (Qgent'*bt_d) ;
							ct_d = ct_d - Qtsparse * (Qgen'*ct_d) ;

							if b_f_semax_s == 0
                                bK_d =  (A\ bt_d(1:n,1)) ; 
                                cK_d =  (A'\ct_d(1:n,1)) ;
                            else
                                bK_d =  f_semax( 'N', bt_d(1:n,1), 0 ) ; 
                                cK_d =  f_semax( 'T', ct_d(1:n,1), 0) ;    
                            end
						end
                           
                        for kk = 1 : size(V,2)
						    tmp = [ X(:,kk) ; SH(kk,kk)*X(:,kk) ] ;
                            tmp = mgs_skew( Qgen, Qtsparse, tmp, 1, 1 ) ;
                            [X(:, kk), vv] =  mgs(X(:,1:kk-1), tmp(1:n,1)) ;
                            X(:, kk) = X(:, kk) / vv(end) ;
						    tmp = [ V(:,kk) ; SHt(kk,kk)*V(:,kk) ] ;
                            tmp = mgs_skew( Qgent, Qsparse, tmp, 1, 1 ) ;
                            [V(:, kk), vv] =  mgs(V(:,1:kk-1), tmp(1:n,1)) ;
                            V(:, kk) = V(:, kk) / vv(end) ;
                        end
                    end                    
                end
            end
            if b_f_ex == 0
                G = V' * E * X ;
            else
                G = V' * f_ex( 'N', X ) ;
            end
            if b_f_ax == 0
                AX = A*X ;
            else
                AX = f_ax( 'N', X ) ;
            end
            H = V' * (AX) ;

            if b_f_mx == 0
                MX = M*X ;
            else
                MX = f_mx( 'N', X ) ;
            end
            K = V' * (MX) ;

            if k > 0
                [SH, SHt, UR, URt] = select_max_res( H, G, K, X, V, b, c, ...
                                 strategy, bt_d, ct_d, tdefl,...
                                 A, E, M, f_ex, f_mx, yEx_scaling) ;
            end

            found = k > 0 ;
        elseif k >= krestart
            %do a restart by selecting the kmin most dom triplets
            k = kmin ;
            K = 1:kmin ;
            X = X * UR(:, K) ; 
            V = V * URt(:, K) ;

            for kk = 1 : size(V,2)
                if tdefl == 0
					tmp = [ X(:,kk) ; SH(kk,kk)*X(:,kk) ] ;
					tmp = mgs_skew( Qgen, Qtsparse, tmp, 1, 1 ) ;
                else
					tmp = X(:,kk) ;
				end
					
				[X(:, kk), vv] =  mgs(X(:,1:kk-1), tmp(1:n,1)) ;
                X(:, kk) = X(:, kk) / vv(end) ;
				if tdefl == 0
					tmp = [ V(:,kk) ; SHt(kk,kk)*V(:,kk) ] ;
					tmp = mgs_skew( Qgent, Qsparse, tmp, 1, 1 ) ;
                else
					tmp = V(:,kk) ;
				end
				[V(:, kk), vv] =  mgs(V(:,1:kk-1), tmp(1:n,1)) ;
                V(:, kk) = V(:, kk) / vv(end) ;
            end
            if b_f_ex == 0
                G = V' * E * X ;
            else
                G = V' * f_ex( 'N', X ) ;
            end
            if b_f_ax == 0
                AX = A*X ;
            else
                AX = f_ax( 'N', X ) ;
            end
            H = V' * (AX) ;

            if b_f_mx == 0
                MX = M*X ;
            else
                MX = f_mx( 'N', X ) ;
            end
            K = V' * (MX) ;

            nrestart = nrestart + 1 ;
            if displ
                fprintf( 'Restart %d, %d eigenvals found.\n', nrestart, nr_converged ) ;
            end
        end %end if found
    end

    if nr_converged == nwanted || nrestart == maxrestarts %nr_converged == nshifts
        rightev = Q ;
        leftev = Qt ;
        poles = R.' ;
        
        for i = 1 : length( poles )
            residues(i,1) = (rightev(:,i).' * c) * (leftev(:,i)' * b) ;
        end
        %sort on largest residual
        if strategy == 'LR'
            [y,idx] = sort(-abs(residues) ./ abs(poles)) ;
        else
            [y,idx] = sort(-abs(residues) ) ;
        end
        residues = residues(idx) ;
        poles = poles(idx) ;
        rightev = rightev(:, idx) ;
        leftev = leftev(:, idx) ;
        break ;
    end
end %end while

if displ
    output = sprintf('Number of LU factorizations: %i', nr_solves ) ;
    disp( output ) ;
    fprintf( 'Dominant poles:\n' ) ;
    disp( poles ) ;
end

function [SH, SHt, UR, URt] = select_max_res( H, G, K, X, V, b, c, strat, ...
                                         b_d, c_d, tdefl, ...
                                         A, E, M, f_ex, f_mx, yEx_scaling )
% select the approximate eigenpair with largest residue

n = size(H, 1) ;

%compute and sort the eigendecomps
[Vs,D] = polyeig( H, G, K ) ;
[y,idx] = sort(D) ;
dd = D ;
dd = dd(idx) ;
Vs = Vs(:, idx) ;

[Vt, Dt] = polyeig(H', G', K') ;
[y,idx] = sort(Dt') ;
ddt = Dt ;
ddt = ddt(idx) ;
Vt = Vt(:, idx) ;

%compute the approximate residues
X = X*Vs ;
V = V*Vt ;

for i = 1 : 2*n
	V(:,i) = V(:,i) / norm(V(:,i)) ;
    X(:,i) = X(:,i) / norm(X(:,i)) ;
    
    if tdefl
        xg = [X(:,i) ; dd(i)*X(:,i)] ; xg = xg / norm(xg) ;
        vg = [V(:,i) ; dd(i)'*V(:,i)] ; vg = vg / norm(vg) ;
    end
    
    %the method appears to be very sensitive to the way the MVs are done
    if yEx_scaling
    %	X(:,i) = X(:,i) / ( (V(:,i)' * ((E + dd(i)^2*M)*X(:,i)) )) ;
        
        if isa( f_ax, 'function_handle' )
            Ax = f_ax( 'N', X(:, i) ) ;
        else
            Ax = A*X(:, i) ;
        end
        if isa( f_mx, 'function_handle' )
            Mx = f_mx( 'N', X(:, i) ) ;
        else
            Mx = M*X(:, i) ;
        end
        X(:,i) = X(:,i) / ( (V(:,i)' * (Ax + dd(i)^2*Mx))) ;
    end
    
    if ~tdefl
        residue(i,1) = abs(( X(:,i).' * c) * (dd(i)*V(:,i)' * b) ) ;
    else
        %compute residues via deflated c and b to avoid convergence to found
        residue(i,1) = abs( (c_d'*xg) * (vg'*b_d) ) ;
    end
end

switch strat
    case 'LS'  
        residue = residue ./ abs(ddt) ;
    case 'LR'  
        residue = residue ./ abs(real(ddt)) ;
    case 'HF' %s*H(s)
        residue = (residue .* abs(imag(ddt.'))) ./ abs(real(ddt.')) ;
    case 'HFN' %s*H(s)
        residue = (residue .* abs(imag(ddt.'))) ;
    case 'LF' %H(s)/s
        residue = residue ./ ( abs(imag(ddt.')) .*abs(real(ddt.')) ) ;
    case 'SI' %H(s)/s
        residue = -abs(imag(ddt.')) ;
end

[y,idx] = sort(-residue) ;
residue = residue(idx) ;
dd = dd(idx) ;
ddt = ddt(idx) ;
UR = Vs(:, idx) ;
URt = Vt(:, idx) ;

SH = diag(dd) ; 
SHt = diag(ddt) ; 

function [n, nwanted, tol, strategy, use_lu, use_lu_w_amd, kmin, kmax, maxrestarts, ...
          displ, dpa_bordered, yEx_scaling, rqitol, tdefl, ...
          f_ax, f_ex, f_mx, f_semax, f_semax_s]  = init( A, options )
% init initializes all variables needed for JD process

[n, m] = size( A ) ;
kmin = min(1, n) ; %min searchspace
kmax = min(10, n) ; %max searchspace
maxrestarts = 100 ; %max nr of restarts
nwanted = 5 ;
fields = fieldnames( options ) ;


if isfield( options, 'nwanted' )
    nwanted = options.nwanted ;
end

if isfield( options, 'kmin' )
    kmin = options.kmin ;
end

if isfield( options, 'kmax' )
    kmax = options.kmax ;
end

if isfield( options, 'maxrestarts' )
    maxrestarts = options.maxrestarts ;
end

if isfield( options, 'tol' ) %Ritz pair tolerance
    tol = options.tol ;
else
    tol = 1e-10 ;
end
if isfield( options, 'displ') %show convergence info
    displ = options.displ ;
else
    displ = 1 ;
end

if isfield( options, 'strategy' ) %target
    strategy = options.strategy ;
else
    strategy = 'LM' ;
end


if isfield( options, 'use_lu' ) %target
    use_lu = options.use_lu ;
else
    use_lu = 1 ;
end

if isfield( options, 'use_lu_w_amd' ) %target
    use_lu_w_amd = options.use_lu_w_amd ;
else
    use_lu_w_amd = 1 ;
end

if isfield( options, 'dpa_bordered' ) %target
    dpa_bordered = options.dpa_bordered ;
else
    if use_lu == 1
        dpa_bordered = 0 ;
    else
        dpa_bordered = 1 ;
    end
end
if use_lu == 0 && dpa_bordered == 0
    fprintf( 'Warning: using \\ for solves and *no* dpa_bordered\n' ) ;
    fprintf( '         might lead to unstable iterations/stagnation\n' ) ;
end

if isfield( options, 'yEx_scaling' )
    yEx_scaling = options.yEx_scaling ;
else
    yEx_scaling = 0 ;
end

if isfield( options, 'rqitol' )
    rqitol = options.rqitol ;
else
    rqitol = 1e-4 ;
end

if isfield( options, 'turbo_deflation' )
    tdefl = options.turbo_deflation ;
else
    tdefl = 1 ;
end

if isfield( options, 'f_ax' )
    f_ax = options.f_ax ;
else
    f_ax = 0 ;
end

if isfield( options, 'f_ex' )
    f_ex = options.f_ex ;
else
    f_ex = 0 ;
end

if isfield( options, 'f_mx' )
    f_mx = options.f_mx ;
else
    f_mx = 0 ;
end

if isfield( options, 'f_semax' )
    f_semax = options.f_semax ;
else
    f_semax = 0 ;
end

if isfield( options, 'f_semax_s' )
    f_semax_s = options.f_semax_s ;
else
    f_semax_s = 0 ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lambda, x, v, rqs, res, nr_its] = ...
             rqi_quad(A, B, C, vialu, use_amd, f_ax, f_ex, f_mx, ...
                      f_semax, f_semax_s, options)
% [lambda, x, rqs,res] = rqi_quad(A, B, C, options) computes an eigenvalue near options.sigma
% using Rayleigh Quotient iteration for quadratic problems
%		Ax = lambda Bx + lambda^2Cx
% If instead options.v0 is given,
% this will used as initial vector.
%
% Input:
%     A : (nxn) matrix 
%     B : (nxn) matrix 
%     C : (nxn) matrix 
%     vialu: if 1 then use LU else use \ solve
%     options: struct with
%         tol  : tolerance for Rayleigh pair residual (1e-8)
%         maxit : maximum number of iterations
%         x0,v0    : initial guess for eigenvector (ones)
%         displ: show convergence information (0)
%         sigma: initial eigenvalue guess
%                
% Output:
%     lambda,x: approximate eigenpair with ||Ax-\lambda x||_2 < tol
%     rqs : intermediate rayleigh quotients
%     res : intermediate residuals
%
%
% Joost Rommes (C) 2006 -- 2008
% rommes@gmail.com
% http://rommes.googlepages.com

if nargin < 4
	vialu = 1 ;
    options = struct([]) ;
end

if nargin < 5
    options = struct([]) ;
end
nr_its = 0 ;
[n, x0, v0, tol, maxit, sigma, displ] = init_qrqi( A, options ) ;

%I = eye(n,n) ;
x = x0 / norm(x0) ;
v = v0 / norm(v0) ;

res = [] ;
rqs = [] ;

mu = sigma ;
tgt = sigma ;
rqs(1) = mu ;

displ = 0 ;
if displ
    fprintf( 'iter %d: mu = %s, nr = %s\n', 0, num2str(mu), '-') ;
end

b_f_ax = isa( f_ax, 'function_handle' ) ;
b_f_ex = isa( f_ex, 'function_handle' ) ;
b_f_mx = isa( f_mx, 'function_handle' ) ;
b_f_semax = isa( f_semax, 'function_handle' ) ;
b_f_semax_s = isa( f_semax_s, 'function_handle' ) ;

for k = 1 : maxit
	nr_its = nr_its + 1 ;
    
	if b_f_ex == 0
        Bx = B*x ;
        Btv = B'*v ;
    else
        Bx = f_ex( 'N', x) ;
        Btv = f_ex( 'T', v) ;
    end
	if b_f_mx == 0
        Cx = C*x ;
        Ctv = C'*v ;
    else
        Cx = f_mx( 'N', x) ;
        Ctv = f_mx( 'T', v) ;
    end
    
    if b_f_semax_s == 0
        ABC = A + mu*B + mu^2*C ;
        if vialu
	    	[x,v] = solve_with_lu( ABC, Bx + 2*mu*Cx, Btv + 2*mu'*Ctv, use_amd ) ;
            %[L,U] = lu(ABC) ;     
		    %x = U \ (L \ ( (Bx + 2*mu*Cx))) ;     
		    %v = L' \ (U' \ ( (Btv + 2*mu'*Ctv))) ;
        else
		    x = ABC \ ( (Bx + 2*mu*Cx)) ; 
    
    		v = ABC' \ ( (Btv + 2*mu'*Ctv)) ;	
    	end
	else
        x = f_semax_s( 'N', (Bx + 2*mu*Cx), mu ) ;
        v = f_semax_s( 'T', (Bx + 2*mu*Cx), mu ) ;    
    end
    %mu = mu - ( vo'*x ) / ( v'*(2*mu*C+B)*x ) 
	x = x / norm(x) ;
	v = v / norm(v) ;

    mu = full( grq(A, B, C, x, v, tgt, f_ax, f_ex, f_mx) ) ;

    rqs(k+1) = mu ;
    if b_f_semax == 0
        r = A*x + mu * (B* x) + mu^2*(C*x) ;
    else
        r = f_semax( 'N', x, mu ) ;
    end
    res(k+1) = norm( r ) ;
    
	if displ
		fprintf( 'iter %d: mu = %s, nr = %s\n', k, num2str(mu), num2str(res(k+1))) ;
    end
	if res(k+1) < tol
        break ;
    end
	if strcmp( num2str(mu), 'NaN' )
		lambda = 0 ;
		fprintf( 'RQI exited with NaN' ) ;
		return ;
	end
end

lambda = mu ;
%x = v ;

function mu = grq(A, B, C, x, v, tgt, f_ax, f_ex, f_mx)
b_f_ax = isa( f_ax, 'function_handle' ) ;
b_f_ex = isa( f_ex, 'function_handle' ) ;
b_f_mx = isa( f_mx, 'function_handle' ) ;

if b_f_ax == 0
    Ax = A*x ;
else
    Ax = f_ax( 'N', x ) ;
end
if b_f_ex == 0
    Bx = B*x ;
else
    Bx = f_ex( 'N', x ) ;
end
if b_f_mx == 0
    Cx = C*x ;
else
    Cx = f_mx( 'N', x ) ;
end
c = v'*(Ax) ;
b = v'*(Bx) ;
a = v'*(Cx) ;
d = sqrt( b^2 - 4*a*c ) ;
mu1 = (-b + d ) / (2*a) ;
mu2 = (-b - d ) / (2*a) ; %, pause;

if abs( tgt -mu1 ) < abs( tgt -mu2)
	mu = mu1 ;
else
    mu = mu2 ;
end

function [n, x0, v0, tol, maxit, sigma, displ] = init_qrqi( A, options )
% init initializes all variables

[n, m] = size( A ) ;
v0 = rand( n, 1 ) ; %initial guess
x0 = rand( n, 1 ) ; %initial guess
maxit = 20 ; %max number of iterations
sigma = 'n' ;
fields = fieldnames( options ) ;

if strmatch( 'v0', fields )
    v0 = options.v0 ;
end

if strmatch( 'x0', fields )
    x0 = options.x0 ;
end

if strmatch( 'maxit', fields )
    maxit = options.maxit ;
end

if strmatch( 'sigma', fields )
    sigma = options.sigma ;
end

if strmatch( 'tol', fields ) %Ritz pair tolerance
    tol = options.tol ;
else
    tol = 1e-12 ;
end

if strmatch( 'displ', fields ) %show convergence info
    displ = options.displ ;
else
    displ = 0 ;
end

function [w, h] = mgs(V, v, iter_ref)
% [w,h] = mgs(V, v, iter_ref) orthogonalizes v to V, with
% optional iterarive refinement (Modified Gram-Schmidt). Note
% that w is not normalized (its norm is in h(length(w)+1).
%
% Input:
%     V       : (n,k) orthogonal matrix (V'*V=I)
%     v       : (n,1) vector
%     iter_ref: 0: no refinement (default)
%               1 : refinement
% Output:
%     w       : v orthogonalized to V
%     h       : (k+1,1) vector with V'*v + ||w||*e_k+1
% Joost Rommes (C) 2004
% rommes@math.uu.nl

if nargin < 3
    iter_ref = 0 ;
end

[n,k] = size(V) ;

if k < 1
    w = v ;
    h = [norm(w)];
    return ;
end

w = v ;
nv = norm(v) ;
kappa = 0.1 ;

for i = 1 : k
    h(i) = V(:,i)' * w ;
    w    = w - h(i)*V(:,i) ;
end

if iter_ref & norm(w) < kappa * nv
    for i = 1 : k
        p    = V(:,i)' * w ;
        h(i) = h(i) + p ;
        w    = w - p*V(:,i) ;
    end
end

h(k+1) = norm(w) ;
h = h.' ;

if k > 0 && abs( V(:,k)'*w ) > 0.1
    fprintf( 'Warning: MGS stability problem: V^T*w = %s\n', num2str(V(:,k)'*w) ) ;
end

function [u, hu] = mgs_skew(Q, Qt, u, iter_ref, ogs)
% [u, hu] = mgs_skew(Q, Qt, u, iter_ref) computes the skew projection
% PI_i (I - q_iqt_i')u with optional iterative refinement
%
% Input:
%     Q, Qt       : (n,k) orthogonal matrices with Qt'Q = I
%     u       : (n,1) vectors
%     iter_ref: 0: no refinement (default)
%               1 : refinement
%     ogs     : 0: mgs (default)
%               1: ogs
% Output:
%     u       : skew projected u ;
%     hu       : (k+1,1) vector with hu(end) = ||u||_2
%
% Joost Rommes (C) 2005
% rommes@math.uu.nl

if nargin < 5
    iter_ref = 0 ;
    ogs = 0 ;
end

if nargin > 6
    ogs = 0 ;
end
error = 0 ;

[n,k] = size(Q) ;
[nu,nk] = size(Qt) ;
if n ~= nu | k ~= nk
    fprintf( 'mgs_skew: unequal basis dimensions\n' ) ;
    return ;
end

if k < 1
    hu = [norm(u)];
    return ;
end

nu = norm(u) ;
kappa = 0.1 ;
hu = [] ;

if ogs
    u = u - Q*(Qt'*u) ; %norm(Qt'*u)
    if 1%norm(u) < kappa * nu 
        u = u - Q*(Qt'*u) ; 
    end
    %nqqu = norm( Qt'*u) ; %nqqt = norm(Qt'*Q)
    if norm(u) < 1e-8 %|| nqqu > 1e-14
        fprintf( 'New random vector in MGS_skew\n' ) ;
        u = rand(n,1) ;
        u = u - Q*(Qt'*u) ;
        u = u - Q*(Qt'*u) ;
    end
    
    return ;
end

for i = 1 : k
    hu(i) = Qt(:,i)' * u  ;
    u    = u - hu(i)*Q(:,i) ;
end

if iter_ref && (norm(u) < kappa * nu )
    for i = 1 : k
        p    = Qt(:,i)' * u ;
        hu(i) = hu(i) + p ;
        u    = u - p*Q(:,i) ;
    end
end

hu(k+1) = norm(u) ;
hu = hu.' ;

if k > 0 && abs( Qt(:,k)'*u ) > 0.1
    fprintf( 'Warning: Skew-MGS stability problem: V^T*u = %s\n', num2str(Qt(:,k)'*u) ) ;
end

function [x,y] = solve_with_lu( sEmA, b, c, use_amd )
%solves sEmA*x=b and sEmA'*y = c
        if use_amd
            [L,U,p,q] = lu(sEmA,'vector') ;
            %b = S\b ;
            x = L\b(p,:) ;
            x(q,:) = U\x ;
            
            if nargout > 1
                y = U'\c(q,:) ;
                y(p,:) = (L'\y) ;
                %y = S' \ y ;
            end
        else
            [L,U] = lu(sEmA) ;
            x = U \ ( L \ b ) ;
            if nargout > 1
                y = L' \ ( U' \ c ) ;
            end
        end    
