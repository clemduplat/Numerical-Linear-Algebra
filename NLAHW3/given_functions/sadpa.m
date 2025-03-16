function [poles, residues, rightev, leftev, nr_solves, ress] = ...
           sadpa(A, E, b, c, d, s0, options)
% [poles, residues, leftev, rightev, nr_solves] = ...
%     sadpa(A, E, b, c, d, s0, options)
% computes the dominant poles of SISO transfer function
%               H(s) = c'*(sE-A)^{-1}b + d
% of the first-order descriptor (DAE) system
%               E * x' = A*x    + b*u
%                   y  = c' * x + d*u
% using the Subspace Accelerated Dominant Pole Algorithm of Joost Rommes and
% Nelson Martins. This algorithm is specifically designed for large sparse 
% matrices A and E that allow *cheap* solves  using lu or \. 
%
% When using this software, please cite
%    Joost Rommes and Nelson Martins: 
%    Efficient computation of transfer function dominant poles 
%    using subspace acceleration.
%    IEEE Trans. on Power Systems 21 (3), pp. 1218-1226, 2006
%
%    and for several improvements and advanced options
%
%    Joost Rommes: 
%    Methods for eigenvalue problems with applications in model order
%    reduction. (Chapter 3)
%    PhD thesis, Utrecht University, 2007.
%
%
% Input:
%      A: (nxn) matrix -- state matrix
%      b: (nx1) matrix -- input vector
%      c: (nx1) matrix -- output vector
%      E: (nxn) matrix -- descriptor matrix
%      d: (1x1) direct term of system
%      s0: (kx1) initial estimates for pole
%          if size(s0,1) > 1, s0(i) will be injected after the
%              the i-th pole has been found (ie, one solve with
%              (s0(i)*E-A) in the next iteration)
%      options: structure with
%         nwanted : number of wanted poles (5)
%         tol  : tolerance for eigenpair residual (1e-10)
%         displ: show convergence information (1)
%         strategy: select pole p and residue with maximal ('LM'):
%                          one of
%                           'LM' : largest |res| magnitude (default)
%                           'LS' : largest |res|/|p| magnitude
%                           'LR' : largest |res|/|re(p)| magnitude
%         kmin: minimal searchspace dimension (min(1,n))       
%         kmax: maximal searchspace dimension (min(10,n))
%         maxrestarts: maximal nr of restarts (100)
%
%         Function handles to user functions for MVs and solves
%         can be used instead of A and E:
%         f_ax    :   Y = f_ax( tr, X) should return
%                         A*X if tr = 'N'
%                         A'*X if tr = 'T'
%         f_ex    :   Y = f_ex( tr, X) should return
%                         E*X if tr = 'N'
%                         E'*X if tr = 'T'
%         f_semax :   Y = f_semax( tr, X, s) should return
%                         (sE - A)*X  if tr = 'N'
%                         (sE - A)'*X if tr = 'T'
%         f_semax_s : Y = f_semax_s( tr, X, s) should return
%                         inv(sE - A) * X  if tr = 'N'
%                         inv((sE - A)') * X if tr = 'T'
%         Note that when using function handles, still dummy A and E need
%         to be supplied as arguments to sadpa_adi, typically just A=E=sparse(n,n).         
%
%         Advanced options:
%         use_lu: if use_lu==1, use LU for solves (default), otherwise \ (0)
%         use_lu_w_amd : if 1 use complete UMFPACK reordering and scaling
%         dpa_bordered: if dpa_bordered==1, use bordered matrices in solves
%                           ie [sE - A, -b ; c', d] (default with use_lu=0)
%                       else use just sE - A (default with use_lu=1) 
%         yEx_scaling: if yEx_scaling == 1, scale approximate eigenvectors
%                         during selection such that y'Ex = 1
%                      else scale ||x|| = ||y|| = 1 (default)
%         rqitol: tolerance on ||r|| = ||Ax-lambda Ex||. If ||r||<rqitol,
%                 refine using twosided Rayleigh quotient iteration (1e-4). Note
%                 rqitol > tol
%         turbo_deflation: if turbo_deflation == 1, use b and c to efficiently deflate
%                             (default, see Chapter 3 of PhD Thesis Joost Rommes)
%                          else use classical deflation via explicit
%                               projections with found eigenvectors
%
%
% Output:
%      poles: converged poles
%      residues: residues needed to construct bode plot
%      leftev: left eigenvectors
%      rightev: right eigenvectors
%      nr_solves: number of LU factorizations used
%      ress : residues during process
%
%
% Optimal SADPA settings may vary per transfer function. 
% If convergence is not as expected/desired, consider tuning 
% options in this order:
%        [ 0) increase options.nwanted ]
%          1) change   options.strategy, for instance to 'LR' or 'LS' in case of
%             convergence to Inf (typically when E*b=0 and/or E*c=0)
%          2) inject (complex) shifts: s0 = [w1, w2, ...]
%          3) increase options.kmin and/or options.kmax
%
%
% Joost Rommes (C) 2005 -- 2010
% rommes@gmail.com
% http://sites.google.com/site/rommes

if nargin < 7
    options = struct([]) ;
end

[n, nwanted, tol, strategy, use_lu, use_lu_w_amd, kmin, krestart, maxrestarts, ...
    displ, dpa_bordered, yEx_scaling, rqitol, tdefl, f_ax, f_ex, f_semax, f_semax_s] = ...
     init( A, options ) ;

if isa( f_semax_s, 'function_handle' ) && dpa_bordered
    error( 'Bordered DPA does not support function handles for solves' ) ;
end

b_f_ax = isa( f_ax, 'function_handle' ) ;
b_f_ex = isa( f_ex, 'function_handle' ) ;
b_f_semax = isa( f_semax, 'function_handle' ) ;
b_f_semax_s = isa( f_semax_s, 'function_handle' ) ;


conjtol = 1e-8 ; %tolerance used for complex-conjugate convergence check
imagtol = 1e-8 ; %tolerance used for complex_or_real check

nr_shifts = length(s0) ;
shifts = s0 ; %pole approximations
shift_nr = 2 ;
just_found = 0 ;

k = 0 ;
rhs = sparse([ zeros(n, 1); 1 ] );

X = zeros(n,1) ;
V = zeros(n,1) ;

nr_converged = 0 ;
poles = [] ;
leftev = [] ;
rightev = [] ;
residues = [] ;
nr_solves = 0 ;
ress = [] ;
res_count = 0 ;
iter = 0 ;

Q = zeros(n,0) ;
Qt= zeros(n,0) ;
Qsparse = zeros(n,0) ;
Qtsparse = zeros(n,0) ;
R = [] ;
Rt= [] ;
G = [] ;
H = [] ;
AX = [] ;
st = shifts(1) ;

nrestart = 0 ;

b_orig = b ;
c_orig = c ;


while nrestart < maxrestarts
    k = k+1 ;
	iter = iter + 1 ;

    if just_found && shift_nr < nr_shifts
        just_found = 0 ;
        st = shifts( shift_nr) ;
        shift_nr = shift_nr + 1 ;
    end
    
    if dpa_bordered %original dpa form
	    sEmA = [st * E - A, -b; c', d ] ;    

        if use_lu
            [xu,vu] = solve_with_lu( sEmA, b, c, use_lu_w_amd ) ;
            %[L,U] = lu( sEmA ) ;
            %xu = U \ ( L \ rhs ) ;    
            %vu = L' \ ( U' \ rhs ) ;
        else
            xu = sEmA \ rhs ;    
            vu = sEmA' \ rhs ;
        end
        x = xu(1:n) ;
        v = vu(1:n) ;
    else %dpa variant without border (exact Newton)
        if b_f_semax_s == 0
            sEmA = st * E - A ;    
    
    	    if use_lu
                [x,v] = solve_with_lu( sEmA, b, c, use_lu_w_amd ) ;
                %[L,U] = lu( sEmA ) ;
                %x = U \ ( L \ b ) ;
                %v = L' \ ( U' \ c ) ;
            else
                x = sEmA \ b ;
                v = sEmA' \ c ;
            end
        else
            x = f_semax_s( 'N', b, st ) ;
            v = f_semax_s( 'T', c, st ) ;
        end
    end
    nr_solves = nr_solves + 1 ;

    if tdefl == 0 %use traditional deflation
        x = mgs_skew( Q, Qtsparse, x, 1, 1 ) ;
        v = mgs_skew( Qt, Qsparse, v, 1, 1) ;
    end
    
	[x,h] = mgs(X, x,1) ;
    X(:, k) = x / h(end) ;

    [v,h] = mgs(V, v,1) ;
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
    
    [SH, SHt, UR, URt] = select_max_res( H, G, X, V, b, c, strategy,E, f_ex, yEx_scaling ) ;
    SG = eye(size(SH)) ;
    
    found = 1 ;
    do_rqi = 1 ;
    while found
        schurvec = X * UR(:,1) ;
        theta = SH(1,1) / SG(1,1) ;
        schurvec = schurvec / norm(schurvec) ;
        lschurvec = V * URt(:,1) ;
        lschurvec = lschurvec / norm(lschurvec) ;
        st = theta ;

        if b_f_semax == 0
            r = A*schurvec - theta * (E* schurvec) ; %TODO more efficient
        else
            r = f_semax( 'N', schurvec, theta ) ;
        end
        nr = norm(r) ;
		res_count = res_count + 1 ;
		ress( iter, res_count ) = nr;
        if displ
            fprintf('Iteration (%d): theta = %s, norm(r) = %s\n', k, num2str(theta),num2str(nr) ) ;
        end    
        %if the direct rqi residual is smaller, then use that approx.
        %this avoids unnecessary stagnation due to rounding errors        
        if nr < rqitol && do_rqi %& nr > tol            
            [schurvec, lschurvec, theta, nr, nr_solvesp] = ...
                twosided_rqi(A, E, f_ax, f_ex, f_semax, f_semax_s, ...
                             schurvec, lschurvec, theta, tol, ...
                             imagtol, use_lu, use_lu_w_amd, nr, dpa_bordered) ;
            nr_solves = nr_solves + nr_solvesp ;

            do_rqi = 0 ; %to avoid duplicate convergence if more than one found
            
            if abs(imag( theta ) )/abs(theta) < imagtol
            %possible real pole
               if b_f_semax == 0
                   rres = A*real(schurvec) - real(theta)*(E*real(schurvec)) ;
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
        end
        
        found = nr < tol && nr_converged < nwanted ;
        
        if found
            Q = [Q, schurvec] ;
            Qt = [Qt, lschurvec] ;
            if b_f_ex == 0
                Esch = E*schurvec ;
                Qsparse = [Qsparse, Esch] ;
                Qtsparse = [Qtsparse, E'*lschurvec] ;
            else
                Esch = f_ex('N', schurvec) ;
                Qsparse = [Qsparse, Esch ] ;
                Qtsparse = [Qtsparse, f_ex( 'T', lschurvec) ] ;
            end
           
            R = [R, theta] ;            
            Rt = [ Rt, conj(theta)] ;
            
            nqqt = lschurvec' * Esch ;
            Q(:,end) = Q(:,end) / nqqt ;
            Qsparse(:, end) = Qsparse(:, end) / nqqt ;

            poles = [ poles, theta ] ;

            if displ
                fprintf( 'Dominant pole %s found in iteration %d\n', num2str(theta), k ) ;
            end
            nr_converged = nr_converged + 1 ;
            just_found = 1 ;

            J = 2 : k ;
            X = X * UR(:, J) ; 
            V = V * URt(:, J) ;


            if abs(imag( theta ) )/abs(theta) < imagtol %relative
            %in case of real pole, do reortho. If a complex pair
            %is detected, do ortho after defl of complex conj
				%X = X - Q(:,end)*(Qtsparse(:,end)'*X) ; %todo not needed?
				%V = V - Qt(:,end)*(Qsparse(:,end)'*V) ;
				for kk = 1 : size(V,2)
					if tdefl == 0
                        X(:,kk) = mgs_skew( Q, Qtsparse, X(:,kk), 1, 1 ) ;
                        V(:,kk) = mgs_skew( Qt, Qsparse, V(:,kk), 1, 1 ) ;
                    end
                    [X(:, kk), vv] =  mgs(X(:,1:kk-1), X(:, kk)) ;
					X(:, kk) = X(:, kk) / vv(end) ;
					[V(:, kk), vv] =  mgs(V(:,1:kk-1), V(:, kk)) ;
					V(:, kk) = V(:, kk) / vv(end) ;
				end
            end
			%deflation from input and output vector
			if tdefl
                if b_f_ex == 0
                    b = b - E  * (Q(:,end)*(Qt(:,end)'*b)) ;
                    c = c - E' * (Qt(:,end)*(Q(:,end)'*c)) ;
                else
                    b = b - f_ex( 'N', Q(:,end)*(Qt(:,end)'*b)) ;
                    c = c - f_ex( 'T', Qt(:,end)*(Q(:,end)'*c)) ;
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
                       r = A*ccv - cce*(E*ccv) ;
                   else
                       r = f_semax( 'N', ccv, cce ) ;                   
                   end
                   if displ
                        fprintf('(conjugate) Iteration (%d): directly deflate conj...\n', k ) ;
                   end

                   if norm(r) < conjtol
                       Q = [Q, ccv] ;
                       ccet = conj(cce) ;
                       ccvt = conj( lschurvec ) ;
                       Qt = [Qt, ccvt];

                       if b_f_ex == 0
                           Esch = E*ccv ;
                           Qsparse = [Qsparse, Esch] ;
                           Qtsparse = [Qtsparse, E'*ccvt] ;
                       else
                           Esch = f_ex( 'N', ccv ) ;
                           Qsparse = [Qsparse, Esch] ;
                           Qtsparse = [Qtsparse, f_ex( 'T', ccvt) ] ;
                       end
                       R = [ R, cce ] ;
                       Rt = [ Rt, ccet ] ;

                       nqqt = ccvt' * Esch ;
                       Q(:,end) = Q(:,end) / nqqt ;
                       Qsparse(:, end) = Qsparse(:, end) / nqqt ;

                       poles = [poles, cce] ;
                       if displ
                           fprintf( '(conjugate) ...dominant pole %s found in iteration %d\n', num2str(cce), k ) ;
                       end
                       X = X - Q(:,end-1:end)*(Qtsparse(:,end-1:end)'*X) ;
                       V = V - Qt(:,end-1:end)*(Qsparse(:,end-1:end)'*V) ;

                       for kk = 1 : size(V,2)
                           if tdefl == 0
                               X(:,kk) = mgs_skew( Q, Qtsparse, X(:,kk), 1, 1 ) ;
                               V(:,kk) = mgs_skew( Qt, Qsparse, V(:,kk), 1, 1 ) ;
                           end
                           [X(:, kk), vv] =  mgs(X(:,1:kk-1), X(:, kk)) ;
                           X(:, kk) = X(:, kk) / vv(end) ;
                           
                           [V(:, kk), vv] =  mgs(V(:,1:kk-1), V(:, kk)) ;
                           V(:, kk) = V(:, kk) / vv(end) ;
                       end

                       if tdefl
                            if b_f_ex == 0
                                b = b - E  * (Q(:,end)*(Qt(:,end)'*b)) ;
                                c = c - E' * (Qt(:,end)*(Q(:,end)'*c)) ;
                            else
                                b = b - f_ex( 'N', Q(:,end)*(Qt(:,end)'*b)) ;
                                c = c - f_ex( 'T', Qt(:,end)*(Q(:,end)'*c)) ;
                            end
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

            if k > 0
                [SH, SHt, UR, URt] = select_max_res( H, G, X, V, b, c, strategy, E, f_ex, yEx_scaling ) ;
                SG = speye(size(SH)) ;
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
                    X(:,kk) = mgs_skew( Q, Qtsparse, X(:,kk), 1, 1 ) ;
                    V(:,kk) = mgs_skew( Qt, Qsparse, V(:,kk), 1, 1 ) ;
                end
                [X(:, kk), vv] =  mgs(X(:,1:kk-1), X(:, kk)) ;
                X(:, kk) = X(:, kk) / vv(end) ;
                [V(:, kk), vv] =  mgs(V(:,1:kk-1), V(:, kk)) ;
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
        b = b_orig ;
		c = c_orig ;
        for i = 1 : length( poles )
            if b_f_ex == 0
                leftev(:,i) = leftev(:, i) / ((leftev(:,i)' * E *rightev(:,i))') ;
            else
                leftev(:,i) = leftev(:, i) / ((leftev(:,i)' * f_ex( 'N', rightev(:,i)))') ;
            end
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

function [SH, SHt, UR, URt] = select_max_res( H, G, X, V, b, c, strat, E, f_ex, yEx_scaling )
% select the approximate eigenpair with largest residue

n = size(H, 1) ;

%compute and sort the eigendecomps
[Vs,D] = eig( H, G ) ;
[y,idx] = sort(diag(D)) ;
dd = diag(D) ;
dd = dd(idx) ;
%D = diag(dd) ;
Vs = Vs(:, idx) ;

[Vt, Dt] = eig(H', G') ;
[y,idx] = sort(diag(Dt')) ;
ddt = diag(Dt) ;
ddt = ddt(idx) ;
Vt = Vt(:, idx) ;

%compute the approximate residues
X = X*Vs ;
V = V*Vt ;

%b_g = V'*b ;
%c_g = X'*c ;
for i = 1 : n
	V(:,i) = V(:,i) / norm(V(:,i)) ;
    X(:,i) = X(:,i) / norm(X(:,i)) ;
    
    if yEx_scaling %can be numerically unstable
        if isa( f_ex, 'function_handle' ) == 0
            X(:,i) = X(:,i) / ( (V(:,i)' * (E*X(:,i)))) ;
        else
            X(:,i) = X(:,i) / ( (V(:,i)' * (f_ex('N',X(:,i)))) ) ;
        end
    end
    
    residue(i) = abs(( X(:,i).' * c) * (V(:,i)' * b) ) ;

%efficient
%	Vs(:,i) = Vs(:,i) / ( Vt(:,i)'*G*Vs(:,i) ) ;
%	Vs(:,i) = Vs(:,i) / ( Vt(:,i)'*Vs(:,i) ) ;
%	residue(i) = abs( (c_g'*Vs(:,i)) * (Vt(:,i)'*b_g) ) ;
end
if isnumeric( strat )
	residue = -abs( strat -dd) ;
else
	switch strat
		case 'LS'  
			residue = residue ./ abs(ddt.') ;
		case 'LR'  
			residue = residue ./ abs(real(ddt.')) ;
		case 'HF' %s*H(s)
			residue = (residue .* abs(imag(ddt.'))) ./ abs(real(ddt.')) ;
		case 'HFN' %s*H(s)
			residue = (residue .* abs(imag(ddt.'))) ;
		case 'LF' %H(s)/s
			residue = residue ./ ( abs(imag(ddt.')) .*abs(real(ddt.')) ) ;
		case 'SI' %H(s)/s
			residue = -abs(imag(ddt.')) ;
	end
end

[y,idx] = sort(-residue) ;
residue = residue(idx) ;
dd = dd(idx) ;
ddt = ddt(idx) ;
UR = Vs(:, idx) ;
URt = Vt(:, idx) ;

SH = diag(dd) ; 
SHt = diag(ddt) ; 


function [n, nwanted, tol, strategy, use_lu, use_lu_w_amd, kmin, kmax, maxrestarts, displ, dpa_bordered, yEx_scaling, rqitol, tdefl, f_ax, f_ex, f_semax, f_semax_s]  = init( A, options )
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


%%%% other functions %%%%%
function [schurvec, lschurvec, theta, nr, nr_solves] =  ...
    twosided_rqi(A, E, f_ax, f_ex, f_semax, f_semax_s, x, y, theta, tol, ...
                 imagtol, use_lu, use_lu_w_amd, r_orig, dpa_bordered)
%twosided Rayleigh quotient iteration
nr_solves = 0 ;
nrq = 1 ;
maxit = 10 ; %usually two iterations are enough
iter = 0 ;

b_f_ax = isa( f_ax, 'function_handle' ) ;
b_f_ex = isa( f_ex, 'function_handle' ) ;
b_f_semax = isa( f_semax, 'function_handle' ) ;
b_f_semax_s = isa( f_semax_s, 'function_handle' ) ;

n = size(A, 1) ;
rhs = sparse(n+1,1) ; 
rhs(n+1,1) = 1 ;

while nrq > tol && iter < maxit
    iter = iter + 1 ;

    if b_f_ex == 0
        Ex = E*x ;
        Ey = E'*y ;
    else
        Ex = f_ex( 'N', x ) ;
        Ey = f_ex( 'T', y ) ;
    end
    if b_f_semax_s == 0
        sEmA = theta*E-A ;
        if use_lu
            [x_rqi, v_rqi] = solve_with_lu( sEmA, Ex, Ey, use_lu_w_amd) ;
            %[L,U] = lu(theta*E-A) ;
            %x_rqi = U \ (L \ (Ex)) ;
            %v_rqi = L' \ (U' \ (Ey)) ;
        else
            if dpa_bordered
                sEmA = [theta * E - A, -Ex; Ey', 0 ] ;
                x_rqi = sEmA \ rhs ; x_rqi = x_rqi(1:n,1) ;
                v_rqi = sEmA' \ rhs ; v_rqi = v_rqi(1:n,1) ;
            else
                x_rqi = (sEmA \ (Ex) ) ;
                v_rqi =  (sEmA' \ (Ey) ) ;
            end
        end
    else
        x_rqi = f_semax_s( 'N', Ex, theta ) ;
        v_rqi = f_semax_s( 'T', Ey, theta ) ;
    end
    nr_solves = nr_solves + 1 ;
    x_rqi = x_rqi / norm(x_rqi) ; 
    v_rqi = v_rqi / norm(v_rqi) ; 
    if b_f_ax == 0
        Ax_rqi = A*x_rqi ;
    else
        Ax_rqi = f_ax( 'N', x_rqi ) ;
    end
    if b_f_ex == 0
        Ex_rqi = E*x_rqi ;
    else
        Ex_rqi = f_ex( 'N', x_rqi ) ;
    end
    x_rq = (v_rqi'*Ax_rqi) / (v_rqi' * Ex_rqi) ;
    
    if ~isfinite( x_rq )
        %it may happen that x_rq is Inf or Nan, if exact pole is used
        %in that case, use previous data and perturb theta
        x_rqi = x ;
        v_rqi = y ;
        x_rq  = theta + 1e-10 ;
    end
    
    v_rq = x_rq' ;
    rqi_res = Ax_rqi - x_rq * (Ex_rqi) ;
    if abs(imag(x_rq)) / abs( x_rq ) < imagtol
        %check for possible real eigenpair
        rx_rqi = real(x_rqi) ;
        
        rx_rqi = rx_rqi / norm(rx_rqi) ;
        if b_f_semax == 0
            rres = A*rx_rqi - real(x_rq)*E*rx_rqi ;
        else
            rres = f_semax( 'N', rx_rqi, real(x_rq) ) ;
        end
        nrr = norm( rres) ;
        if nrr < norm(rqi_res)
            x_rqi = rx_rqi ;
            v_rqi = real(v_rqi ) ; 
            v_rqi = v_rqi / norm(v_rqi) ;
            x_rq = real(x_rq ) ;
            rqi_res = rres ;
        end
    end    
    x = x_rqi ;
    y = v_rqi ;
    theta = x_rq ;
    nrq = norm(rqi_res) ;
    if ~isfinite( nrq )
        nrq = 1 ;
    end
end

if nrq < norm( r_orig )
%fprintf( 'RQI estimation used\n' ) ;
    schurvec = x_rqi ;
    lschurvec = v_rqi ;
    theta = x_rq ;
    r = rqi_res ;
    nr = nrq ;
else
    schurvec = x ;
    lschurvec = y ;
    r = r_orig ;
    nr = norm(r) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            x = L\b(p) ;
            x(q) = U\x ;
            
            y = U'\c(q) ;
            y(p) = (L'\y) ;
            %y = S' \ y ;
        else
            [L,U] = lu(sEmA) ;
            x = U \ ( L \ b ) ;
            y = L' \ ( U' \ c ) ;
        end    
