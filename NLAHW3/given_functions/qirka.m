function [Mr, Dr, Kr, Br, Cr, sigma, Vr] = qirka(M, D, K, B, C, sigma0, tol)
    err = Inf;
    d = numel(sigma0);
    n = size(M,1);
    Vr = zeros(n,d);
    Wr = zeros(n,d);

    sigma = sigma0;

    while err > tol
        disp(['Err: ', num2str(err)])
        for i=1:d
            % taking into account the second-order system
            Vi = (sigma(i)^2*M + sigma(i)*D + K) \ B;
            Wi = (sigma(i)^2*M + sigma(i)*D + K)' \ C;

      
            Vi = Vi - Vr(:,1:i-1)*(Vr(:,1:i-1)'*Vi);
            Vi = Vi - Vr(:,1:i-1)*(Vr(:,1:i-1)'*Vi);
            Wi = Wi - Wr(:,1:i-1)*(Wr(:,1:i-1)'*Wi);
            Wi = Wi - Wr(:,1:i-1)*(Wr(:,1:i-1)'*Wi);

            Vr(:,i) = Vi / norm(Vi);
            Wr(:,i) = Wi / norm(Wi);
        end

       
        Mr = Wr'*M*Vr;
        Dr = Wr'*D*Vr;
        Kr = Wr'*K*Vr;

      
        lambda = sort(polyeig(Kr, Dr, Mr));

        
        diffs = zeros(d,1);
        for i=1:d
            diffs(i) = min(abs(lambda(i) + sigma));
        end
        
        err = norm(diffs) / norm(lambda);
        sigma = -lambda;
    end

    
    Br = Wr'*B;
    Cr = Vr'*C;
end
