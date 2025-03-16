function [Ahat, Ehat, bhat, chat, V] = irka(A, E, b, c, sigma0, tol)
    err = Inf;
    d = numel(sigma0);
    n = size(A,1);
    V = zeros(n,d);
    W = zeros(n,d);
    
    sigma = sigma0;

    while err > tol
        disp(strcat('Err: ', num2str(err)))
        for i=1:d
            Vi = (sigma(i)*E - A) \ b;
            Wi = (sigma(i)*E - A)' \ c;
            
            Vi = Vi - V(:,1:i-1)*(V(:,1:i-1)'*Vi);
            Vi = Vi - V(:,1:i-1)*(V(:,1:i-1)'*Vi);
            Wi = Wi - W(:,1:i-1)*(W(:,1:i-1)'*Wi);
            Wi = Wi - W(:,1:i-1)*(W(:,1:i-1)'*Wi);

            V(:,i) = Vi / norm(Vi);
            W(:,i) = Wi / norm(Wi);
        end
        Ehat = W'*E*V;
        Ahat = W'*A*V;        
        lambda = sort(eig(Ahat, Ehat));
        diffs = zeros(d,1);
        for i=1:d
            diffs(i) = min(abs(lambda(i) + sigma));
        end
        
        err = norm(diffs) / norm(lambda);
        sigma = -lambda;
    end
    
    bhat = W'*b;
    chat = V'*c;
end
