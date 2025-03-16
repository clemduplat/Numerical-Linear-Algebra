function [lambda_hat, x, y] = DPA(E, A, b, c, s0, tol)
    k = 0;
    err = inf;
    sk = s0;
    
    while err > tol
        %solve
        vk = (sk*E - A) \ b; 
        %solve
        wk = (sk*E - A)' \ c; 
        %new pole estimate
        lambda_new = sk - (c'*vk) / (wk'*E*vk); 
        x = vk / norm(vk); 
        y = wk / norm(wk); 
        
        err = norm(A*x - lambda_new*E*x); 
        sk = lambda_new; 
        k = k + 1; 
    end
    
    lambda_hat = sk; 
end
