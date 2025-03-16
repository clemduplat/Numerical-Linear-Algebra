function [Ahat, Ehat, bhat, chat, i] = grka(A, E, b, c, sigma1, smin, smax, scount, tol)
    V = zeros(size(A,1),0);
    W = zeros(size(A,1),0);
    Ehat = zeros(0);
    Ahat = zeros(0);
    bhat = zeros(0,1);
    
    ss = linspace(smin,smax,scount);
    i = 1;
    while true
        if i == 1
            sigma = sigma1;
        elseif mod(i,2)
            sigma = -sigma;
        else
            AV = A*V(:,1:i-1);
            EV = E*V(:,1:i-1);
            ress = arrayfun(@(s) res(s,AV,EV,b,Ahat(1:i-1,1:i-1),Ehat(1:i-1,1:i-1),bhat(1:i-1)),ss);
            [mx, imax] = max(ress);
            if mx / norm(bhat(1:i-1)) < tol
                chat = V'*c;
                i = i - 1;
                return
            end
            sigma = ss(imax);
        end
        Vi = (sigma*E - A) \ b;
        Wi = (sigma*E - A)' \ c;
        
        Vi = Vi - V(:,1:i-1)*(V(:,1:i-1)'*Vi);
        Vi = Vi - V(:,1:i-1)*(V(:,1:i-1)'*Vi);
        Wi = Wi - W(:,1:i-1)*(W(:,1:i-1)'*Wi);
        Wi = Wi - W(:,1:i-1)*(W(:,1:i-1)'*Wi);
        V(:,i) = Vi / norm(Vi); W(:,i) = Wi / norm(Wi);
        
        Ehat(i,1:i) = W(:,i)'*E*V(:,1:i);
        Ehat(1:i-1,i) = W(:,1:i-1)'*E*V(:,i);
        Ahat(i,1:i) = W(:,i)'*A*V(:,1:i);
        Ahat(1:i-1,i) = W(:,1:i-1)'*A*V(:,i);
        bhat(i,1) = W(:,i)'*b;
        i = i + 1;
    end
end

function res = res(s,AV,EV,b,Ahat,Ehat,bhat)
    xhat = (s*Ehat - Ahat) \ bhat;
    res = norm((s*EV - AV)*xhat - b);
end
