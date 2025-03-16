%% System x' = Ax+Bu, y = Cx + Du
% computes the frequency response of a system with transfer function
% H(s) = c'*(sE - A)^(-1)*b
function FRF = bode_from_system(A,E,b,c,s)
    n = length(s);
    FRF = zeros(n,size(c,2));
    parfor j = 1:n
        FRF(j,:) = c'*((s(j)*E - A)\b);
    end
end