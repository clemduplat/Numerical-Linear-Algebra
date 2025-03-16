function FRF = bode_from_function(fun,b,c,s)
    n = length(s);
    FRF = zeros(n,size(c,2));
    parfor j = 1:n
        FRF(j,:) = c'*(fun(s(j))\b);
    end
end