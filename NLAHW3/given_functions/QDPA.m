function [lambda, x, y] = QDPA(M, C, K, b, c, s0, epsilon)

    sk = s0;
    k = 0;
    
    while true
    
        vk = (sk^2 * M + sk * C + K) \ b;
        wk = (sk^2 * M + sk * C + K) \ c;


        numerator = c' * vk;
        denominator = wk' * (2 * sk * M + C) * vk;
        skp1 = sk - numerator / denominator;
        
    
        x = vk / norm(vk, 2);
        y = wk / norm(wk, 2);

      
        if max(norm((sk^2 * M + sk * C + K) * x, 2), norm((skp1^2 * M + skp1 * C + K) * y, 2)) < epsilon
            lambda = skp1;
            break;
        end
        
   
        sk = skp1;
        k = k + 1;
    end
end
