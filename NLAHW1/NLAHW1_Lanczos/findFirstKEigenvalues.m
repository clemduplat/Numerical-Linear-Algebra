function k_eigenvalues = findFirstKEigenvalues(T, k, tol)
    n = size(T, 1);
    %-------This time stock th k eigenvalues in an array----%
    k_eigenvalues=zeros(k, 1);  
        
    %----------Gershorgin Theorem-----------%
    % Search for the biggest disk to start
    % For the first row we only have one off-diagonal element
    a=T(1,1)-abs(T(2,1)); 
    b=T(1,1)+ abs(T(2,1)); 
    
    % for the rows before te last one, they have 2 off-diagonals elements
    for i= 2:n-1
        R_i=abs(T(i, i-1)) +abs(T(i, i+1));
        a=min(a, T(i, i)-R_i);
        b=max(b, T(i, i)+R_i);
    end
    % last row as the first one only contain one off-diagonal element
    a= min(a,T(n, n)-abs(T(n, n-1))); 
    b= max(b,T(n, n)+abs(T(n, n-1)));
   
    %----Implementation of Bisection method-----%
    for i=1:k
        low=a;
        high=b;
        count=0;
        while (high-low)>tol
            mid=(low+high)/2;
            if countEigenvaluesBelow(T, mid)<i
                low=mid;
            else
                high=mid;
            end
            if(count>10000)
                disp("While loop is too long");
                break;
            end
            count=count+1;
        end
        %---- stock the value of the i-th eigenvalue----%
        k_eigenvalues(i)=(low+high)/2;
        
        % --- Upgarde search interval-----%
        a=k_eigenvalues(i)+tol;
    end
end
