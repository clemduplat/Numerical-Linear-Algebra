function kth_eigenvalue = findKthEigenvalue(T, k, tol)

    %----------Gershorgin Theorem-----------%
    % First we will initialize the search interval [a,b] using Gershorgin
    % Theorem, each diagonal elements define a G circle, union of all will
    % contain all eigenvalues
    %We  know that each eigenvalue of A lies in at least one of the G.
    %disc D(a_ii,R_i) where R_i=sum(a_ij) for j dirrent from i
    % So the interval is [a_ii-R_i;a_ii+R_i]
    % Suppose now that a= lower bound and b=upper bound
    n=size(T,1);
    
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
    
    %----Implementation of Classic Bisection method-----%
    while (b-a)>tol
        m= (a+b)/2;
        if countEigenvaluesBelow(T,m)<k
            a=m;
        else
            b=m;
        end
    end
    
    kth_eigenvalue=(a+b)/2;
end

