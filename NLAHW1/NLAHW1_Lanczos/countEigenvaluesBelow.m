function numEigenvalues = countEigenvaluesBelow(T, x)
%---------This is a function counting the number of eigenvalues below x in T using sturm sequence property----%
    n= size(T, 1); 
    p0=1; %always 1 b
    %------calculate d the determinant of the 1*1 submatrix on the left(top)---%
    p1=T(1,1)- x;  
    numEigenvalues=0; 
    
    if p1< 0 || p1==0
        % first sign change
        numEigenvalues = numEigenvalues + 1;
   
    end
    
    for i= 2:n
        %---------p_i(x) formula for characteristic polynomial--------%
        pi=(T(i,i) - x)*p1 - T(i-1,i)^2 *p0; 
        
        if pi<0
            if p1>0
                numEigenvalues=numEigenvalues+1; % Sign changed
            end
        elseif pi>0
            if p1<0
                numEigenvalues=numEigenvalues+1; % Sign changed
            end
        elseif(pi==0)
                numEigenvalues=numEigenvalues+1; % Sign changed
        end
        % -------Update polynomials for the next iteration-----%
        p0=p1;
        p1=pi;
    end
    
end
