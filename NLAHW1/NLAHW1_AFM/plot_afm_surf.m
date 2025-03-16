function plot_afm_surf(F,surf)


% function that plots (filled in) top of cantilever with displacement
% as 1D representation of cantilever

% INPUT: F= nodes+dU with dU displacement (can be zero)
%           FORMAT= 3xtotalDOFs
%        surf = given surface


l = ones(3*size(F,1),1);
  
hold on

 ord=[1,2,3,4,1];


for e=1:size(surf,1)
  
   x=F(surf(e,ord),1);
   y=F(surf(e,ord),2);      
   z=F(surf(e,ord),3);
   
   f=l(surf(e,ord));
   
   fill3(x,y,z,f)
end

shading interp
axis equal
hold off
end