function plot_afm_frame(F,surf,linespec)

% function that plots (frame of) top of cantilever with displacement
% as 2D representation of cantilever

% INPUT: F= nodes+dU with dU displacement (can be zero)
%           FORMAT= 3xtotalDOFs
%           surf = given surface

if ( nargin < 3 )
   linespec='k-';
end
hold on
if (size(F,2) < 3)
   for c=size(F,2)+1:3
      F(:,c)=[zeros(size(F,1),1)];
   end
end
ord=[1,2,3,4,1];
x = zeros(size(ord,2));
y = zeros(size(ord,2));
z = zeros(size(ord,2));
for e=1:size(surf,1)
   for n=1:size(ord,2)
      x(n)=F(surf(e,ord(n)),1);
      y(n)=F(surf(e,ord(n)),2);      
      z(n)=F(surf(e,ord(n)),3);
   end
   plot3(x,y,z,linespec)
end
rotate3d on
axis equal
hold off
end