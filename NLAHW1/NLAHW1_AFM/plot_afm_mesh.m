function plot_afm_mesh(nodes,elements)
%Plot AFM mesh

%Input: nodes at rest, in format [3xn]


% Initialize
nel = length(elements) ;                  % number of elements
nnel = size(elements,2);                % number of nodes per element

X = zeros(nnel,nel) ;
Y = zeros(nnel,nel) ;
Z = zeros(nnel,nel) ;


c = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];
XYZ = cell(1,nel) ;
for e=1:nel
    nd=elements(e,:);
    X(:,e) = nodes(nd,1) ;
    Y(:,e) = nodes(nd,2) ;
    Z(:,e) = nodes(nd,3) ;
    XYZ{e} = [X(:,e)  Y(:,e) Z(:,e)] ;
end
% Plot FEM mesh 
figure
set(gcf,'color','w')
axis off 
axis equal
cellfun(@patch,repmat({'Vertices'},1,nel),XYZ,...
 repmat({'Faces'},1,nel),repmat({c},1,nel),repmat({'Facecolor'},1,nel),repmat({[0 1 1]},1,nel));
view(3)
end