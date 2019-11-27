function LM_stancil = get_lagrange_multiplier_stancil(TR,cell,X,size_A)
%GET_LAGRANGE_MULTIPLIER_STANCIL Summary of this function goes here
% return the stancil for the lagrange multiplier at point X

dx=X(1)-TR.points(TR.connect(cell,1),1);
dy=X(2)-TR.points(TR.connect(cell,1),2);
hx=TR.points(TR.connect(cell,2),1)-TR.points(TR.connect(cell,1),1);
hy=TR.points(TR.connect(cell,4),2)-TR.points(TR.connect(cell,1),2);

LM_stancil=zeros(1,size_A);
%u0
LM_stancil(TR.connect(cell,1))=(dy/hy-1)*dx/hx-(dy/hy)+1;

%u1
LM_stancil(TR.connect(cell,2))=(-dy/hy+1)*dx/hx;

%u2
LM_stancil(TR.connect(cell,3))=dy/hy*dx/hx;

%u3
LM_stancil(TR.connect(cell,4))=dy/hy*(-dx/hx+1);


end

