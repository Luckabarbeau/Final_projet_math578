function LM_stancil = get_lagrange_multiplier_stancil_weak_avreage(TR,cell,X,size_A)
%GET_LAGRANGE_MULTIPLIER_STANCIL Summary of this function goes here
% return the stancil for the lagrange multiplier at point X

dx=X(1)-TR.points(TR.connect(cell,1),1);
dy=X(2)-TR.points(TR.connect(cell,1),2);
hx=TR.points(TR.connect(cell,2),1)-TR.points(TR.connect(cell,1),1);
hy=TR.points(TR.connect(cell,4),2)-TR.points(TR.connect(cell,1),2);

% k0=TR.connect(cell,1);
% k1=TR.connect(cell,2);
% k2=TR.connect(cell,3);
% k3=TR.connect(cell,4);


LM_stancil=zeros(1,size_A);
%u0
LM_stancil(TR.connect(cell,1))=1/4;

%u1
LM_stancil(TR.connect(cell,2))=1/4;

%u2
LM_stancil(TR.connect(cell,3))=1/4;

%u3
LM_stancil(TR.connect(cell,4))=1/4;


end