function LM_stancil = get_lagrange_multiplier_stancil_weak_grad(TR,cell,X,size_A)
%GET_LAGRANGE_MULTIPLIER_STANCIL Summary of this function goes here
% return the stancil for the lagrange multiplier at point X

dx=X(1)-TR.points(TR.connect(cell,1),1);
dy=X(2)-TR.points(TR.connect(cell,1),2);
hx=TR.points(TR.connect(cell,2),1)-TR.points(TR.connect(cell,1),1);
hy=TR.points(TR.connect(cell,4),2)-TR.points(TR.connect(cell,1),2);

k0=TR.connect(cell,1);
k1=TR.connect(cell,2);
k2=TR.connect(cell,3);
k3=TR.connect(cell,4);


alpha=0.001;
LM_stancil=zeros(1,size_A);
% %u0
% LM_stancil=LM_stancil+1/alpha*(get_stancil_dx(TR,k0,1)+get_stancil_dy(TR,k0,1))*((dy/hy-1)*dx/hx-(dy/hy)+1);
% 
% %u1
% LM_stancil=LM_stancil+1/alpha*(get_stancil_dx(TR,k1,2)+get_stancil_dy(TR,k1,2))*((-dy/hy+1)*dx/hx);
% 
% %u2
% LM_stancil=LM_stancil+1/alpha*(get_stancil_dx(TR,k2,3)+get_stancil_dy(TR,k2,3))*dy/hy*dx/hx;
% 
% %u3
% LM_stancil=LM_stancil+1/alpha*(get_stancil_dx(TR,k3,4)+get_stancil_dy(TR,k3,4))*dy/hy*(-dx/hx+1);

%u0
LM_stancil=LM_stancil+1/alpha*(-get_stancil_L(TR,k0,1))*((dy/hy-1)*dx/hx-(dy/hy)+1);

%u1
LM_stancil=LM_stancil+1/alpha*(-get_stancil_L(TR,k1,1))*((-dy/hy+1)*dx/hx);

%u2
LM_stancil=LM_stancil+1/alpha*(-get_stancil_L(TR,k2,1))*dy/hy*dx/hx;

%u3
LM_stancil=LM_stancil+1/alpha*(-get_stancil_L(TR,k3,1))*dy/hy*(-dx/hx+1);


end