function LM_stancil = get_lagrange_multiplier_stancil_weak_sum(TR,cell,X_all,cells_sum)
%GET_LAGRANGE_MULTIPLIER_STANCIL Summary of this function goes here
% return the stancil for the lagrange multiplier at point X

LM_stancil=zeros(1,length(TR.points));
line=find(cells_sum.index==cell);
if isempty(line)
    error('the cell does not containt immersed boundary point')
end


constant_factor=1;
epsilone=0.00001;
order=02;


for i=1:length(cells_sum.boundary{line})
X=X_all(cells_sum.boundary{line}(i),:);
dx=X(1)-TR.points(TR.connect(cell,1),1);
dy=X(2)-TR.points(TR.connect(cell,1),2);
hx=TR.points(TR.connect(cell,2),1)-TR.points(TR.connect(cell,1),1);
hy=TR.points(TR.connect(cell,4),2)-TR.points(TR.connect(cell,1),2);

%u0
LM_stancil(TR.connect(cell,1))=LM_stancil(TR.connect(cell,1))+(dy/hy-1)*dx/hx-(dy/hy)+1;

%u1
LM_stancil(TR.connect(cell,2))=LM_stancil(TR.connect(cell,2))+(-dy/hy+1)*dx/hx;

%u2
LM_stancil(TR.connect(cell,3))=LM_stancil(TR.connect(cell,3))+dy/hy*dx/hx;

%u3
LM_stancil(TR.connect(cell,4))=LM_stancil(TR.connect(cell,4))+dy/hy*(-dx/hx+1);

% 
%  l1=((dx/hx)^2+(dy/hy)^2)^(order/2)+epsilone;
%  l2=((1-dx/hx)^2+(dy/hy)^2)^(order/2)+epsilone;
%  l3=((1-dx/hx)^2+(1-dy/hy)^2)^(order/2)+epsilone;
%  l4=((dx/hx)^2+(1-dy/hy)^2)^(order/2)+epsilone;
% 
% totalnorm=1/l1+1/l2+1/l3+1/l4;
% 
% 
% 
% 
%     %u0
% LM_stancil(TR.connect(cell,1))=LM_stancil(TR.connect(cell,1))+constant_factor/4+(1-constant_factor)*1/l1/totalnorm;
% 
% %u1
% LM_stancil(TR.connect(cell,2))=LM_stancil(TR.connect(cell,2))+constant_factor/4+(1-constant_factor)*1/l2/totalnorm;
% 
% %u2
% LM_stancil(TR.connect(cell,3))=LM_stancil(TR.connect(cell,3))+constant_factor/4+(1-constant_factor)*1/l3/totalnorm;
% 
% %u3
% LM_stancil(TR.connect(cell,4))=LM_stancil(TR.connect(cell,4))+constant_factor/4+(1-constant_factor)*1/l4/totalnorm;
% 
% 
% 


end


end