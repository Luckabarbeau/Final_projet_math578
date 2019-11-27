%code by Lucka Barbeau
function F_X = eval_solution_on_mesh(TR,X,Solution)

%GET_LAGRANGE_MULTIPLIER_STANCIL Summary of this function goes here
% return the stancil for the lagrange multiplier at point X


cell=which_cell(TR,X);

dx=X(1)-TR.points(TR.connect(cell,1),1);
dy=X(2)-TR.points(TR.connect(cell,1),2);
hx=TR.points(TR.connect(cell,2),1)-TR.points(TR.connect(cell,1),1);
hy=TR.points(TR.connect(cell,4),2)-TR.points(TR.connect(cell,1),2);

epsilone=0.00001;

order=0;
 l1=((dx/hx)^2+(dy/hy)^2)^(order/2)+epsilone;
 l2=((1-dx/hx)^2+(dy/hy)^2)^(order/2)+epsilone;
 l3=((1-dx/hx)^2+(1-dy/hy)^2)^(order/2)+epsilone;
 l4=((dx/hx)^2+(1-dy/hy)^2)^(order/2)+epsilone;


totalnorm=1/l1+1/l2+1/l3+1/l4;

LM_stancil=zeros(1,length(Solution));
constant_factor=0;

    %u0
LM_stancil(TR.connect(cell,1))=constant_factor/4+(1-constant_factor)*1/l1/totalnorm;

%u1
LM_stancil(TR.connect(cell,2))=constant_factor/4+(1-constant_factor)*1/l2/totalnorm;

%u2
LM_stancil(TR.connect(cell,3))=constant_factor/4+(1-constant_factor)*1/l3/totalnorm;

%u3
LM_stancil(TR.connect(cell,4))=constant_factor/4+(1-constant_factor)*1/l4/totalnorm;

F_X=LM_stancil*Solution;
% F_X=((dy/hy-1)*dx/hx-(dy/hy)+1)*Solution((TR.connect(cell,1)))+((-dy/hy+1)*dx/hx)*Solution((TR.connect(cell,2)))+(dy/hy*dx/hx)*Solution((TR.connect(cell,3)))+dy/hy*(-dx/hx+1)*Solution((TR.connect(cell,3)));

end

