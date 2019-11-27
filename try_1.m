U=[1,2,3,4]
i=1
for X=0:0.01:1
    j=1
    for Y=0:0.01:1
        A=get_lagrange_multiplier_stancil([X,Y]);
        Z(i,j)=A*U';
        X_1(i,j)=X;
        Y_1(i,j)=Y;
        j=j+1
    end
    i=i+1
end



surf(X_1,Y_1,Z)


function LM_stancil = get_lagrange_multiplier_stancil(X)
%GET_LAGRANGE_MULTIPLIER_STANCIL Summary of this function goes here
% return the stancil for the lagrange multiplier at point X
% cell=which_cell(TR,X);

dx=X(1);
dy=X(2);
hx=1;
hy=1;

epsilone=0.00001;

<<<<<<< 4b7ac153599ed2555be0cbaad76fd8b6752704b3
order=5;
=======
order=16;
>>>>>>> Final_projet_stage
 l1=((dx/hx)^2+(dy/hy)^2)^(order/2)+epsilone;
 l2=((1-dx/hx)^2+(dy/hy)^2)^(order/2)+epsilone;
 l3=((1-dx/hx)^2+(1-dy/hy)^2)^(order/2)+epsilone;
 l4=((dx/hx)^2+(1-dy/hy)^2)^(order/2)+epsilone;


totalnorm=1/l1+1/l2+1/l3+1/l4;

LM_stancil=zeros(1,4);
constant_factor=0;

    %u0
LM_stancil(1)=constant_factor/4+(1-constant_factor)*1/l1/totalnorm;

%u1
LM_stancil(2)=constant_factor/4+(1-constant_factor)*1/l2/totalnorm;

%u2
LM_stancil(3)=constant_factor/4+(1-constant_factor)*1/l3/totalnorm;

%u3
LM_stancil(4)=constant_factor/4+(1-constant_factor)*1/l4/totalnorm;



end