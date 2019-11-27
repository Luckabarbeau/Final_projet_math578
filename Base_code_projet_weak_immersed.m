% code by Lucka Barbeau
% try one  heat equation Laplacien U=g(x)whit boundary condition=0

%g(x)=1 good approximation of the solution at the center is  0.073670467524337
%theorical convergence order of the method 2 for a full  grid mesh
clc ; clear all; close all

X_domain=[0 1];
Y_domain=[0 1];

% nb point 
% number on grid N by N grid

<<<<<<< 4b7ac153599ed2555be0cbaad76fd8b6752704b3
initial_Raffienement=3
nb_refinement=7

nbpoint_immersed=100;
radius=0.1;
center=[0.3 0.5];
=======
initial_Raffienement=4
initial_immersed_refinement=1
nb_refinement=5
nbpoint_immersed=20000;
radius=0.35;
center=[0.5 0.5];
X=[];
>>>>>>> Final_projet_stage
for i=1:nbpoint_immersed

        X(i,1)=radius*cos((i-1)*2*pi/(nbpoint_immersed))+center(1);
        X(i,2)=radius*sin((i-1)*2*pi/(nbpoint_immersed))+center(2);
end

for i=1:nbpoint_immersed

<<<<<<< 4b7ac153599ed2555be0cbaad76fd8b6752704b3
        X2(i,1)=radius*cos((i-1)*2*pi/(nbpoint_immersed))+center(1)+0.4;
        X2(i,2)=radius*sin((i-1)*2*pi/(nbpoint_immersed))+center(2);
=======
        X2(i,1)=1*radius*0.2*cos((i-1)*2*pi/(nbpoint_immersed))+center(1);
        X2(i,2)=1*radius*0.2*sin((i-1)*2*pi/(nbpoint_immersed))+center(2);
>>>>>>> Final_projet_stage
end

X=[X ;X2];
nbpoint_immersed=length(X);

for i=1:nbpoint_immersed
    immersed_value(i)=immersed_boundary_value_function(X(i,:));
end



convergence=[];
% base grid definition
mesh=cells_define(X_domain(1),X_domain(2),Y_domain(1),Y_domain(2),initial_Raffienement);
for j=1:nb_refinement

disp('size of the problem')
disp('refinement step')
disp(j)
disp('nb elements : ');
disp(length(mesh.points))
disp('nb elements unconstraint : ');
disp(length(mesh.points)-sum(mesh.hanging))
A_Sparse_index=[];
n_dof_domain=length(mesh.points);

parfor i=1:length(mesh.points)
    a_p=zeros(1,length(mesh.points));
    a_p(i)=1;
    if mesh.points(i,1)==0 | mesh.points(i,1)==1 | mesh.points(i,2)==0 | mesh.points(i,2)==1% | mesh.points(i,:)==[0.5 0.5]
        a_p=a_p;
    elseif mesh.hanging(i)==true
        a_p=get_stancil_hanging(mesh,i,1);
    else
        a_p=-get_stancil_L(mesh,i,1);
    end
    A_Sparse_index=[A_Sparse_index ;sparsematrix(a_p,i)];
end


cells_sum=eval_immersed_boundary_sum_in_cells(mesh,X);

B_immersed=[]
C_Sparse_index=[]

for i=1:length(cells_sum.boundary)
    
    a_p=zeros(1,n_dof_domain);
    a_p=get_lagrange_multiplier_stancil_weak_sum(mesh,cells_sum.index(i),X,cells_sum);
    B_immersed=[B_immersed ; cells_sum.sum(i)];    
    C_Sparse_index=[C_Sparse_index ;sparsematrix(a_p,i+n_dof_domain)];
    
end

C_Sparse_index=[C_Sparse_index;C_Sparse_index(:,2) C_Sparse_index(:,1) C_Sparse_index(:,3)];


A_Sparse_index=[A_Sparse_index ;C_Sparse_index];

A=sparse(A_Sparse_index(:,1),A_Sparse_index(:,2),A_Sparse_index(:,3));



B=zeros(length(mesh.points),1);
Reel=zeros(length(mesh.points),1);
for i=1:length(mesh.points)
    if mesh.points(i,1)==0 | mesh.points(i,1)==1 | mesh.points(i,2)==0 | mesh.points(i,2)==1
        B(i)=0;
    elseif mesh.hanging(i)==true
        B(i)=0;
%     elseif mesh.points(i,:)==[0.5 0.5]
%         B(i)=1;
    else
        B(i)=0;
    end
    
end

B=[B ;B_immersed];

B=sparse(B);


U=A\B;

U=full(U);
Solution=U(1:n_dof_domain);
T=delaunay(mesh.points);

trisurf(mesh.connect_active,mesh.points(:,1),mesh.points(:,2),Solution);
hold on
plot3(X(:,1),X(:,2),immersed_value,'.r','markersize',30);
 pause(0.001);
hold off
mesh_order=max(mesh.order);
mesh_size=1/(2^(mesh_order-1));
<<<<<<< 4b7ac153599ed2555be0cbaad76fd8b6752704b3
convergence=[convergence; 1/mesh_size max(Solution)-1 ];
=======


rout=0.35;
rin=rout*0.2;

Tin=1;
Tout=0;
%error estimation
error=[];
Mesh_position=[];
for i=1:length(mesh.points)
    if norm(mesh.points(i,:)-center)>rout*0.2 & norm(mesh.points(i,:)-center)<rout
        error=[error; Solution(i)-(Tout+log(norm(mesh.points(i,:)-center)/rout)/log(rin/rout)*(Tin-Tout))];
        Mesh_position=[ Mesh_position ; mesh.points(i,:)];

    end
    
end
figure()
plot3(Mesh_position(:,1) ,Mesh_position(:,2),error,'.')

relative_error_L2=norm(error,2)/sqrt(length(error))
relative_error_L1=mean(abs(error))
L_inf=max(abs(error))




convergence=[convergence; 1/sqrt(length(error)) relative_error_L2 relative_error_L1 L_inf ];
>>>>>>> Final_projet_stage

mesh_uniformity_index(j)=std(mesh.order(mesh.active==true));


if(j<nb_refinement)
mesh=raffine(mesh,Solution,0.3,0,cells_sum.index);
end


end





figure()
hold on
T=delaunay(mesh.points);
trisurf(T,mesh.points(:,1),mesh.points(:,2),Solution);
plot3(X(:,1),X(:,2),immersed_value,'.r','markersize',30);
hold off

<<<<<<< 4b7ac153599ed2555be0cbaad76fd8b6752704b3
figure()
convergence_order=zeros(length(convergence)-2,2);
for i=1:length(convergence)-2
convergence_order(i,2)=log((convergence(i,2)-convergence(i+1,2))/(convergence(i+1,2)-convergence(i+2,2)))/log(1.01);
convergence_order(i,1)=convergence(i+1,1);
end
plot(convergence_order(:,1),convergence_order(:,2));

figure()
plot(convergence(:,2));

=======
% figure()
% convergence_order=zeros(length(convergence)-2,3);
% for i=1:length(convergence)-2
%     nb_point_ratio=convergence(i+1,1)/convergence(i,1)*0.5+convergence(i+2,1)/convergence(i+1,1)*0.5
% convergence_order(i,2)=log((convergence(i,2)-convergence(i+1,2))/(convergence(i+1,2)-convergence(i+2,2)))/log(nb_point_ratio);
% convergence_order(i,3)=log((convergence(i,3)-convergence(i+1,3))/(convergence(i+1,3)-convergence(i+2,3)))/log(nb_point_ratio);
% convergence_order(i,1)=convergence(i+1,1);
% end
% plot(convergence_order(:,1),convergence_order(:,2),convergence_order(:,1),convergence_order(:,3));
% 

iter_X=log(convergence(:,1))';

iter_Y_L1=log(convergence(:,3))';
iter_Y_L2=log(convergence(:,2))';
iter_Y_Linf=log(convergence(:,4))';
line_converge_L1=polyfit(iter_X,iter_Y_L1,1)
line_converge_L2=polyfit(iter_X,iter_Y_L2,1)
line_converge_Linf=polyfit(iter_X,iter_Y_Linf,1)
approx_converge_L1=polyval(line_converge_L1,iter_X);
approx_converge_L2=polyval(line_converge_L2,iter_X);
approx_converge_Linf=polyval(line_converge_Linf,iter_X);
figure()
hold on
plot(iter_X,iter_Y_L1);
plot(iter_X,iter_Y_L2);
plot(iter_X,iter_Y_Linf);
plot(iter_X,approx_converge_L1);
plot(iter_X,approx_converge_L2);
plot(iter_X,approx_converge_Linf);
legend('L1_error','L2_error','Linf_error');
xlabel('ln(dx)');
ylabel('ln(error)');

hold off
>>>>>>> Final_projet_stage
 
    
    