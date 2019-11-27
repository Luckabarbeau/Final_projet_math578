% code by Lucka Barbeau
% try one  heat equation Laplacien U=g(x)whit boundary condition=0

%g(x)=1 good approximation of the solution at the center is  0.073670467524337
%theorical convergence order of the method 2 for a full  grid mesh
clc ; clear all; close all

X_domain=[0 1];
Y_domain=[0 1];

% nb point 
% number on grid N by N grid

initial_Raffienement=5
initial_immersed_refinement=4
nb_refinement=5

alpha=1

nbpoint_immersed=10000;
radius=0.1;
center=[0.5 0.5];
X=[]
X2=[]
for i=1:nbpoint_immersed

        X(i,1)=radius*cos((i-1)*2*pi/(nbpoint_immersed))+center(1);
        X(i,2)=radius*sin((i-1)*2*pi/(nbpoint_immersed))+center(2);
end

% for i=1:nbpoint_immersed
%     if i<=nbpoint_immersed
%         X2(i,1)=0;
%         X2(i,2)=(i-1)/(nbpoint_immersed-1);
%     end
% end
% 
% for i=1:nbpoint_immersed
%     if i<=nbpoint_immersed
%         X3(i,1)=1-0.0001;
%         X3(i,2)=(i-1)/(nbpoint_immersed-1);
%     end
% end
% X=[X;X2;X3]



% base grid definition
mesh=cells_define(X_domain(1),X_domain(2),Y_domain(1),Y_domain(2),initial_Raffienement);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,1,X)
mesh=update_neighbors_v(mesh);







disp('size of the problem')
disp('refinement step')
disp(j)
disp('nb elements : ');
disp(length(mesh.points))
disp('nb elements unconstraint : ');
disp(length(mesh.points)-sum(mesh.hanging))
A_Sparse_index=[];
n_dof_domain=length(mesh.points);

for i=1:length(mesh.points)
    a_p=alpha*get_stancil_L_v2(mesh,i,1);
    A_Sparse_index=[A_Sparse_index ;sparsematrix(a_p,i)];
end

cells_sum=eval_immersed_boundary_sum_in_cells(mesh,X);
B_immersed=[];
C_Sparse_index=[];

k=0
% %%%boundary
for i=1:length(mesh.points)
    
    if mesh.points(i,1)==X_domain(1) | mesh.points(i,1)==X_domain(2) | mesh.points(i,1)==Y_domain(1) | mesh.points(i,1)==Y_domain(2)
        a_p=zeros(1,length(mesh.points));
        a_p(i)=10;
        k=k+1;
        C_Sparse_index=[C_Sparse_index ;sparsematrix(a_p,k)]; 
    end
    
end

%%% immersed boundary
for i=1:length(cells_sum.boundary)
    a_p=get_lagrange_multiplier_stancil_weak_sum(mesh,cells_sum.index(i),X,cells_sum);
    C_Sparse_index=[C_Sparse_index ;sparsematrix(a_p,i+k)];
end



A=sparse(A_Sparse_index(:,1),A_Sparse_index(:,2),A_Sparse_index(:,3));
C=sparse(C_Sparse_index(:,1),C_Sparse_index(:,2),C_Sparse_index(:,3),k+length(cells_sum.boundary),length(mesh.points));
A_global=[A C'; -C*A -C*C'];

for i=1:length(mesh.points)+length(cells_sum.boundary)+k
%     if mesh.points(i,1)>0.3 & mesh.points(i,1)<0.7 & mesh.points(i,2)>0.3 & mesh.points(i,2)<0.7 
if i<= length(mesh.points)
    if norm(mesh.points(i,:)-[0.5 0.5])<0.3
        U(i)=1;
    else
        U(i)=0;
    end
else
        U(i)=0;
end

end



Solution=U';
mindt=0.00001;
relative_tol=0.001;
max_time_step_modification=0.05
dt=mindt;
t=0;
tmax=2
Boundary=[X_domain Y_domain]
while t<tmax

[Solution, max_du_dt]=Runge_kutta(Solution,t,dt,A_global,mesh,Boundary)  ;
 

T=delaunay(mesh.points);
trisurf(T,mesh.points(:,1),mesh.points(:,2),Solution(1:length(mesh.points)));
hold on
plot3(X(:,1),X(:,2),ones(length(X),1),'.r');
hold off
pause(0.0000001)
t=t+dt
delta_ratio=max_du_dt/(max(abs(Solution)+0.001));
dt_modifier=delta_ratio/relative_tol;
if delta_ratio/relative_tol<max_time_step_modification*2
dt=max(mindt,dt*(1+max_time_step_modification-dt_modifier))
else
dt=max(mindt,dt*(1-max_time_step_modification))
end

end











 
    
    