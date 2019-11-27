% resoltion of the steady state navier stocke equation

%boundary condition
 clear all; %close all;

X_domain=[0 1];
Y_domain=[0 0.5];

rho=1000;   
mu=0.001;
gx=0;
gy=0;

vx_in=0.015
vx_out=1

nbpoint_immersed=0000;
radius=0.02;
center=[0.25 0.25];
X=[];
randfactor=0.0 
%cercles
for i=1:nbpoint_immersed

        X(i,1)=radius*cos((i-1)*2*pi/(nbpoint_immersed))+center(1);
        X(i,2)=radius*sin((i-1)*2*pi/(nbpoint_immersed))+center(2)+rand()*randfactor-randfactor/2;
end

for i=1:nbpoint_immersed

        X2(i,1)=1.01*radius*cos((i-1)*2*pi/(nbpoint_immersed))+center(1);
        X2(i,2)=1.01*radius*sin((i-1)*2*pi/(nbpoint_immersed))+center(2)+rand()*randfactor-randfactor/2;
end



for i=1:nbpoint_immersed

        X3(i,1)=1.01*radius*cos((i-1)*2*pi/(nbpoint_immersed))+center(1)+1;
        X3(i,2)=1.01*radius*sin((i-1)*2*pi/(nbpoint_immersed))+center(2)+rand()*randfactor-randfactor/2;
end
% % theta=-pi/16
% X=X*[cos(theta) -sin(theta); sin(theta) cos(theta)]


renold_number=vx_in*rho*2*radius/mu


initial_Raffienement=5
initial_immersed_refinement=1
nb_refinement=1

mesh_size_suggestion=radius*renold_number^(-3/4);
mesh_order_suggestion=log2(renold_number^(3/4));

load_mesh=0
load_mesh_file='mesh_circles_65560_1_05_C_025_025';
load_mesh_file=['Meshs/' load_mesh_file];
if load_mesh==1
    load(load_mesh_file)
else
mesh=cells_define(X_domain(1),X_domain(2),Y_domain(1),Y_domain(2),initial_Raffienement);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,36,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,24,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,12,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,6,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,6,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,5,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,5,X);
mesh=update_neighbors_v(mesh);
end
% mesh.points=mesh.points*0.1
% % mesh=raffine_immersed(mesh,initial_immersed_refinement,6,X);
% % mesh=raffine_immersed(mesh,initial_immersed_refinement,6,X);
% mesh=update_neighbors_v(mesh);
max_iter=1000

% X=[]
%re=60 cd cercle =1.13 mesh order =8  mesh order =7 cd =1.38
% 40 =ok 200 reponse osiclle
% min(10/renold_number,1);
relaxation=1;
figure(1);
figure(2);

figure(3);
% figure(4);
% figure(5);
for j=1:nb_refinement
error=1
iter=1

U_x=ones(length(mesh.points),1)*vx_in;


U_y=ones(length(mesh.points),1)*0;


P=zeros(length(mesh.points),1)*1;
Solution=[U_x; U_y; P];



A_Sparse_index=[];
A_Sparse_index_dx=[];
A_Sparse_index_dy=[];
A_Sparse_index_L=[];
parfor i=1:length(mesh.points)*2
    if i<=length(mesh.points)
        mesh_position=i;
    else
        mesh_position=i-length(mesh.points);
    end
        a_p=rho*get_stancil_dx(mesh,mesh_position,2);
        a_p2=rho*get_stancil_dy(mesh,mesh_position,2);
        a_p3=-mu*get_stancil_L_v2(mesh,mesh_position,2);
        
    if i<=length(mesh.points)
        a_p=[a_p zeros(1,length(mesh.points))];
        a_p2=[a_p2 zeros(1,length(mesh.points))];
        a_p3=[a_p3 zeros(1,length(mesh.points))];
    else
        a_p=[zeros(1,length(mesh.points)) a_p];
        a_p2=[zeros(1,length(mesh.points)) a_p2];
        a_p3=[zeros(1,length(mesh.points)) a_p3];
        
    end
    A_p=Solution(mesh_position)*a_p+Solution(mesh_position+length(mesh.points))*a_p2+a_p3
    A_Sparse_index_dx=[A_Sparse_index_dx ;sparsematrix(a_p,i)];
    A_Sparse_index_dy=[A_Sparse_index_dy ;sparsematrix(a_p2,i)];
    A_Sparse_index_L=[A_Sparse_index_L ;sparsematrix(a_p3,i)];
    A_Sparse_index=[A_Sparse_index;sparsematrix(A_p,i)];
end


C_Sparse_index=[];
parfor i=1:length(mesh.points)*2
    if i<=length(mesh.points)
        mesh_position=i;
    else
        mesh_position=i-length(mesh.points);
    end
        a_p=zeros(1,length(mesh.points));
    if i<=length(mesh.points) 
            a_p=get_stancil_dx(mesh,mesh_position,2);
    else
            a_p=get_stancil_dy(mesh,mesh_position,2);
    end
%        end
    a_p=[zeros(1,length(mesh.points)*2) a_p];
    C_Sparse_index=[C_Sparse_index ; sparsematrix(a_p,i)];  
end
 P_Sparse_index=C_Sparse_index;
parfor i=1:length(mesh.points)
        a_p1=get_stancil_dx(mesh,i,2);
        a_p2=get_stancil_dy(mesh,i,2);
        a_p=[a_p1 a_p2];
        C_Sparse_index=[C_Sparse_index ; sparsematrix(a_p,i+2*length(mesh.points))]; 
end


% if isempty(C_Sparse_index)
% else
% P_Sparse_index=[C_Sparse_index(:,2) C_Sparse_index(:,1) C_Sparse_index(:,3)];
% C_Sparse_index=[C_Sparse_index;C_Sparse_index(:,2) C_Sparse_index(:,1) C_Sparse_index(:,3)];
% end



cells_sum=eval_immersed_boundary_sum_in_cells(mesh,X);

B_immersed=[];
D_Sparse_index=[];

parfor i=1:length(cells_sum.boundary)
    a_p=zeros(1,length(mesh.points));
    a_p=get_lagrange_multiplier_stancil_weak_sum(mesh,cells_sum.index(i),X,cells_sum);
    B_immersed=[B_immersed ; cells_sum.sum(i)];
    a_p2=[zeros(1,length(mesh.points)) a_p];
    D_Sparse_index=[D_Sparse_index ;sparsematrix(a_p,i+length(mesh.points)*3)];
    D_Sparse_index=[D_Sparse_index ;sparsematrix(a_p2,i+length(mesh.points)*3+length(cells_sum.boundary))];
end

if isempty(D_Sparse_index)
else
D_Sparse_index=[D_Sparse_index;D_Sparse_index(:,2) D_Sparse_index(:,1) D_Sparse_index(:,3)];
end




B=zeros(length(mesh.points)*3,1);
for i=1:length(mesh.points)

         B(i)=0;
    
end

hanging_rhs=[];
%%hanging nods
parfor i=1:length(mesh.points)
    if mesh.hanging(i)==true & mesh.points(i,1)~=X_domain(1) & mesh.points(i,1)~=X_domain(2)& mesh.points(i,2)~=Y_domain(1)& mesh.points(i,2)~=Y_domain(2)
        hanging_rhs=[hanging_rhs ;0] ; % hanging ux 
        hanging_rhs=[hanging_rhs ;0] ;% hangign uy
        hanging_rhs=[hanging_rhs ;0];% hangign p
    end
end
 H_Sparse_index=[];
 l=0;
for i=1:length(mesh.points)
    if mesh.hanging(i)==true & mesh.points(i,1)~=X_domain(1) & mesh.points(i,1)~=X_domain(2)& mesh.points(i,2)~=Y_domain(1)& mesh.points(i,2)~=Y_domain(2)
        a_p=get_stancil_hanging(mesh,i);
        l=l+1;    
        a_p2=[zeros(1,length(a_p)),a_p];
        
        H_Sparse_index=[H_Sparse_index ; sparsematrix(a_p,l+3*length(mesh.points)+2*length(cells_sum.boundary))];
        l=l+1;
        H_Sparse_index=[H_Sparse_index ; sparsematrix(a_p2,l+3*length(mesh.points)+2*length(cells_sum.boundary))];
        l=l+1;
        a_p3=[zeros(1,length(a_p)*2),a_p];
        H_Sparse_index=[H_Sparse_index ; sparsematrix(a_p3,l+3*length(mesh.points)+2*length(cells_sum.boundary))];
    end
end

if isempty(H_Sparse_index)
else
H_Sparse_index=[H_Sparse_index;H_Sparse_index(:,2) H_Sparse_index(:,1) H_Sparse_index(:,3)];
end
%boundary condition

boundary=[];
for i=1:length(mesh.points)
     if mesh.points(i,1)==X_domain(1) & (mesh.points(i,2)==Y_domain(1) | mesh.points(i,2)==Y_domain(2))
         boundary=[boundary; vx_in]; %speed value X
         boundary=[boundary; 0]; %speed value Y
%          boundary=[boundary; vx_in];
%     elseif mesh.points(i,1)==X_domain(2) & (mesh.points(i,2)==Y_domain(1) | mesh.points(i,2)==Y_domain(2))
%          boundary=[boundary; 0]; %speed value X
%          boundary=[boundary; 0]; %speed value Y


%     elseif mesh.points(i,1)==X_domain(2) 
%          boundary=[boundary; vx_in] %speed value X
%          boundary=[boundary; 0] %speed value Y
    elseif mesh.points(i,2)==Y_domain(1) 
         boundary=[boundary; vx_in]; %speed value X
         boundary=[boundary; 0]; %speed value Y
    elseif mesh.points(i,2)==Y_domain(2) 
         boundary=[boundary; vx_in]; %speed value X
         boundary=[boundary; 0]; %speed value Y    
     elseif mesh.points(i,1)==X_domain(1)
         boundary=[boundary; vx_in]; %speed value X
%          boundary=[boundary; (-mesh.points(i,2)^2+mesh.points(i,2))*vx_in*6]; %speed value X
         boundary=[boundary; 0]; %speed value Y
%          boundary=[boundary; vx_in]; %speed value P
    elseif  mesh.points(i,1)==X_domain(2)
         boundary=[boundary; 0]; %speed value P
     end
end
E_Sparse_index=[];
k=0;
for i=1:length(mesh.points)
    a_p=zeros(1,2*length(mesh.points));
    a_p2=zeros(1,2*length(mesh.points));
    if mesh.points(i,1)==X_domain(1) & (mesh.points(i,2)==Y_domain(1) | mesh.points(i,2)==Y_domain(2))
        a_p(1,i)=1; %speed value X
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
        a_p2(1,i+length(mesh.points))=1; %speed value Y;
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
%         a_p=zeros(1,3*length(mesh.points));
%         a_p(1,i+length(mesh.points)*2)=1; %speed value Y
%         k=k+1;
%         E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
%         
%     elseif mesh.points(i,1)==X_domain(2) & (mesh.points(i,2)==Y_domain(1) | mesh.points(i,2)==Y_domain(2))
%         a_p(1,i)=1; %speed value X
%         k=k+1;
%         E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
%         a_p2(1,i+length(mesh.points))=1; %speed value Y
%         k=k+1;
%         E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
%         
%         %     elseif mesh.points(i,1)==X_domain(2)
        %          a_p(1,i)=1; %speed value X
        %          k=k+1;
        %          E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
        %          a_p2(1,i+length(mesh.points))=1; %speed value Y
        %          k=k+1;
        %          E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
    elseif mesh.points(i,2)==Y_domain(1)
        a_p(1,i)=1; %speed value X
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
        a_p2(1,i+length(mesh.points))=1; %speed value Y
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
    elseif mesh.points(i,2)==Y_domain(2)
        a_p(1,i)=1; %speed value X
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
        a_p2(1,i+length(mesh.points))=1; %speed value Y
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
    elseif mesh.points(i,1)==X_domain(1)
        a_p(1,i)=1; %speed value X
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
        a_p2(1,i+length(mesh.points))=1; %speed value Y
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
%         a_p=zeros(1,3*length(mesh.points));
%         a_p(1,i+length(mesh.points)*2)=1; %speed value Y
%         k=k+1;
%         E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
    elseif  mesh.points(i,1)==X_domain(2)
        a_p=zeros(1,3*length(mesh.points));
        a_p(1,i+length(mesh.points)*2)=1; %speed value Y
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
     end
    
end

if isempty(E_Sparse_index)
else
E_Sparse_index=[E_Sparse_index;E_Sparse_index(:,2) E_Sparse_index(:,1) E_Sparse_index(:,3)];
end

A_Sparse_m_dx=sparse(A_Sparse_index_dx(:,1),A_Sparse_index_dx(:,2),A_Sparse_index_dx(:,3));
A_Sparse_m_dy=sparse(A_Sparse_index_dy(:,1),A_Sparse_index_dy(:,2),A_Sparse_index_dy(:,3));



A_Sparse_index_global=[A_Sparse_index ;C_Sparse_index;D_Sparse_index;E_Sparse_index;H_Sparse_index];
%A_Sparse_index=[A_Sparse_index ;C_Sparse_index;D_Sparse_index;E_Sparse_index];
A=sparse(A_Sparse_index_global(:,1),A_Sparse_index_global(:,2),A_Sparse_index_global(:,3));
A_L=sparse(A_Sparse_index_L(:,1),A_Sparse_index_L(:,2),A_Sparse_index_L(:,3));
B=[B; B_immersed ; B_immersed ;hanging_rhs ;boundary];
B=sparse(B);
if iter==1
    Solution=[Solution ; zeros(length(cells_sum.boundary)*2,1); zeros(k,1); zeros(l,1)];
end


matrix_size=size(A)
A_Sparse_m_L=sparse(A_Sparse_index_L(:,1),A_Sparse_index_L(:,2),A_Sparse_index_L(:,3),matrix_size(1),matrix_size(2));
while error>0.0000001 & iter<max_iter
U_x_a=U_x;
U_y_a=U_y;
Solution_A=Solution;
U=A\B;
U(isnan(U))=0;
U(isnan(U))=0;
U(isinf(U))=0;
U(isinf(U))=0;
Solution=(1-relaxation)*Solution+relaxation*U;

U_x=Solution(1:length(mesh.points));
U_y=Solution(length(mesh.points)+1:length(mesh.points)*2);
P=Solution(length(mesh.points)*2+1:length(mesh.points)*3);
% U_x_av=avreage_whit_neighbors(mesh,U_x);
% U_y_av=avreage_whit_neighbors(mesh,U_y);
Pa=avreage_whit_neighbors(mesh,P);
% Solution(1:length(mesh.points))=U_x_av;
% Solution(length(mesh.points)+1:length(mesh.points)*2)=U_y_av;
% Solution(length(mesh.points)*2+1:length(mesh.points)*3)=Pa;
% U_x=U_x_av;
% U_y=U_y_av;
% P=Pa;

error(iter)=max(norm(U_x-U_x_a,2),norm(U_x-U_x_a,2))/sqrt(length(error));
figure(2)
semilogy(linspace(1,iter,iter),error);
speed=sqrt(U_x.^2+U_y.^2);

U_x_u=U_x./sqrt(U_x.^2+U_y.^2);
U_y_u=U_y./sqrt(U_x.^2+U_y.^2);
figure(1)
quiver(mesh.points(:,1),mesh.points(:,2),U_x_u,U_y_u,'color',[1 0 1]);
hold on

if isempty(X)==0
plot(X(:,1),X(:,2),'.r');
end
trisurf(mesh.connect_active,mesh.points(:,1),mesh.points(:,2),-speed);
hold off


% P_sparse_m=sparse(P_Sparse_index(:,1),P_Sparse_index(:,2)-2*length(mesh.points),P_Sparse_index(:,3));
% U_dx_sparse_m=sparse(A_Sparse_index_dx(:,1),A_Sparse_index_dx(:,2),A_Sparse_index_dx(:,3));
% U_dy_sparse_m=sparse(A_Sparse_index_dy(:,1),A_Sparse_index_dy(:,2),A_Sparse_index_dy(:,3));
% U=[U_x;U_y];
% %%%%%% a valider
% %%%%%% %%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% dU_dx=U_dx_sparse_m*U;
% dU_dy=U_dy_sparse_m*U;
% dux_dy=dU_dy(1:length(mesh.points));
% duy_dx=dU_dx(length(mesh.points)+1:length(mesh.points)*2);
% 
% [F_v_x, F_v_y]=force_on_boundary_viscosity(mesh,dux_dy,duy_dx,mu,X);
% [F_p_x, F_p_y]=force_on_boundary_pressure(mesh,Pa,X2);
% 
% CD(iter)=(F_v_x+F_p_x)*2/(vx_in^2*radius*rho);
% CL(iter)=(F_v_y+F_p_y)*2/(vx_in^2*radius*rho);

% if iter>=2
%     error(iter)=max(norm(CD(iter-1)-CD(iter),2),norm(CL(iter-1)-CL(iter),2))/sqrt(length(error));
% else
%     error(iter)=1
% end

% figure(3)
% plot(CD)
% hold on
% plot(CL)
% hold off
% P_sparse_m=sparse(P_Sparse_index(:,1),P_Sparse_index(:,2)-2*length(mesh.points),P_Sparse_index(:,3));
% pressure_grad=P_sparse_m*P;
% pressure_grad_dx=pressure_grad(1:length(mesh.points));
% pressure_grad_dy=pressure_grad(length(mesh.points)+1:length(mesh.points)*2);
% figure(3)
% trisurf(T,mesh.points(:,1),mesh.points(:,2),pressure_grad_dx);
% figure(4)
% trisurf(T,mesh.points(:,1),mesh.points(:,2),pressure_grad_dy);
% figure(5)
% trisurf(T,mesh.points(:,1),mesh.points(:,2),P);
pause(0.0001)

A_Sparse_index_dx_iter=sparse_matrix_multiply_NS(mesh,A_Sparse_index_dx,U_x);
A_Sparse_index_dy_iter=sparse_matrix_multiply_NS(mesh,A_Sparse_index_dy,U_y);
A_dx=sparse(A_Sparse_index_dx_iter(:,1),A_Sparse_index_dx_iter(:,2),A_Sparse_index_dx_iter(:,3),matrix_size(1),matrix_size(2));
A_dy=sparse(A_Sparse_index_dy_iter(:,1),A_Sparse_index_dy_iter(:,2),A_Sparse_index_dy_iter(:,3),matrix_size(1),matrix_size(2));
A_momentum=A_dx+A_dy+A_Sparse_m_L;

A_Sparse_index_global=[C_Sparse_index;D_Sparse_index;E_Sparse_index;H_Sparse_index];
A=sparse(A_Sparse_index_global(:,1),A_Sparse_index_global(:,2),A_Sparse_index_global(:,3),matrix_size(1),matrix_size(2));
A=A+A_momentum;


iter=iter+1


end





if(j<nb_refinement)
mesh=raffine(mesh,speed,0.0,0,cells_sum.index);
end
mesh=update_neighbors_v(mesh);

% relaxation=relaxation

end
