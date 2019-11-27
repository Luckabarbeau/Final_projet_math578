% resoltion of the steady state navier stocke equation

%boundary condition
 clear all; %close all;

X_domain=[0 0.5];
Y_domain=[0 0.3];

rho=4;   
mu=0.16;
gx=0;
gy=0;

vx_in=100
vx_out=1

nbpoint_immersed=10000;
radius=0.02;
center=[0.15 0.15];
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


initial_Raffienement=6
initial_immersed_refinement=1
nb_refinement=1

mesh_size_suggestion=radius*renold_number^(-3/4);
mesh_order_suggestion=log2(renold_number^(3/4));

load_mesh=0
load_mesh_file='mesh_circle_27540_elipse_0_05_0_03_c_015_015';
load_mesh_file=['Meshs/' load_mesh_file];
if load_mesh==1
    load(load_mesh_file)
else
mesh=cells_define(X_domain(1),X_domain(2),Y_domain(1),Y_domain(2),initial_Raffienement);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,36,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,24,X);
mesh=raffine_immersed(mesh,initial_immersed_refinement,12,X);
mesh=raffine_immersed(mesh,initial_immersed_refinement,6,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,6,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,5,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,5,X);
mesh=update_neighbors_v(mesh);
end

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
cells_sum=eval_immersed_boundary_sum_in_cells(mesh,X);
U_x=ones(length(mesh.points),1)*vx_in;
U_y=ones(length(mesh.points),1)*0;

for i = 1:length(mesh.points)
    if isempty(find(mesh.connect(cells_sum.index,:)==i))==0
%     if norm(mesh.points(i,:)-center)<radius*2
        U_x(i)=0;
        U_y(i)=0;
    end
end



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
        a_p=-get_stancil_dx(mesh,mesh_position,2);
        a_p2=-get_stancil_dy(mesh,mesh_position,2);
        a_p3=mu/rho*get_stancil_L_v2(mesh,mesh_position,2);
        
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
            a_p=-get_stancil_dx(mesh,mesh_position,2)/rho;
    else
            a_p=-get_stancil_dy(mesh,mesh_position,2)/rho;
    end

%     a_p=[zeros(1,length(mesh.points)*2) a_p];
    C_Sparse_index=[C_Sparse_index ; sparsematrix(a_p,i)];  
end

P_Sparse_index=C_Sparse_index;
C_Sparse_index=[];
parfor i=1:length(mesh.points)
        a_p1=get_stancil_dx(mesh,i,2)/rho;
        a_p2=get_stancil_dy(mesh,i,2)/rho;
        a_p=[a_p1 a_p2];
        C_Sparse_index=[C_Sparse_index ; sparsematrix(a_p,i)]; 
end
% 
% if isempty(C_Sparse_index)
% else
% P_Sparse_index=[C_Sparse_index(:,2) C_Sparse_index(:,1) C_Sparse_index(:,3)];
% end


B_immersed=[];
D_Sparse_index=[];

parfor i=1:length(cells_sum.boundary)
    a_p=zeros(1,length(mesh.points));
    a_p=get_lagrange_multiplier_stancil_weak_sum(mesh,cells_sum.index(i),X,cells_sum);
    B_immersed=[B_immersed ; cells_sum.sum(i)];
    a_p2=[zeros(1,length(mesh.points)) a_p];
    D_Sparse_index=[D_Sparse_index ;sparsematrix(a_p,i)];
    D_Sparse_index=[D_Sparse_index ;sparsematrix(a_p2,i+length(cells_sum.boundary))];
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
        H_Sparse_index=[H_Sparse_index ; sparsematrix(a_p,l+2*length(cells_sum.boundary))];
        l=l+1;
        H_Sparse_index=[H_Sparse_index ; sparsematrix(a_p2,l+2*length(cells_sum.boundary))];
        l=l+1;
        a_p3=[zeros(1,length(a_p)*2),a_p];
        H_Sparse_index=[H_Sparse_index ; sparsematrix(a_p3,l+2*length(cells_sum.boundary))];
    end
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
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+2*length(cells_sum.boundary))];
        a_p2(1,i+length(mesh.points))=1; %speed value Y;
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+2*length(cells_sum.boundary))];
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
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+2*length(cells_sum.boundary))];
        a_p2(1,i+length(mesh.points))=1; %speed value Y
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+2*length(cells_sum.boundary))];
    elseif mesh.points(i,2)==Y_domain(2)
        a_p(1,i)=1; %speed value X
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+2*length(cells_sum.boundary))];
        a_p2(1,i+length(mesh.points))=1; %speed value Y
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+2*length(cells_sum.boundary))];
    elseif mesh.points(i,1)==X_domain(1)
        a_p(1,i)=1; %speed value X
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+2*length(cells_sum.boundary))];
        a_p2(1,i+length(mesh.points))=1; %speed value Y
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+2*length(cells_sum.boundary))];
%         a_p=zeros(1,3*length(mesh.points));
%         a_p(1,i+length(mesh.points)*2)=1; %speed value Y
%         k=k+1;
%         E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
    elseif  mesh.points(i,1)==X_domain(2)
        a_p=zeros(1,3*length(mesh.points));
        a_p(1,i+length(mesh.points)*2)=1; %speed value Y
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+2*length(cells_sum.boundary))];
      end
    
end



A_Sparse_m_dx=sparse(A_Sparse_index_dx(:,1),A_Sparse_index_dx(:,2),A_Sparse_index_dx(:,3));
A_Sparse_m_dy=sparse(A_Sparse_index_dy(:,1),A_Sparse_index_dy(:,2),A_Sparse_index_dy(:,3));


A_sparse_matrix_index=[A_Sparse_index ];
A_sparse_matrix=sparse(A_sparse_matrix_index(:,1),A_sparse_matrix_index(:,2),A_sparse_matrix_index(:,3));

Constraint_sparse_matrix_index=[D_Sparse_index;E_Sparse_index;H_Sparse_index];
Constraint_sparse_matrix=sparse(Constraint_sparse_matrix_index(:,1),Constraint_sparse_matrix_index(:,2),Constraint_sparse_matrix_index(:,3),l+k+2*length(cells_sum.boundary),3*length(mesh.points));

Constraint_sparse_conservation=sparse(C_Sparse_index(:,1),C_Sparse_index(:,2),C_Sparse_index(:,3),1*length(mesh.points),2*length(mesh.points));
Constraint_sparse_pressure=sparse(P_Sparse_index(:,1),P_Sparse_index(:,2),P_Sparse_index(:,3),2*length(mesh.points),length(mesh.points));



A_1=[A_sparse_matrix Constraint_sparse_pressure; Constraint_sparse_conservation sparse(length(mesh.points),length(mesh.points))];

A_Sparse_global=[A_1 Constraint_sparse_matrix' ; Constraint_sparse_matrix sparse(l+k+2*length(cells_sum.boundary),l+k+2*length(cells_sum.boundary))];


if iter==1
    Solution=[Solution ; zeros(length(cells_sum.boundary)*2,1); zeros(k,1); zeros(l,1)];
end
% 
% load('Solution')
% Solution(1:length(mesh.points)*3)=S2.Solution(1:length(mesh.points)*3)
time_sol=Solution;


mindt=0.000001;
relative_tol=0.001;
max_time_step_modification=0.1
dt=mindt;
t=0;
tmax=2
Boundary=[X_domain(1) Y_domain(1) Y_domain(2)]
I=speye(length(A_Sparse_global));
back_euler_M=(I-A_Sparse_global*dt);
A_Sparse_m_L=sparse(A_Sparse_index_L(:,1),A_Sparse_index_L(:,2),A_Sparse_index_L(:,3));

B=[zeros(length(mesh.points),1);B_immersed;B_immersed ;hanging_rhs ;boundary];
nb_zeros=l+k+2*length(cells_sum.boundary)

while t<tmax
    
U_x=Solution(1:length(mesh.points));
U_y=Solution(length(mesh.points)+1:length(mesh.points)*2);
P=Solution(length(mesh.points)*2+1:length(mesh.points)*3);
% U_x_a=avreage_whit_neighbors(mesh,U_x);
% U_y_a=avreage_whit_neighbors(mesh,U_y);
% Pa=avreage_whit_neighbors(mesh,P);
% Solution(1:length(mesh.points))=U_x_a;
% Solution(length(mesh.points)+1:length(mesh.points)*2)=U_y_a;
% Solution(length(mesh.points)*2+1:length(mesh.points)*3)=Pa;
% Solution(length(mesh.points)*2+1:length(mesh.points)*3)=Pa;

if mod(iter,2)==0
speed=sqrt(U_x.^2+U_y.^2);

U_x_u=U_x./sqrt(U_x.^2+U_y.^2);
U_y_u=U_y./sqrt(U_x.^2+U_y.^2);
figure(1)
quiver(mesh.points(:,1),mesh.points(:,2),U_x,U_y,'color',[1 0 1]);
hold on

if isempty(X)==0
plot(X(:,1),X(:,2),'.r');
end
trisurf(mesh.connect_active,mesh.points(:,1),mesh.points(:,2),-speed);
hold off
pause(0.0001)
end



D_Solution=Runge_kutta_dudt(Solution,t,dt,A_Sparse_global);
Solution(1:length(mesh.points)*2)=Solution(1:length(mesh.points)*2)+D_Solution(1:length(mesh.points)*2);
U_x=Solution(1:length(mesh.points));
U_y=Solution(length(mesh.points)+1:length(mesh.points)*2);
P=Solution(length(mesh.points)*2+1:length(mesh.points)*3);

D_Solution_S=A_Sparse_global*Solution;
A_Sparse_global=update_matrix_NS(U_x,U_y,mesh,A_Sparse_index_dx,A_Sparse_index_dy,A_Sparse_m_L,Constraint_sparse_pressure,Constraint_sparse_conservation,Constraint_sparse_matrix,nb_zeros);

RHS=[D_Solution_S(1:length(mesh.points)*2);B];
Solution=A_Sparse_global\RHS;
Solution(isnan(Solution))=0;
Solution(isinf(Solution))=0;

time_sol=[time_sol Solution];
t=t+dt



A_Sparse_global=update_matrix_NS(U_x,U_y,mesh,A_Sparse_index_dx,A_Sparse_index_dy,A_Sparse_m_L,Constraint_sparse_pressure,Constraint_sparse_conservation,Constraint_sparse_matrix,nb_zeros);



iter=iter+1


end





if(j<nb_refinement)
mesh=raffine(mesh,speed,0.0,0,cells_sum.index);
end
mesh=update_neighbors_v(mesh);

% relaxation=relaxation

end
