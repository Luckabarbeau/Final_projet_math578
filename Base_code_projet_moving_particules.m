% resoltion of the steady state navier stocke equation

%boundary condition
 clear all; %close all;

X_domain=[0 0.2];
Y_domain=[0 0.2];
rho_p=1500;
rho=1000;   
mu=0.001;
gx=0;
gy=0;
g=9.81
Contact_stiff=10000 
epsilone_contact=0.0005

vx_in=0.001;
vx_out=1;

nbpoint_immersed=1000;
radius=0.005;
center_p1=[0.1 0.1];
X=[];
randfactor=0.0 ;
%cercles
for i=1:nbpoint_immersed
        X(i,1)=1*radius*cos((i-1)*2*pi/(nbpoint_immersed))+center_p1(1);
        X(i,2)=1*radius*sin((i-1)*2*pi/(nbpoint_immersed))+center_p1(2)+rand()*randfactor-randfactor/2;
end
% center_p2=[0.11 0.185];
% for i=nbpoint_immersed+1:nbpoint_immersed*2
%         X(i,1)=1*radius*cos((i-1)*2*pi/(nbpoint_immersed))+center_p2(1);
%         X(i,2)=1*radius*sin((i-1)*2*pi/(nbpoint_immersed))+center_p2(2)+rand()*randfactor-randfactor/2;
% end

% particule initial speed
v=[0 0; 0 0];
m_p=rho_p*radius^2*pi;

    



renold_number=vx_in*rho*2*radius/mu;
disp('the reynolds number for that simulation is :')
disp(renold_number)
initial_Raffienement=7  ;
initial_immersed_refinement=1;
nb_refinement=1;

mesh_size_suggestion=radius*renold_number^(-3/4);
mesh_order_suggestion=log2(renold_number^(3/4));

load_mesh=0
load_mesh_file='mesh_circle_27540_elipse_0_05_0_03_c_015_015';
load_mesh_file=['Meshs/' load_mesh_file];
if load_mesh==1
    load(load_mesh_file)
else
disp('mesh initialisation...')    
mesh=cells_define(X_domain(1),X_domain(2),Y_domain(1),Y_domain(2),initial_Raffienement);
disp('mesh initialisation complete')
disp('mesh zones refinement...')
% mesh=raffine_zone(mesh,[0.05 0.15 0 0.4]);  
% mesh=raffine_zone(mesh,[0.0 0.7 0.13 0.37]);
% mesh=raffine_zone(mesh,[0.0 0.5 0.16 0.34]);
% mesh=raffine_zone(mesh,[0.08 0.13 0.18 0.32]);
% mesh=raffine_zone(mesh,[0.1 0.13 0.2 0.3]);
% mesh=raffine_zone(mesh,[0.15 0.35 0.20 0.3]);
% mesh=raffine_zone(mesh,[0.28 0.3 0.23 0.245]);
% mesh=raffine_zone(mesh,[0.28 0.3 0.23 0.245]);
% mesh=raffine_zone(mesh,[0.28 0.3 0.23 0.245]);
% mesh=raffine_zone(mesh,[0.28 0.3 0.23 0.245]);

disp('mesh zones refinement complete')
disp('mesh refinement arrond immersed boundary...')
% mesh=raffine_immersed(mesh,initial_immersed_refinement,12,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,10,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,8,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,6,X);
% mesh=raffine_zone(mesh,[0.297 0.299 0.232 0.234]);
% mesh=raffine_zone(mesh,[0.297 0.299 0.232 0.234]);
% mesh=raffine_zone(mesh,[0.298 0.29855 0.2326 0.2328]);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,2,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,2,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,2,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,2,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,2,X);
disp('mesh refinement arrond immersed boundary complete')
disp('update mesh connection...')
mesh=update_neighbors_v(mesh);
disp('update mesh connection complete')
end




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
error=1;
iter=1;

disp('define immersed boundary value...')

cells_sum=eval_immersed_boundary_sum_in_cells_NS(mesh,X,1);
cells_sum_uy=eval_immersed_boundary_sum_in_cells_NS(mesh,X,2);

disp('define immersed boundary value complet')
disp('define initial solution...')
U_x=ones(length(mesh.points),1)*vx_in;
U_y=ones(length(mesh.points),1)*0;

for i = 1:length(mesh.points)
    if isempty(find(mesh.connect(cells_sum.index,:)==i))==0
%     if norm(mesh.points(i,:)-center)<radius*2
        U_x(i)=0;
        U_y(i)=0;
    end
end



P=zeros(length(mesh.center),1)*1;
Solution=[U_x; U_y; P];
disp('define initial solution complete')


A_Sparse_index=[];
A_Sparse_index_dx=[];
A_Sparse_index_dy=[];
A_Sparse_index_L=[];
disp('define momentum matrix...')
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
disp('define momentum matrix complete')


disp('define pressure matrix...')
C_Sparse_index=[];
pressure_map_Sparse_index=[];
parfor i=1:length(mesh.points)*2
    if i<=length(mesh.points)
        mesh_position=i;
    else
        mesh_position=i-length(mesh.points);
    end
    
    
        
    if i<=length(mesh.points)
        a_p=get_stancil_dx_stag_p(mesh,mesh_position)/rho;
        
    else
        a_p=get_stancil_dy_stag_p(mesh,mesh_position)/rho;
        
    end
    if i<=length(mesh.points)
    pressure_map_Sparse_index=[pressure_map_Sparse_index ; sparsematrix(get_stancil_eval_from_stag(mesh,mesh_position),i)];
    end

    C_Sparse_index=[C_Sparse_index ; sparsematrix(a_p,i)];  
end

disp('define pressure matrix  complete')
P_Sparse_index=C_Sparse_index;
C_Sparse_index=[];
disp('define pressure matrix  conservation matrix...')
for i=1:length(mesh.active_list)
    if mesh.center(i,1)==max(mesh.center(:,1))
        disp('here')
    end
        a_p1=get_stancil_dx_stag_u(mesh,i);
        a_p2=get_stancil_dy_stag_u(mesh,i);
        a_p=[a_p1 a_p2];
        C_Sparse_index=[C_Sparse_index ; sparsematrix(a_p,i)]; 
end
disp('define pressure matrix  conservation matrix complete')
% 
% if isempty(C_Sparse_index)
% else
% P_Sparse_index=[C_Sparse_index(:,2) C_Sparse_index(:,1) C_Sparse_index(:,3)];
% end


B_immersed_ux=[];
B_immersed_uy=[];
D_Sparse_index=[];
disp('define immersed boundary condition...')
parfor i=1:length(cells_sum.boundary)
    a_p=zeros(1,length(mesh.points));
    a_p=get_lagrange_multiplier_stancil_weak_sum(mesh,cells_sum.index(i),X,cells_sum);
    B_immersed_ux=[B_immersed_ux ; cells_sum.sum(i)];
    B_immersed_uy=[B_immersed_uy ; cells_sum_uy.sum(i)];
    a_p2=[zeros(1,length(mesh.points)) a_p];
    D_Sparse_index=[D_Sparse_index ;sparsematrix(a_p,i)];
    D_Sparse_index=[D_Sparse_index ;sparsematrix(a_p2,i+length(cells_sum.boundary))];
end
disp('define immersed boundary condition complete')



B=zeros(length(mesh.points)*3,1);
for i=1:length(mesh.points)

         B(i)=0;

end
disp('define hanging node constraint...')
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
disp('define hanging node constraint complete')

%boundary condition
Boundary_for_time_position=[]
boundary=[];
disp('define boundary condition...')
for i=1:length(mesh.points)
     if mesh.points(i,1)==X_domain(1) & (mesh.points(i,2)==Y_domain(1) | mesh.points(i,2)==Y_domain(2))
%          boundary=[boundary; 0]; %speed value X
%          Boundary_for_time_position=[Boundary_for_time_position length(boundary)];
%          boundary=[boundary; 0]; %speed value Y
         boundary=[boundary; 0]; %speed value P
%          boundary=[boundary; vx_in];
%     if mesh.points(i,:)==[0.1 0.1] %norm(mesh.points(i,:)-center_p1)<radius %& mesh.points(i,1)>=X_domain(2)/4 &  mesh.points(i,1)<=X_domain(2)*3/4
%          boundary=[boundary; 0]; %speed value X
%          boundary=[boundary; 0]; %speed value Y


%     elseif mesh.points(i,1)==X_domain(2) 
%          boundary=[boundary; vx_in] %speed value X
%          boundary=[boundary; 0] %speed value Y
    elseif mesh.points(i,2)==Y_domain(1) 
         boundary=[boundary; vx_in]; %speed value X
         Boundary_for_time_position=[Boundary_for_time_position length(boundary)];
         boundary=[boundary; 0]; %speed value Y
    elseif mesh.points(i,2)==Y_domain(2) 
         boundary=[boundary; vx_in]; %speed value X
         Boundary_for_time_position=[Boundary_for_time_position length(boundary)];
         boundary=[boundary; 0]; %speed value Y    
     elseif mesh.points(i,1)==X_domain(1)
         boundary=[boundary; vx_in]; %speed value X
         Boundary_for_time_position=[Boundary_for_time_position length(boundary)];
%          boundary=[boundary; (-(mesh.points(i,2)/0.5)^2+mesh.points(i,2)/0.5)*vx_in*6]; %speed value X
         boundary=[boundary; 0]; %speed value Y
%          boundary=[boundary; vx_in]; %speed value P
    elseif  mesh.points(i,1)==X_domain(2)
         boundary=[boundary; vx_in]; %speed value X
         boundary=[boundary; 0]; %speed value Y   
%          boundary=[boundary; 0]; %speed value P
     end
end

E_Sparse_index=[];
k=0;
for i=1:length(mesh.points)
    a_p=zeros(1,2*length(mesh.points));
    a_p2=zeros(1,2*length(mesh.points));
    if mesh.points(i,1)==X_domain(1) & (mesh.points(i,2)==Y_domain(1) || mesh.points(i,2)==Y_domain(2))
%         a_p(1,i)=1; %speed value X
%         k=k+1;
%         E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+2*length(cells_sum.boundary))];
%         a_p2(1,i+length(mesh.points))=1; %speed value Y;
%         k=k+1;
%         E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+2*length(cells_sum.boundary))];
        a_p=zeros(1,2*length(mesh.points));
        a_p2=get_stancil_eval_from_stag(mesh,i);
        a_p=[a_p a_p2];%speed value Y
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+2*length(cells_sum.boundary))];
        
%     if mesh.points(i,:)==[0.1 0.1] %norm(mesh.points(i,:)-center_p1)<radius% & mesh.points(i,1)>=X_domain(2)/4 &  mesh.points(i,1)<=X_domain(2)*3/4
%         a_p(1,i)=1; %speed value X
%         k=k+1;
%         E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+2*length(cells_sum.boundary))];
%         a_p2(1,i+length(mesh.points))=1; %speed value Y
%         k=k+1;
%         E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+2*length(cells_sum.boundary))];
% %         
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
         a_p(1,i)=1; %speed value X
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+2*length(cells_sum.boundary))];
        a_p2(1,i+length(mesh.points))=1; %speed value Y
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+2*length(cells_sum.boundary))];
%         a_p=zeros(1,2*length(mesh.points));
%         a_p2=get_stancil_eval_from_stag(mesh,i);
%         a_p=[a_p a_p2];%speed value Y
%         k=k+1;
%         E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+2*length(cells_sum.boundary))];
       end
    
end
disp('define boundary condition complete')

disp('assembling matrix from sub matrix...')
A_Sparse_m_dx=sparse(A_Sparse_index_dx(:,1),A_Sparse_index_dx(:,2),A_Sparse_index_dx(:,3));
A_Sparse_m_dy=sparse(A_Sparse_index_dy(:,1),A_Sparse_index_dy(:,2),A_Sparse_index_dy(:,3));


A_sparse_matrix_index=[A_Sparse_index ];
A_sparse_matrix=sparse(A_sparse_matrix_index(:,1),A_sparse_matrix_index(:,2),A_sparse_matrix_index(:,3));

Constraint_sparse_matrix_index=[D_Sparse_index;E_Sparse_index;H_Sparse_index];
Constraint_sparse_matrix=sparse(Constraint_sparse_matrix_index(:,1),Constraint_sparse_matrix_index(:,2),Constraint_sparse_matrix_index(:,3),l+k+2*length(cells_sum.boundary),2*length(mesh.points)+length(mesh.active_list));

Constraint_sparse_conservation=sparse(C_Sparse_index(:,1),C_Sparse_index(:,2),C_Sparse_index(:,3),length(mesh.active_list),2*length(mesh.points));
Constraint_sparse_pressure=sparse(P_Sparse_index(:,1),P_Sparse_index(:,2),P_Sparse_index(:,3),2*length(mesh.points),length(mesh.active_list));



A_1=[A_sparse_matrix Constraint_sparse_pressure; Constraint_sparse_conservation sparse(length(mesh.active_list),length(mesh.active_list))];

A_Sparse_global=[A_1 Constraint_sparse_matrix' ; Constraint_sparse_matrix sparse(l+k+2*length(cells_sum.boundary),l+k+2*length(cells_sum.boundary))];
disp('assembling matrix from sub matrix complete')
U_dx_sparse_m=sparse(A_Sparse_index_dx(:,1),A_Sparse_index_dx(:,2),A_Sparse_index_dx(:,3));
U_dy_sparse_m=sparse(A_Sparse_index_dy(:,1),A_Sparse_index_dy(:,2),A_Sparse_index_dy(:,3));


pressure_map=sparse(pressure_map_Sparse_index(:,1),pressure_map_Sparse_index(:,2),pressure_map_Sparse_index(:,3));

if iter==1
    Solution=[Solution ; zeros(length(cells_sum.boundary)*2,1); zeros(k,1); zeros(l,1)];
end
% 
% load('Solution')
% Solution(1:length(mesh.points)*3)=S2.Solution(1:length(mesh.points)*3)
time_sol{1}=Solution

mindt=0.0025;
relative_tol=0.001;
max_time_step_modification=0.1;
dt=mindt;
t=0;
tmax=2
Boundary=[X_domain(1) Y_domain(1) Y_domain(2)];
I=speye(length(mesh.points)*2);
% back_euler_M=(I-A_Sparse_global*dt);
if isempty(A_Sparse_index_L)
    A_Sparse_m_L=sparse(length(mesh.points)*2);
else
A_Sparse_m_L=sparse(A_Sparse_index_L(:,1),A_Sparse_index_L(:,2),A_Sparse_index_L(:,3));
end
B=[zeros(length(mesh.active_list),1);B_immersed_ux;B_immersed_uy ;hanging_rhs ;boundary];
Boundary_time_offset=length(B)-length(boundary);
nb_zeros=l+k+2*length(cells_sum.boundary);
disp('Simulation will run for ')
disp(tmax)
disp('seconds')
disp('solving the equation in time...')
disp('time step being solve...')
disp(0)






while t<tmax




RHS=[Solution(1:2*length(mesh.points))/dt; B];
Solution=A_Sparse_global\RHS;

Solution(isnan(Solution))=0;
Solution(isinf(Solution))=0;

time_sol{length(time_sol)+1}=Solution;
t=t+dt;

I2=I/dt;

    
U_x=Solution(1:length(mesh.points));
U_y=Solution(length(mesh.points)+1:length(mesh.points)*2);
P=Solution(length(mesh.points)*2+1:length(mesh.points)*2+length(mesh.active_list));
P=pressure_map*P;
U_x_a=avreage_whit_neighbors(mesh,U_x);
U_y_a=avreage_whit_neighbors(mesh,U_y);
% P=avreage_whit_neighbors(mesh,P);
% Solution(1:length(mesh.points))=U_x_a;
% Solution(length(mesh.points)+1:length(mesh.points)*2)=U_y_a;
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
U=[U_x;U_y];
dU_dx=U_dx_sparse_m*U;
dU_dy=U_dy_sparse_m*U;
dux_dy=dU_dy(1:length(mesh.points));
duy_dx=dU_dx(length(mesh.points)+1:length(mesh.points)*2);
% X_p1=X(1:nbpoint_immersed,:);
% X_p2=X(nbpoint_immersed+1:nbpoint_immersed*2,:);
% 
% %contact force
% if norm(center_p1(iter,:)-center_p2(iter,:))<2*radius+epsilone_contact
%     s_l=2*radius+epsilone_contact-norm(center_p1(iter,:)-center_p2(iter,:))
%     FC=Contact_stiff*s_l*(center_p1(iter,:)-center_p2(iter,:))/norm(center_p1(iter,:)-center_p2(iter,:))
% else
%     FC=[0 0]
% end
% [F_v_x, F_v_y]=force_on_boundary_viscosity(mesh,dux_dy,duy_dx,mu,X_p1);
% [F_p_x, F_p_y]=force_on_boundary_pressure(mesh,P,X_p1);
% f_g=(rho-rho_p)*radius^2*pi*g;
% f_p=[F_p_x+F_v_x F_v_y+F_p_y+f_g];
% a=(f_p+FC)/m_p;
% 
% v(1,:)=v(1,:)+a*dt;
% 
% [F_v_x, F_v_y]=force_on_boundary_viscosity(mesh,dux_dy,duy_dx,mu,X_p2);
% [F_p_x, F_p_y]=force_on_boundary_pressure(mesh,P,X_p2);
% f_g=(rho-rho_p)*radius^2*pi*g;
% f_p=[F_p_x+F_v_x F_v_y+F_p_y+f_g];
% a=(f_p-FC)/m_p;
% v(2,:)=v(2,:)+a*dt;
% % v=[0 -0.1];
% center_p1(iter+1,:)=center_p1(iter,:)+v(1,:)*dt;
% center_p2(iter+1,:)=center_p2(iter,:)+v(2,:)*dt;
% X_p1=X_p1+[ones(nbpoint_immersed,1)*v(1,1) ones(nbpoint_immersed,1)*v(1,2)]*dt;
% X_p2=X_p2+[ones(nbpoint_immersed,1)*v(2,1) ones(nbpoint_immersed,1)*v(2,2)]*dt;
% X=[X_p1;X_p2];


B_immersed_ux=[];
B_immersed_uy=[];
D_Sparse_index=[];

nb_immersed_boundary=length(cells_sum.boundary);
cells_sum=eval_immersed_boundary_sum_in_cells_particule(mesh,X,1,v,0,nbpoint_immersed);
cells_sum_uy=eval_immersed_boundary_sum_in_cells_particule(mesh,X,2,v,0,nbpoint_immersed);

parfor i=1:length(cells_sum.boundary)
    a_p=zeros(1,length(mesh.points));
    a_p=get_lagrange_multiplier_stancil_weak_sum(mesh,cells_sum.index(i),X,cells_sum);
    B_immersed_ux=[B_immersed_ux ; cells_sum.sum(i)];
    B_immersed_uy=[B_immersed_uy ; cells_sum_uy.sum(i)];
    a_p2=[zeros(1,length(mesh.points)) a_p];
    D_Sparse_index=[D_Sparse_index ;sparsematrix(a_p,i)];
    D_Sparse_index=[D_Sparse_index ;sparsematrix(a_p2,i+length(cells_sum.boundary))];
end

nb_immersed_boundary=nb_immersed_boundary-length(cells_sum.boundary);
E_Sparse_index(:,1)=E_Sparse_index(:,1)-nb_immersed_boundary*2;
if isempty(H_Sparse_index)==0
H_Sparse_index(:,1)=H_Sparse_index(:,1)-nb_immersed_boundary*2;
end

Constraint_sparse_matrix_index=[D_Sparse_index;E_Sparse_index;H_Sparse_index];
Constraint_sparse_matrix=sparse(Constraint_sparse_matrix_index(:,1),Constraint_sparse_matrix_index(:,2),Constraint_sparse_matrix_index(:,3),l+k+2*length(cells_sum.boundary),2*length(mesh.points)+length(mesh.active_list));

nb_zeros=l+k+2*length(cells_sum.boundary);
A_Sparse_global=update_matrix_NS(U_x,U_y,mesh,A_Sparse_index_dx,A_Sparse_index_dy,A_Sparse_m_L,Constraint_sparse_pressure,Constraint_sparse_conservation,Constraint_sparse_matrix,nb_zeros,I2);

B=[zeros(length(mesh.active_list),1);B_immersed_ux;B_immersed_uy ;hanging_rhs ;boundary];

iter=iter+1;
disp('time step being solve...')
disp(t)

end





if(j<nb_refinement)
mesh=raffine(mesh,speed,0.0,0,cells_sum.index);
end
mesh=update_neighbors_v(mesh);

% relaxation=relaxation
end

