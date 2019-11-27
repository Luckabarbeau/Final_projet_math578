% resoltion of the steady state navier stocke equation

%boundary condition
clear all

X_domain=[0 5];
Y_domain=[0 1];

rho=1;   
mu=0.1;
gx=0;
gy=0;

vx_in=100
vx_out=1

nbpoint_immersed=10000;
radius=0.05;
center=[2 .5];
X=[];
randfactor=0.0 
%cercles
for i=1:nbpoint_immersed

        X(i,1)=3*radius*cos((i-1)*2*pi/(nbpoint_immersed))+center(1);
        X(i,2)=radius*sin((i-1)*2*pi/(nbpoint_immersed))+center(2)+rand()*randfactor-randfactor/2;
end
%boundary side

% for i=1:nbpoint_immersed
%         X2(i,1)=(i-1)*X_domain(2)/(nbpoint_immersed);
%         X2(i,2)=Y_domain(1)+0.01;
% end
% for i=1:nbpoint_immersed
%         X3(i,1)=(i-1)*X_domain(2)/(nbpoint_immersed);
%         X3(i,2)=Y_domain(2)-0.01;
% end

% X=[X;X2;X3];
% theta=-pi/16
% X=X*[cos(theta) -sin(theta); sin(theta) cos(theta)]
renold_number=vx_in*rho*2*radius/mu




initial_Raffienement=7
initial_immersed_refinement=1
nb_refinement=1

mesh_size_suggestion=radius*renold_number^(-3/4);
mesh_order_suggestion=log2(renold_number^(3/4));

% initial_Raffienement=min(ceil(mesh_order_suggestion),initial_Raffienement);
mesh=cells_define(X_domain(1),X_domain(2),Y_domain(1),Y_domain(2),initial_Raffienement);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,3,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,3,X);
% mesh=raffine_immersed(mesh,initial_immersed_refinement,3,X);
mesh=update_neighbors_v(mesh);



max_iter=1000

% X=[]
% 40 =ok 200 reponse osiclle
% min(10/renold_number,1);
relaxation=1;
figure(1);
figure(2);

figure(3);
figure(4);
figure(5);
for j=1:nb_refinement
error=1
iter=1

U_x=ones(length(mesh.points),1)*vx_in;


U_y=ones(length(mesh.points),1)*0;


P=zeros(length(mesh.points),1)*1;
Solution=[U_x; U_y; P];

while error>0.0000001 & iter<max_iter

A_Sparse_index=[];
parfor i=1:length(mesh.points)*2
    a_p=zeros(1,length(mesh.points));
    
    if i<=length(mesh.points)
        mesh_position=i;
    else
        mesh_position=i-length(mesh.points);
    end
    
%     if mesh.hanging(mesh_position)==true & mesh.points(mesh_position,1)~=X_domain(1) & mesh.points(mesh_position,1)~=X_domain(2)& mesh.points(mesh_position,2)~=Y_domain(1)& mesh.points(mesh_position,2)~=Y_domain(2)
%         a_p=get_stancil_hanging(mesh,mesh_position);
%     else
        a_p=rho*(Solution(mesh_position)*get_stancil_dx(mesh,mesh_position,3)+Solution(mesh_position+length(mesh.points))*get_stancil_dy(mesh,mesh_position,3))-mu*get_stancil_L_v2(mesh,mesh_position,2);
%     end
%     
    if i<=length(mesh.points)
        a_p=[a_p zeros(1,length(mesh.points))];
    else
        a_p=[zeros(1,length(mesh.points)) a_p];
    end
    A_Sparse_index=[A_Sparse_index ;sparsematrix(a_p,i)];
end

C_Sparse_index=[];
parfor i=1:length(mesh.points)*2
    if i<=length(mesh.points)
        mesh_position=i;
    else
        mesh_position=i-length(mesh.points);
    end
        a_p=zeros(1,length(mesh.points));
%        if mesh.hanging(mesh_position)==true & mesh.points(mesh_position,1)~=X_domain(1) & mesh.points(mesh_position,1)~=X_domain(2)& mesh.points(mesh_position,2)~=Y_domain(1)& mesh.points(mesh_position,2)~=Y_domain(2)
%        else
    if i<=length(mesh.points) 
            a_p=get_stancil_dx(mesh,mesh_position,3);
    else
            a_p=get_stancil_dy(mesh,mesh_position,3);
    end
%        end
    a_p=[zeros(1,length(mesh.points)*2) a_p];
    C_Sparse_index=[C_Sparse_index ; sparsematrix(a_p,i)];  
end

P_Sparse_index=C_Sparse_index;
parfor i=1:length(mesh.points)
        a_p1=-get_stancil_dx(mesh,i,3);
        a_p2=-get_stancil_dy(mesh,i,3);
        a_p=[a_p1 a_p2];
        C_Sparse_index=[C_Sparse_index ; sparsematrix(a_p,i+2*length(mesh.points))]; 
end

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
% for i=1:length(mesh.points)
%     if mesh.hanging(i)==true & mesh.points(i,1)~=X_domain(1) & mesh.points(i,1)~=X_domain(2)& mesh.points(i,2)~=Y_domain(1)& mesh.points(i,2)~=Y_domain(2)
%         hanging_rhs=[hanging_rhs ;0] ; % hanging ux 
%         hanging_rhs=[hanging_rhs ;0] ;% hangign uy
%     end
% end
 H_Sparse_index=[];
 l=0;
% for i=1:length(mesh.points)
%     if mesh.hanging(i)==true & mesh.points(i,1)~=X_domain(1) & mesh.points(i,1)~=X_domain(2)& mesh.points(i,2)~=Y_domain(1)& mesh.points(i,2)~=Y_domain(2)
%         a_p=get_stancil_hanging(mesh,mesh_position);
%         l=l+1;    
%         a_p2=[zeros(1,length(a_p)),a_p];
%         H_Sparse_index=[H_Sparse_index ; sparsematrix(a_p,l+3*length(mesh.points)+2*length(cells_sum.boundary))];
%         l=l+1;
%         H_Sparse_index=[H_Sparse_index ; sparsematrix(a_p2,l+3*length(mesh.points)+2*length(cells_sum.boundary))];
%     end
% 
% end
% if isempty(H_Sparse_index)
% else
% H_Sparse_index=[H_Sparse_index;H_Sparse_index(:,2) H_Sparse_index(:,1) H_Sparse_index(:,3)];
% end
% %boundary condition

boundary=[];
for i=1:length(mesh.points)
     if mesh.points(i,1)==X_domain(1) & (mesh.points(i,2)==Y_domain(1) | mesh.points(i,2)==Y_domain(2))
         boundary=[boundary; 0]; %speed value X
         boundary=[boundary; 0]; %speed value Y
%          boundary=[boundary; vx_in];
    elseif mesh.points(i,1)==X_domain(2) & (mesh.points(i,2)==Y_domain(1) | mesh.points(i,2)==Y_domain(2))
         boundary=[boundary; 0]; %speed value X
         boundary=[boundary; 0]; %speed value Y


%     elseif mesh.points(i,1)==X_domain(2) 
%          boundary=[boundary; vx_in] %speed value X
%          boundary=[boundary; 0] %speed value Y
    elseif mesh.points(i,2)==Y_domain(1) 
         boundary=[boundary; 0]; %speed value X
         boundary=[boundary; 0]; %speed value Y
    elseif mesh.points(i,2)==Y_domain(2) 
         boundary=[boundary; 0]; %speed value X
         boundary=[boundary; 0]; %speed value Y    
     elseif mesh.points(i,1)==X_domain(1)
         boundary=[boundary; (-mesh.points(i,2)^2+mesh.points(i,2))*vx_in*6]; %speed value X
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
    elseif mesh.points(i,1)==X_domain(2) & (mesh.points(i,2)==Y_domain(1) | mesh.points(i,2)==Y_domain(2))
        a_p(1,i)=1; %speed value X
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
        a_p2(1,i+length(mesh.points))=1; %speed value Y
        k=k+1;
        E_Sparse_index=[E_Sparse_index ; sparsematrix(a_p2,l+k+3*length(mesh.points)+2*length(cells_sum.boundary))];
        
        %     elseif mesh.points(i,1)==X_domain(2)
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



 A_Sparse_index=[A_Sparse_index ;C_Sparse_index;D_Sparse_index;E_Sparse_index;H_Sparse_index];
%A_Sparse_index=[A_Sparse_index ;C_Sparse_index;D_Sparse_index;E_Sparse_index];
A=sparse(A_Sparse_index(:,1),A_Sparse_index(:,2),A_Sparse_index(:,3));

B=[B; B_immersed ; B_immersed ;hanging_rhs ;boundary];
B=sparse(B);
if iter==1
    Solution=[Solution ; zeros(length(cells_sum.boundary)*2,1); zeros(k,1); zeros(l,1)];
end
    
U_x_a=U_x;
U_y_a=U_y;
Solution_A=Solution;
U=A\B;
% U(isnan(U))=0;
% U(isnan(U))=0;
Solution=(1-relaxation)*Solution+relaxation*U;

U_x=Solution(1:length(mesh.points));
U_y=Solution(length(mesh.points)+1:length(mesh.points)*2);
P=Solution(length(mesh.points)*2+1:length(mesh.points)*3);
Pa=avreage_whit_neighbors(mesh,P);
error(iter)=max(norm(U_x-U_x_a,2),norm(U_x-U_x_a,2))/sqrt(length(error));
figure(2)
semilogy(linspace(1,iter,iter),error);
speed=sqrt(U_x.^2+U_y.^2);

figure(1)
quiver(mesh.points(:,1),mesh.points(:,2),U_x,U_y,'color',[1 0 1]);
hold on
T=delaunay(mesh.points);
if isempty(X)==0
plot(X(:,1),X(:,2),'.r');
end
trisurf(T,mesh.points(:,1),mesh.points(:,2),-speed);
hold off

iter=iter+1

P_sparse_m=sparse(P_Sparse_index(:,1),P_Sparse_index(:,2)-2*length(mesh.points),P_Sparse_index(:,3));
pressure_grad=P_sparse_m*P;
pressure_grad_dx=pressure_grad(1:length(mesh.points));
pressure_grad_dy=pressure_grad(length(mesh.points)+1:length(mesh.points)*2);
figure(3)
trisurf(T,mesh.points(:,1),mesh.points(:,2),pressure_grad_dx);
figure(4)
trisurf(T,mesh.points(:,1),mesh.points(:,2),pressure_grad_dy);
figure(5)
trisurf(T,mesh.points(:,1),mesh.points(:,2),P);
pause(0.0001)


end





if(j<nb_refinement)
mesh=raffine(mesh,speed,0.1,0,cells_sum.index);
end
mesh=update_neighbors_v(mesh);

% relaxation=relaxation

end

