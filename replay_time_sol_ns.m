% solution Replay
%   load('time_sol_circle_reynold1000')
 nb=length(mesh.points);
%  nb=26213
nb_time=size(time_sol,2)
for i=1:nb_time
    Solution=time_sol(:,i);
    U_x=Solution(1:nb);
    U_y=Solution(nb+1:nb*2);
    P=Solution(nb*2+1:nb*3);
U_x_a=avreage_whit_neighbors(mesh,U_x);
U_y_a=avreage_whit_neighbors(mesh,U_y);
P=avreage_whit_neighbors(mesh,P);
% P=avreage_whit_neighbors(mesh,P);
speed=sqrt(U_x_a.^2+U_y_a.^2);

i=i




% U_x_u=U_x./sqrt(U_x.^2+U_y.^2);
% U_y_u=U_y./sqrt(U_x.^2+U_y.^2);
select_rand=rand(nb,1);
Select=select_rand<0.02;

figure(1)
quiver(mesh.points(Select,1),mesh.points(Select,2),U_x(Select),U_y(Select),'color',[1 0 0]);
% quiver(mesh.points(:,1),mesh.points(:,2),U_x(:),U_y(:),'color',[1 0 0]);
hold on

% if isempty(X)==0
% plot(X(:,1),X(:,2),'.r');
% 
% plot(center_p1(i,1),center_p1(i,2),'or');
% plot(center_p2(i,1),center_p2(i,2),'or');
% end

trisurf(mesh.connect_active,mesh.points(:,1),mesh.points(:,2),-speed);
colorbar

daspect([1 1 1])
xlim([min(mesh.points(:,1)) max(mesh.points(:,1))])
ylim([min(mesh.points(:,2)) max(mesh.points(:,2))])
% xlim([0.14 0.32])
% ylim([0.2 0.3])
shading interp;
hold off
pause(0.0001)


% U_dx_sparse_m=sparse(A_Sparse_index_dx(:,1),A_Sparse_index_dx(:,2),A_Sparse_index_dx(:,3));
% U_dy_sparse_m=sparse(A_Sparse_index_dy(:,1),A_Sparse_index_dy(:,2),A_Sparse_index_dy(:,3));
% U=[U_x;U_y];
% %%%%%% a valider
% %%%%%% %%%%%%%%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% dU_dx=U_dx_sparse_m*U;
% dU_dy=U_dy_sparse_m*U;
% dux_dy=dU_dy(1:length(mesh.points));
% duy_dx=dU_dx(length(mesh.points)+1:length(mesh.points)*2);
% select=rand(length(X),1)<0.02;
% [F_v_x, F_v_y]=force_on_boundary_viscosity(mesh,dux_dy,duy_dx,mu,X(select,:));
% [F_p_x, F_p_y]=force_on_boundary_pressure(mesh,P,X(select,:));
% 
% CD(i)=(F_v_x+F_p_x)*2/((exp(i/50)-1)^2*radius*rho);
% CL(i)=(F_v_y+F_p_y)*2/((exp(i/50)-1)^2*radius*rho);
% reynolds(i)=(exp(i/50)-1)*radius*rho/mu;
% 
% figure(3)
% loglog(reynolds,CD)
% hold on
% loglog(reynolds,CL)
% hold off




end
