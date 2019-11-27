function [F_v_x, F_v_y]= force_on_boundary_viscosity(TR,dux_dy,duy_dx,mu,boundary)
%FORCE_ON_BOUNDARY_VISCOSITY Summary of this function goes here
%   return the viscosity force on the boundary

parfor i=1:length(boundary)-1
    da(i,:)=[boundary(i+1,:)-boundary(i,:) (boundary(i+1,:)+boundary(i,:))/2];
end
da(length(boundary),:)=[boundary(1,:)-boundary(length(boundary),:) (boundary(1,:)+boundary(length(boundary),:))/2];
F_v_x=0;
F_v_y=0;
parfor  i=1:length(da)
    F_v_x=eval_solution_on_mesh(TR,da(i,3:4),dux_dy)*da(i,1)*-1*mu+F_v_x;
    F_v_y=eval_solution_on_mesh(TR,da(i,3:4),duy_dx)*da(i,2)*mu+F_v_y;
end


end

