function [F_p_x, F_p_y]= force_on_boundary_pressure(TR,P,boundary)
%FORCE_ON_BOUNDARY_VISCOSITY Summary of this function goes here
%   return the viscosity force on the boundary

parfor i=1:length(boundary)-1
    da(i,:)=[boundary(i+1,:)-boundary(i,:) (boundary(i+1,:)+boundary(i,:))/2];
end
da(length(boundary),:)=[boundary(1,:)-boundary(length(boundary),:) (boundary(1,:)+boundary(length(boundary),:))/2];
F_p_x=0;
F_p_y=0;
parfor  i=1:length(da)
    F_p_x=eval_solution_on_mesh(TR,da(i,3:4),P)*da(i,2)*-1+F_p_x;
    F_p_y=eval_solution_on_mesh(TR,da(i,3:4),P)*da(i,1)+F_p_y;
end


end

