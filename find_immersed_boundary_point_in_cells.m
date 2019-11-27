function cells_boundary_containt= find_immersed_boundary_point_in_cells(cells,X)
%GENERATE_LAMBDA_INTERPOLATION creat a object whit all cells and the
%boundary points that they containt

k=1;
for i=1:length(cells.active_list)
    boundary_containt=find(X(:,1)>=cells.boundary(cells.active_list(i),1) & X(:,1)<cells.boundary(cells.active_list(i),2) & X(:,2)>=cells.boundary(cells.active_list(i),3) & X(:,2)<cells.boundary(cells.active_list(i),4));
    if isempty(boundary_containt)~=1
        cells_boundary_containt.index(k,1)=cells.active_list(i);
        cells_boundary_containt.boundary{k}=boundary_containt';
        k=k+1;
    end
end


    
end

