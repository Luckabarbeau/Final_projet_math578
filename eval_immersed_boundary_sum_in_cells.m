function cells_boundary_sum= eval_immersed_boundary_sum_in_cells(cells,X,type)
%GENERATE_LAMBDA_INTERPOLATION creat a object whit all cells and the
%boundary points that they containt


if isempty(X)==0
cells_boundary_containt=find_immersed_boundary_point_in_cells(cells,X);
cells_boundary_sum=cells_boundary_containt;
cells_boundary_sum.sum=zeros(length(cells_boundary_sum.boundary),1);
for i=1:length(cells_boundary_sum.boundary)
        for j=1:length(cells_boundary_sum.boundary{i})
            cells_boundary_sum.sum(i,1)=cells_boundary_sum.sum(i,1)+immersed_boundary_value_function(X(cells_boundary_sum.boundary{i}(j),:));
        end
end
else
    cells_boundary_sum.sum=[];
    cells_boundary_sum.boundary=[];
    cells_boundary_sum.index=[];
end

end