% code by Lucka Barbeau

function active_cell_error_estimation = error_estimator(cells,U)
%ERROR_ESTIMATOR Summary of this function goes here

for  i=1:length(cells.connect_active)
    point_cell=[U(cells.connect_active(i,1)),U(cells.connect_active(i,2)),U(cells.connect_active(i,3)),U(cells.connect_active(i,4))];
    area=norm(cells.points(cells.connect_active(i,1),:)-cells.points(cells.connect_active(i,2),:))*norm(cells.points(cells.connect_active(i,2),:)-cells.points(cells.connect_active(i,3),:));
    active_cell_error_estimation(i,:)=[std(point_cell)*area,cells.active_list(i)];
end
active_cell_error_estimation=sortrows(active_cell_error_estimation,'descend');

end

