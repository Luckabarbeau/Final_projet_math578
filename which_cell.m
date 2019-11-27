function cell_index = which_cell(cells,X)
%WHICH_CELL Summary of this function goes here
% return the cell in which the point X is found in the cells mesh
cell_correspond=find(cells.boundary(:,1)<=X(1) & cells.boundary(:,2)>=X(1) & cells.boundary(:,3)<=X(2) & cells.boundary(:,4)>=X(2));
cell_index=cell_correspond(find(cells.active(cell_correspond)==true));

cell_index=cell_index(1);
end

