% code by Lucka Barbeau
function center = neighbors(cells,center,order_of_neighbors)
%NEIGHBORS Summary of this function goes here
% find all cells that share a node whit the current cell Take the global
% index of the cell not the index in the active cells caterogie

neighbors_index=[];
points_from_centers=[];
if order_of_neighbors>=1
for i=1:order_of_neighbors
    for j=1:length(center)
        for k=1:4
            if isempty(find(points_from_centers==cells.connect(center(j),k)))
                points_from_centers=[points_from_centers ;cells.connect(center(j),k)];
            end
        end
    end
    
    for j=1:length(points_from_centers)
        check=mod(find(cells.connect_active==points_from_centers(j)),length(cells.connect_active));
        check(check==0)=length(cells.connect_active);

        for k=1:length(check)
            if isempty(find(neighbors_index==cells.active_list(check(k))))
                neighbors_index=[neighbors_index ; cells.active_list(check(k))];
            end
        end
    end
    neighbors_index(neighbors_index==0)=length(cells.connect);
    center=[center ;neighbors_index];

end
else
    center=center;
end
end

