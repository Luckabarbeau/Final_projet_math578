% code by Lucka Barbeau
function cells=raffine_immersed(cells,order,neighbor_order,X)
if isempty(X)==0
for k=1:order
    cell_immersed=find_immersed_boundary_point_in_cells(cells,X)
    refine_flag=[];
    if rand()<1
        for i=1:length(cell_immersed.index)
            to_refine=neighbors(cells,cell_immersed.index(i),neighbor_order);
            for j=1:length(to_refine)
                if isempty(find(refine_flag==to_refine(j)))==1
                    refine_flag=[refine_flag ; to_refine(j) ];
                end
            end
        end
    end
    
    for i=1:length(refine_flag)
        cells=raffine_cell(cells,refine_flag(i));
    end





cells.immersed_iteration=zeros(length(cells.immersed_iteration),1);
active=find(cells.active==1);
cells.connect_active=[];
for i=1:length(active)
    cells.connect_active(i,:)=cells.connect(active(i),:);
end
cells.active_list=active;
end
end
end

