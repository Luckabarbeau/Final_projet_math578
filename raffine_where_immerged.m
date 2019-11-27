% code by Lucka Barbeau

function cells=raffine_where_immerged(cells,cells_immerged,ratio,neighbor_order)


for i=1:length(cells_immerged)
    cells=raffine_cell(cells,cells_immerged(i));
end

cells.immersed_iteration=zeros(length(cells.immersed_iteration),1);
active=find(cells.active==1);
cells.connect_active=[];
for i=1:length(active)
    cells.connect_active(i,:)=cells.connect(active(i),:);
end
cells.active_list=active;

end

