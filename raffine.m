% code by Lucka Barbeau

function cells=raffine(cells,U,ratio,neighbor_order,cell_immersed)

active_cell_error_estimation=error_estimator(cells,U);

refine_flag=[];
nb_refine=ceil(length(cells.connect_active)*ratio);
for i=1:nb_refine
    to_refine=neighbors(cells,active_cell_error_estimation(i,2),neighbor_order);
    for j=1:length(to_refine)
        if isempty(find(refine_flag==to_refine(j)))==1
            refine_flag=[refine_flag;to_refine(j)];
        end
    end
end


if rand()<1
    for i=1:length(cell_immersed)
        to_refine=neighbors(cells,cell_immersed(i),neighbor_order+1);
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

