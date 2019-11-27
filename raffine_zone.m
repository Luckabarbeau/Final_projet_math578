% code by Lucka Barbeau
function cells=raffine_zone(cells,Boundary)
refine_flag=[];
for i=1:length(cells.active_list)
    for j=1:4
        point=cells.points(cells.connect(cells.active_list(i),j),:);
        if point(1)>=Boundary(1) & point(1)<=Boundary(2) & point(2)>=Boundary(3)& point(2)<=Boundary(4)
                    refine_flag=[refine_flag ; cells.active_list(i) ];
        end
        break
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

