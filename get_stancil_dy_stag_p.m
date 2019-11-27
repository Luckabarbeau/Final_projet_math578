% code by Lucka Barbeau
function dy_stancil=get_stancil_dy_stag_p(TR,k)
%give dy stancil for the point k in triangulation TR

neighbors=TR.neighbors_v{k,1};

dy_stancil=zeros(1,length(TR.active_list));
N=length(neighbors);
dy_total=0;
for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    if dx_dy(2)~=0 && dx_dy(1)==0
        dy_total=(abs(dx_dy(2)))+dy_total;
    end
end
cells_parts_of=mod(find(TR.connect_active==k),length(TR.connect_active));
cells_parts_of(cells_parts_of==0)=length(TR.connect_active);
if length(cells_parts_of)~=4
    
for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    dy=dx_dy(2);
    if dy~=0
        dy_stancil=dy_stancil-get_stancil_eval_from_stag(TR,k)/dy*(abs(dy))/dy_total;
        dy_stancil=dy_stancil+get_stancil_eval_from_stag(TR,neighbors(i))/dy*(abs(dy))/dy_total;
    end
end
else
    for i=1:4
        dy=TR.boundary(TR.active_list(cells_parts_of(i)),4)-TR.boundary(TR.active_list(cells_parts_of(i)),3);
        if TR.center(TR.active_list(cells_parts_of(i)),1)>TR.points(k,1) & TR.center(TR.active_list(cells_parts_of(i)),2)>TR.points(k,2)
            dy_stancil(cells_parts_of(i))=0.5/dy;
        elseif TR.center(TR.active_list(cells_parts_of(i)),1)<TR.points(k,1) & TR.center(TR.active_list(cells_parts_of(i)),2)>TR.points(k,2)   
            dy_stancil(cells_parts_of(i))=0.5/dy;
        elseif TR.center(TR.active_list(cells_parts_of(i)),1)<TR.points(k,1) & TR.center(TR.active_list(cells_parts_of(i)),2)<TR.points(k,2)
            dy_stancil(cells_parts_of(i))=-0.5/dy;
        elseif TR.center(TR.active_list(cells_parts_of(i)),1)>TR.points(k,1) & TR.center(TR.active_list(cells_parts_of(i)),2)<TR.points(k,2)
            dy_stancil(cells_parts_of(i))=-0.5/dy;
        end
    end
end

end