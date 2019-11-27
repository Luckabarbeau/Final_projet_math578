% code by Lucka Barbeau
function dx_stancil=get_stancil_dx_stag_p(TR,k)
%give dy stancil for the point k in triangulation TR

neighbors=TR.neighbors_v{k,1};

dx_stancil=zeros(1,length(TR.active_list));
N=length(neighbors);
dx_total=0;
for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    if dx_dy(1)~=0 && dx_dy(2)==0
        dx_total=(abs(dx_dy(1)))+dx_total;
    end
end

cells_parts_of=mod(find(TR.connect_active==k),length(TR.connect_active));
cells_parts_of(cells_parts_of==0)=length(TR.connect_active);
if length(cells_parts_of)~=4
    
for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    dx=dx_dy(1);
    if dx~=0
        dx_stancil=dx_stancil-get_stancil_eval_from_stag(TR,k)/dx*(abs(dx))/dx_total;
        dx_stancil=dx_stancil+get_stancil_eval_from_stag(TR,neighbors(i))/dx*(abs(dx))/dx_total;
    end
end
else
    for i=1:4
        dx=TR.boundary(TR.active_list(cells_parts_of(i)),2)-TR.boundary(TR.active_list(cells_parts_of(i)),1);
        if TR.center(TR.active_list(cells_parts_of(i)),1)>TR.points(k,1) & TR.center(TR.active_list(cells_parts_of(i)),2)>TR.points(k,2)
            dx_stancil(cells_parts_of(i))=0.5/dx;
        elseif TR.center(TR.active_list(cells_parts_of(i)),1)<TR.points(k,1) & TR.center(TR.active_list(cells_parts_of(i)),2)>TR.points(k,2)   
            dx_stancil(cells_parts_of(i))=-0.5/dx;
        elseif TR.center(TR.active_list(cells_parts_of(i)),1)<TR.points(k,1) & TR.center(TR.active_list(cells_parts_of(i)),2)<TR.points(k,2)
            dx_stancil(cells_parts_of(i))=-0.5/dx;
        elseif TR.center(TR.active_list(cells_parts_of(i)),1)>TR.points(k,1) & TR.center(TR.active_list(cells_parts_of(i)),2)<TR.points(k,2)
            dx_stancil(cells_parts_of(i))=0.5/dx;
        end
    end
end


end