% code by Lucka Barbeau
function stag_stancil=get_stancil_eval_from_stag(TR,k)
%give dy stancil for the point k in triangulation TR

cells_parts_of=mod(find(TR.connect_active==k),length(TR.connect_active));
cells_parts_of(cells_parts_of==0)=length(TR.connect_active);
total_area=0;
stag_stancil=zeros(1,length(TR.active_list));
order=1;
if length(cells_parts_of)<=4
for i= 1:length(cells_parts_of)
total_area=total_area+TR.area(TR.active_list(cells_parts_of(i)))^order;
stag_stancil(cells_parts_of(i))=TR.area(TR.active_list(cells_parts_of(i)))^order;
% stag_stancil(cells_parts_of(i))=1;
end
stag_stancil=stag_stancil/total_area;
else


end
% stag_stancil=stag_stancil/i;
end