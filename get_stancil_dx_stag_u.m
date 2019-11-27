% code by Lucka Barbeau
function dx_stancil=get_stancil_dx_stag_u(TR,k)
%give dy stancil for the point k in triangulation TR


dx=TR.boundary(TR.active_list(k),2)-TR.boundary(TR.active_list(k),1);

dx_stancil=zeros(1,length(TR.points));

dx_stancil(TR.connect_active(k,1))=-0.5/dx;
dx_stancil(TR.connect_active(k,2))=0.5/dx;
dx_stancil(TR.connect_active(k,3))=0.5/dx;
dx_stancil(TR.connect_active(k,4))=-0.5/dx;


end