function dy_stancil=get_stancil_dy_stag_u(TR,k)
%give dy stancil for the point k in triangulation TR


dy=TR.boundary(TR.active_list(k),4)-TR.boundary(TR.active_list(k),3);

dy_stancil=zeros(1,length(TR.points));

dy_stancil(TR.connect_active(k,1))=-0.5/dy;
dy_stancil(TR.connect_active(k,2))=-0.5/dy;
dy_stancil(TR.connect_active(k,3))=0.5/dy;
dy_stancil(TR.connect_active(k,4))=0.5/dy;


end