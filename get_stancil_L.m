% code by Lucka Barbeau
function L_stancil=get_stancil_L(TR,k,order)
%give laplacien stancil for the point k in triangulation TR
neighbors=TR.neighbors_v{k,1};

L_stancil=zeros(1,length(TR.points));
N=length(neighbors);
dx_total=0;
dy_total=0;
m=1;
l=1;
order=2;
for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
        if dx_dy(1)~=0
    neighbors_egal_y(m)=neighbors(i);
    m=m+1;
    dx_total=abs(dx_dy(1))+dx_total;
        elseif dx_dy(2)~=0
    neighbors_egal_x(l)=neighbors(i);
    l=l+1;
    dy_total=abs(dx_dy(2))+dy_total;
        end
      
end

% l=1;
% m=1;
if l>=3 
stancil_dx=get_stancil_dx(TR,neighbors_egal_x(1),order)/2+get_stancil_dx(TR,neighbors_egal_x(2),order)/2;
else
stancil_dx=get_stancil_dx(TR,k,order);
end
if m>=3 
stancil_dy=get_stancil_dy(TR,neighbors_egal_y(1),order)/2+get_stancil_dy(TR,neighbors_egal_y(2),order)/2;
else
stancil_dy=get_stancil_dy(TR,k,order);
end


for i=1:length(neighbors)
        dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
        dx=dx_dy(1);
        dy=dx_dy(2);

        if dx~=0
            L_stancil=L_stancil+(get_stancil_dx(TR,neighbors(i),order)-stancil_dx)/dx*(abs(dx))/dx_total;
        end
        if dy~=0
            L_stancil=L_stancil+(get_stancil_dy(TR,neighbors(i),order)-stancil_dy)/dy*(abs(dy))/dy_total;
        end
end

end


