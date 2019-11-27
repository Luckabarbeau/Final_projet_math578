% code by Lucka Barbeau
function hanging_stancil=get_stancil_hanging(TR,k,square)
%give laplacien stancil for the point k in triangulation TR
if square==1
neighbors=V_neighbors_square(TR,k);
else
neighbors=V_neighbors(TR,k);
end 
hanging_stancil=zeros(1,length(TR.points));
N=length(neighbors);

m=1;
l=1;
for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
        if dx_dy(1)~=0
    neighbors_egal_y(m)=neighbors(i);
    m=m+1;
        elseif dx_dy(2)~=0
    neighbors_egal_x(l)=neighbors(i);
    l=l+1;
        end
      
end


if l>=3
hanging_stancil(k)=1   ;
hanging_stancil(neighbors_egal_x(1))=-1/2;
hanging_stancil(neighbors_egal_x(2))=-1/2;

else
hanging_stancil(k)=1   ;
hanging_stancil(neighbors_egal_y(1))=-1/2;
hanging_stancil(neighbors_egal_y(2))=-1/2;
end





end


