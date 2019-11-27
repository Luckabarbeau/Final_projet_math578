% code by Lucka Barbeau
function dx_stancil=get_stancil_dx(TR,k,square)
%give dx stancil for the point k in triangulation TR
if square==1
neighbors=V_neighbors_square(TR,k);
else
neighbors=V_neighbors(TR,k);
end  
dx_stancil=zeros(1,length(TR.points));
N=length(neighbors);
dx_total=0;
for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    if dx_dy(1)~=0
    dx_total=abs(dx_dy(1))+dx_total;
    end
end

for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    dx=dx_dy(1);
    r=norm(dx_dy);
<<<<<<< 4b7ac153599ed2555be0cbaad76fd8b6752704b3
    if dx~=0
    dx_stancil(k)=dx_stancil(k)+-1/dx*abs(dx)/dx_total;
    dx_stancil(neighbors(i))=dx_stancil(neighbors(i))+1/dx*abs(dx)/dx_total;
=======
    if dx>0
    dx_stancil(k)=dx_stancil(k)+-1/dx;
    dx_stancil(neighbors(i))=dx_stancil(neighbors(i))+1/dx;
    break
    elseif dx<0
    dx_stancil(k)=dx_stancil(k)+-1/dx;
    dx_stancil(neighbors(i))=dx_stancil(neighbors(i))+1/dx;
    break
>>>>>>> Final_projet_stage
    end
end

end