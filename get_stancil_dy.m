% code by Lucka Barbeau
function dy_stancil=get_stancil_dy(TR,k,square)
%give dy stancil for the point k in triangulation TR
if square==1
neighbors=V_neighbors_square(TR,k);
else
neighbors=V_neighbors(TR,k);
end 
dy_stancil=zeros(1,length(TR.points));
N=length(neighbors);
dy_total=0;
for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    if dx_dy(2)~=0
        dy_total=(abs(dx_dy(2)))+dy_total;
    end
end

for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    dy=dx_dy(2);
    r=norm(dx_dy);
<<<<<<< 4b7ac153599ed2555be0cbaad76fd8b6752704b3
    if dy~=0
    dy_stancil(k)=dy_stancil(k)+-1/dy*(abs(dy))/dy_total;
    dy_stancil(neighbors(i))=dy_stancil(neighbors(i))+1/dy*(abs(dy))/dy_total;
=======
    if dy>0
    dy_stancil(k)=dy_stancil(k)+-1/dy;
    dy_stancil(neighbors(i))=dy_stancil(neighbors(i))+1/dy;
    break
    elseif dy<0
    dy_stancil(k)=dy_stancil(k)+-1/dy;
    dy_stancil(neighbors(i))=dy_stancil(neighbors(i))+1/dy;
    break
>>>>>>> Final_projet_stage
    end
end

end