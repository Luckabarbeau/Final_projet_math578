% code by Lucka Barbeau
function dy_stancil=get_stancil_dy_cross(TR,k,order)
%give dy stancil for the point k in triangulation TR

neighbors=TR.neighbors_v{k,1};


dy_stancil=zeros(1,length(TR.points));
N=length(neighbors);
dy_total=0;
for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    if dx_dy(2)~=0 & dx_dy(1)==0
        dy_total=(abs(dx_dy(2)))+dy_total;
    end
end





for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    dy=dx_dy(2);
    r=norm(dx_dy);
    if dy~=0
    dy_stancil(k)=dy_stancil(k)+-1/dy*(abs(dy))/dy_total;
    dy_stancil(neighbors(i))=dy_stancil(neighbors(i))+1/dy*(abs(dy))/dy_total;
    end
end

if N==65
    neighbors_2=zeros(4,1);
    j=0;
for i=3:4
    if length(TR.neighbors_v{neighbors(i),1})==4
        j=j+1;
        neighbors_2(j,1)= TR.neighbors_v{neighbors(i),1}(1); 
        j=j+1;
        neighbors_2(j,1)= TR.neighbors_v{neighbors(i),1}(2);
    end
end
if j==4
   dy_stancil2=zeros(1,length(TR.points));
   dy1=TR.points(neighbors(3),2)-TR.points(k,2);
   dy2=TR.points(neighbors(4),2)-TR.points(k,2);
   dy_tot=abs(dy1)+abs(dy2);
   
%    dy_stancil2(neighbors_2(1))=1/2/dy1*abs(dy1)/dy_tot;
%    dy_stancil2(neighbors_2(2))=1/2/dy1*abs(dy1)/dy_tot;
   
   dy_stancil2(neighbors_2(3))=1/2/dy2*abs(dy2)/dy_tot;
   dy_stancil2(neighbors_2(4))=1/2/dy2*abs(dy2)/dy_tot;
   
   dy_stancil2(neighbors(1))=-1/2/dy2*abs(dy2)/dy_tot;
   dy_stancil2(neighbors(2))=-1/2/dy2*abs(dy2)/dy_tot;
   
   dy_stancil2(neighbors(3))=1/dy1*abs(dy1)/dy_tot;
%    dy_stancil2(neighbors(4))=1/3/dy2*abs(dy2)/dy_tot;
   
   dy_stancil2(k)=-1/dy1*abs(dy1)/dy_tot;%+-1/dy2*abs(dy2)/dy_tot;
   
   dy_stancil=dy_stancil2;
    
end

end

end