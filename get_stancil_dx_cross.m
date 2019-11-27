% code by Lucka Barbeau
function dx_stancil=get_stancil_dx_cross(TR,k,order)
%give dx stancil for the point k in triangulation TR

neighbors=TR.neighbors_v{k,1};
dx_stancil=zeros(1,length(TR.points));
N=length(neighbors);
dx_total=0;



for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    if dx_dy(1)~=0 & dx_dy(2)==0
        dx_total=abs(dx_dy(1))+dx_total;
    end
end



for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    dx=dx_dy(1);
    r=norm(dx_dy);
    if dx~=0 
    dx_stancil(k)=dx_stancil(k)+-1/dx*abs(dx)/dx_total;
    dx_stancil(neighbors(i))=dx_stancil(neighbors(i))+1/dx*abs(dx)/dx_total;
    end
end
if N==66
    neighbors_2=zeros(4,1);
    j=0;
for i=1:2
    if length(TR.neighbors_v{neighbors(i),1})==4
        j=j+1;
        neighbors_2(j,1)=TR.neighbors_v{neighbors(i),1}(3); 
        j=j+1;
        neighbors_2(j,1)= TR.neighbors_v{neighbors(i),1}(4);
    end
end
if j==4
   dx_stancil2=zeros(1,length(TR.points));
   dx1=TR.points(neighbors(1),1)-TR.points(k,1);
   dx2=TR.points(neighbors(2),1)-TR.points(k,1);
   dx_tot=abs(dx1)+abs(dx2);
   
%    dx_stancil2(neighbors_2(1))=1/dx1*abs(dx1)/dx_tot;
%    dx_stancil2(neighbors_2(2))=1/3/dx1*abs(dx1)/dx_tot;
   
   dx_stancil2(neighbors_2(3))=1/2/dx2*abs(dx2)/dx_tot;
   dx_stancil2(neighbors_2(4))=1/2/dx2*abs(dx2)/dx_tot;
   
    dx_stancil2(neighbors(1))=1/dx1*abs(dx1)/dx_tot;
%    dx_stancil2(neighbors(2))=1/3/dx2*abs(dx2)/dx_tot;
   
   dx_stancil2(neighbors(3))=-1/2/dx2*abs(dx2)/dx_tot;
   dx_stancil2(neighbors(4))=-1/2/dx2*abs(dx2)/dx_tot;
   
   dx_stancil2(k)=-1/dx1*abs(dx1)/dx_tot;%+-1/dx2*abs(dx2)/dx_tot;
   
   dx_stancil=dx_stancil2;
    
end

end
end