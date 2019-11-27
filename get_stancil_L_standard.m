% code by Lucka Barbeau
function L_stancil=get_stancil_L_standard(TR,k,order)
%give laplacien stancil for the point k in triangulation TR

neighbors=V_neighbors_square(TR,k);

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



stancil_dyy=L_stancil;
stancil_dxx=L_stancil;


hy=norm(TR.points(neighbors_egal_x(1),:)-TR.points(k,:));
hx=norm(TR.points(neighbors_egal_y(1),:)-TR.points(k,:));


if l>=3 
stancil_dyy=L_stancil;
stancil_dyy(neighbors_egal_x(1))=1;
stancil_dyy(neighbors_egal_x(2))=1;
stancil_dyy(k)=-2;
else
    neighbors_x=V_neighbors_square(TR,neighbors_egal_x(1));
    for i=1:length(neighbors_x)
        dx_dy_2=TR.points(neighbors_x(i),:)-TR.points(neighbors_egal_x(1),:);
        if dx_dy_2(1)~=0
            if neighbors_x(i)~=k
                neighbors_egal_x(2)=neighbors_x(i);
                neighbors_x_2=V_neighbors_square(TR,neighbors_egal_x(2));
                for j=1:length(neighbors_x_2)
                    dx_dy_2=TR.points(neighbors_x_2(j),:)-TR.points(neighbors_egal_x(2),:);
                    if dx_dy_2(1)~=0
                        if neighbors_x(j)~=neighbors_egal_x(1)
                            neighbors_egal_x(3)=neighbors_x_2(j);
                        end
                    end
                end
            end
        end
    end
stancil_dyy(neighbors_egal_x(1))=-5;
stancil_dyy(neighbors_egal_x(2))=4;
stancil_dyy(neighbors_egal_x(3))=-1;
stancil_dyy(k)=2 ;   
end
if m>=3 
stancil_dxx=L_stancil;
stancil_dxx(neighbors_egal_y(1))=1;
stancil_dxx(neighbors_egal_y(2))=1;
stancil_dxx(k)=-2;
else
    neighbors_y=V_neighbors_square(TR,neighbors_egal_y(1));
    for i=1:length(neighbors_y)
        dx_dy_2=TR.points(neighbors_y(i),:)-TR.points(neighbors_egal_y(1),:);
        if dx_dy_2(1)~=0
            if neighbors_y(i)~=k
                neighbors_egal_y(2)=neighbors_y(i);
                neighbors_y_2=V_neighbors_square(TR,neighbors_egal_y(2));
                for j=1:length(neighbors_y_2)
                    dx_dy_2=TR.points(neighbors_y_2(i),:)-TR.points(neighbors_egal_y(2),:);
                    if dx_dy_2(1)~=0
                        if neighbors_y(i)~=neighbors_egal_y(1)
                            neighbors_egal_y(3)=neighbors_y(j);
                        end
                    end
                end
            end
        end
    end
stancil_dxx(neighbors_egal_y(1))=-5;
stancil_dxx(neighbors_egal_y(2))=4;
stancil_dxx(neighbors_egal_y(3))=-1;
stancil_dxx(k)=2 ; 

end

L_stancil=stancil_dxx/hx^2+stancil_dyy/hy^2;

end


