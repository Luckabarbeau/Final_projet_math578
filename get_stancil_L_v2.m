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
        
        neighbors_X(m)=neighbors(i);
        m=m+1;
        dx_total=abs(dx_dy(1))+dx_total;
    elseif dx_dy(2)~=0
        neighbors_Y(l)=neighbors(i);
        l=l+1;
        dy_total=abs(dx_dy(2))+dy_total;
    end
    
end



% l=1;
% m=1;
if m>=3
    
    
dx1=TR.points(neighbors_X(1),:)-TR.points(k,:);
dx1=dx1(1);
dx2=TR.points(neighbors_X(2),:)-TR.points(k,:);
dx2=dx2(1);
% 
% stancil_dx=get_stancil_dx(TR,neighbors_X(1),order)/2+get_stancil_dx(TR,neighbors_X(2),order)/2;
 stancil_dx=get_stancil_dx(TR,k,order);

stancil_dudx1=L_stancil;
stancil_dudx1(neighbors_X(1))=1/dx1;
stancil_dudx1(k)=-1/dx1;

stancil_dudx2=L_stancil;
stancil_dudx2(neighbors_X(2))=1/dx2;
stancil_dudx2(k)=-1/dx2;

stancil_dxx=(stancil_dudx1-stancil_dx)/(0.5*dx1)*abs(dx1)/dx_total+(stancil_dudx2-stancil_dx)/(0.5*dx2)*abs(dx2)/dx_total;
else
stancil_dx=get_stancil_dx(TR,k,order);
end

if l>=3 
dy1=TR.points(neighbors_Y(1),:)-TR.points(k,:);
dy1=dy1(2);
dy2=TR.points(neighbors_Y(2),:)-TR.points(k,:);
dy2=dy2(2);

% stancil_dy=get_stancil_dy(TR,neighbors_Y(1),order)/2+get_stancil_dy(TR,neighbors_Y(2),order)/2;
stancil_dy=get_stancil_dy(TR,k,order);

stancil_dudy1=L_stancil;
stancil_dudy1(neighbors_Y(1))=1/dy1;
stancil_dudy1(k)=-1/dy1;

stancil_dudy2=L_stancil;
stancil_dudy2(neighbors_Y(2))=1/dy2;
stancil_dudy2(k)=-1/dy2;

stancil_dyy=(stancil_dudy1-stancil_dy)/(0.5*dy1)*abs(dy1)/dy_total+(stancil_dudy2-stancil_dy)/(0.5*dy2)*abs(dy2)/dy_total;
else
stancil_dy=get_stancil_dy(TR,k,order);
end

L_stancilX=L_stancil;
L_stancilY=L_stancil;

% l=1;
% m=1;
for i=1:length(neighbors)
        dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
        dx=dx_dy(1);
        dy=dx_dy(2);

        if dx~=0 & m<3
            L_stancilX=L_stancilX+(get_stancil_dx(TR,neighbors(i),order)-stancil_dx)/dx*(abs(dx))/dx_total;
        end
        if dy~=0 & l<3
            L_stancilY=L_stancilY+(get_stancil_dy(TR,neighbors(i),order)-stancil_dy)/dy*(abs(dy))/dy_total;
        end
end


if l>=3 & m>=3
    L_stancil=stancil_dyy+stancil_dxx;
elseif m>=3 & l<3
    L_stancil=L_stancilY+stancil_dxx;
elseif m<3 & l>=3
    L_stancil=L_stancilX+stancil_dyy;
else
    L_stancil=L_stancilY+L_stancilX;
end

end


