function dY = clear_dudt_boundary_NS(dY,Boundary,TR)
%CLEAR_DUDT_BOUNDARY Summary of this function goes here
%   Detailed explanation goes here
nb=length(TR.points);
for i=1:nb
    if TR.points(i,1)==Boundary(1) | TR.points(i,2)==Boundary(2) | TR.points(i,2)==Boundary(3)
        dY(i)=0;
        dY(nb+i)=0;
%         dY(2*nb+i)=0;
    elseif norm(TR.points(i,:)-[0.15 0.15])<0.02
        dY(i)=0;
        dY(nb+i)=0;
        dY(2*nb+i)=0;
    end
end

end

