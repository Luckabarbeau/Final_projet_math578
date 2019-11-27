function dY = clear_dudt_boundary(dY,Boundary,TR)
%CLEAR_DUDT_BOUNDARY Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(TR.points)
    if TR.points(i,1)==Boundary(1) | TR.points(i,1)==Boundary(2) | TR.points(i,2)==Boundary(3)| TR.points(i,2)==Boundary(4)
        dY(i)=0;
    end
end

end

