function New_U = avreage_whit_neighbors(TR,U)
New_U=zeros(length(U),1);
for i=1:length(TR.points)
    neighbors=TR.neighbors_v{i,1};
    if length(neighbors)==4| length(neighbors)==2
    for j=1:length(neighbors)
        New_U(i)=New_U(i)+U(neighbors(j));
    end
        New_U(i)=(New_U(i))/(length(neighbors));
        
    elseif length(neighbors)==3
    for j=1:length(neighbors)
        New_U(i)=New_U(i)+U(neighbors(j));
    end
        New_U(i)=(New_U(i)+U(i))/(length(neighbors)+1);    
    end
        
end


end

