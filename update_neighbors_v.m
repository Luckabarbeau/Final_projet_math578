function TR = update_neighbors_v(TR)
%UPDATE_NEIGHBORS Summary of this function goes here
% redefine the vertex neigbors

parfor i= 1: length(TR.points)
    neighbors_v{i,1}=V_neighbors_square(TR,i);
end
TR.neighbors_v=neighbors_v;
% define edge 
TR.edge=[];
for i=1:length(TR.points)
    neighbors=TR.neighbors_v{i,1};
    for j=1:length(neighbors)
        if length(TR.edge)>0
        if isempty(find((TR.edge(:,1)==neighbors(j) & TR.edge(:,2)==i) | (TR.edge(:,2)==neighbors(j) & TR.edge(:,1)==i)))
        TR.edge=[TR.edge ; i  neighbors(j)];
        end
        else
         TR.edge=[TR.edge ; i  neighbors(j)];   
        end
    end
end

end

