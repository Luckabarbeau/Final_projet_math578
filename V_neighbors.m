% code by Lucka Barbeau
function neighbors=V_neighbors(TR,k)
%%% find the vertex neighbors of a specific vertex k in the trianulation TR
neighbors_i=mod(find(TR.connect_active==k),length(TR.connect_active));
neighbors_i(neighbors_i==0)=length(TR.connect_active);


neighbors_value=zeros(length(neighbors_i),4);
for i=1:length(neighbors_i)
    neighbors_value(i,:)=TR.connect_active(neighbors_i(i),:);
end
neighbors=[];
for i=1:length(neighbors_i)*3
    if neighbors_value(i)~=k & isempty(find(neighbors==neighbors_value(i)))
    neighbors=[neighbors neighbors_value(i)];
    end
end

end