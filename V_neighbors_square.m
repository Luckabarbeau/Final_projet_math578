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
px=TR.points(k,1);
py=TR.points(k,2);
for i=1:length(neighbors_i)*4
    if neighbors_value(i)~=k & isempty(find(neighbors==neighbors_value(i))) & (TR.points(neighbors_value(i),1)==px | TR.points(neighbors_value(i),2)==py)
    neighbors=[neighbors neighbors_value(i)];
    end
end

end