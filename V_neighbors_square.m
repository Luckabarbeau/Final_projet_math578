% code by Lucka Barbeau
function neighbors=V_neighbors_square(TR,k)
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
neighbors_idex=neighbors;
neighbors_p=TR.points(neighbors,:);

P=TR.points(k,:);
Candidat_Y_P=find(neighbors_p(:,1)==P(1) & neighbors_p(:,2)>P(2));
Candidat_X_P=find(neighbors_p(:,2)==P(2) & neighbors_p(:,1)>P(1));
Candidat_Y_M=find(neighbors_p(:,1)==P(1) & neighbors_p(:,2)<P(2));
Candidat_X_M=find(neighbors_p(:,2)==P(2) & neighbors_p(:,1)<P(1));
neighbors=[];
if isempty(Candidat_X_P)==0
[n,idex_XP]=min(neighbors_p(Candidat_X_P,1)-P(1));
neighbors=[neighbors neighbors_idex(Candidat_X_P(idex_XP))];
end

if isempty(Candidat_X_M)==0
[n,idex_XM]=min(-neighbors_p(Candidat_X_M,1)+P(1));
neighbors=[neighbors neighbors_idex(Candidat_X_M(idex_XM))];
end

if isempty(Candidat_Y_P)==0
[n,idex_YP]=min(neighbors_p(Candidat_Y_P,2)-P(2));
neighbors=[neighbors neighbors_idex(Candidat_Y_P(idex_YP))];
end

if isempty(Candidat_Y_M)==0
[n,idex_YM]=min(-neighbors_p(Candidat_Y_M,2)+P(2));
neighbors=[neighbors neighbors_idex(Candidat_Y_M(idex_YM))];
end




end