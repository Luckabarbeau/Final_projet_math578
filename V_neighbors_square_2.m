% code by Lucka Barbeau
function [neighbors dist]=V_neighbors_2(TR,k)
%%% find the vertex neighbors of a specific vertex k in the trianulation TR
P=TR.points(k,:);
Candidat_Y_P=find(TR.points(:,1)==P(1) & TR.points(:,2)>P(2));
Candidat_X_P=find(TR.points(:,2)==P(2) & TR.points(:,1)>P(1));
Candidat_Y_M=find(TR.points(:,1)==P(1) & TR.points(:,2)<P(2));
Candidat_X_M=find(TR.points(:,2)==P(2) & TR.points(:,1)<P(1));
neighbors=[];
if isempty(Candidat_X_P)==0
[n,idex_XP]=min(TR.points(Candidat_X_P,1)-P(1));
neighbors=[neighbors Candidat_X_P(idex_XP)];
end

if isempty(Candidat_X_M)==0
[n,idex_XM]=min(-TR.points(Candidat_X_M,1)+P(1));
neighbors=[neighbors Candidat_X_M(idex_XM)];
end

if isempty(Candidat_Y_P)==0
[n,idex_YP]=min(TR.points(Candidat_Y_P,2)-P(2));
neighbors=[neighbors Candidat_Y_P(idex_YP)];
end

if isempty(Candidat_Y_M)==0
[n,idex_YM]=min(-TR.points(Candidat_Y_M,2)+P(2));
neighbors=[neighbors Candidat_Y_M(idex_YM)]
end

end