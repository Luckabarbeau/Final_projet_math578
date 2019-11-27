% code by Lucka Barbeau
function dy_stancil=get_stancil_dy(TR,k,order)
%give dy stancil for the point k in triangulation TR

neighbors=TR.neighbors_v{k,1};

dy_stancil=zeros(1,length(TR.points));
N=length(neighbors);
dy_total=0;
if order==3
    Nb_neighbors=length(neighbors);
    for i=1:Nb_neighbors
        P_neighbors=TR.neighbors_v{neighbors(i),1};
        for j=1:length(P_neighbors)
            if isempty(find(neighbors==P_neighbors(j)))& P_neighbors(j)~=k
                neighbors=[neighbors P_neighbors(j)];
            end
        end
    end
    
    dy_candidate=[];
    N=length(neighbors);
    for i=1:N
        dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
        if dx_dy(2)~=0 & dx_dy(1)==0
            dy_total=abs(dx_dy(2))+dy_total;
            dy_candidate=[dy_candidate; dx_dy(2) neighbors(i)  ];
        end
    end
    dy_candidate=sortrows(dy_candidate,'descend');
    
    
    
elseif order==2
    dy_stancil=zeros(1,length(TR.points));
    N=length(neighbors);
    dy_total=0;
    for i=1:N
        dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
        if dx_dy(2)~=0 & dx_dy(1)==0
            dy_total=(abs(dx_dy(2)))+dy_total;
        end
    end
end



if order==2
for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    dy=dx_dy(2);
    r=norm(dx_dy);
    if dy~=0
    dy_stancil(k)=dy_stancil(k)+-1/dy*(abs(dy))/dy_total;
    dy_stancil(neighbors(i))=dy_stancil(neighbors(i))+1/dy*(abs(dy))/dy_total;
    end
end

elseif order==3
    if length(dy_candidate)==4
            dy_stancil(dy_candidate(1,2))=-1/12*1/abs(dy_candidate(2,1));
            dy_stancil(dy_candidate(2,2))=2/3*1/abs(dy_candidate(2,1));
            dy_stancil(dy_candidate(3,2))=-2/3*1/abs(dy_candidate(2,1));
            dy_stancil(dy_candidate(4,2))=1/12*1/abs(dy_candidate(2,1));
    elseif length(dy_candidate)==3
        if dy_candidate(1,1)>0 & dy_candidate(2,1)>0
            dy_stancil(dy_candidate(2,2))=1/2*1/abs(dy_candidate(2,1));
            dy_stancil(dy_candidate(3,2))=-1/2*1/abs(dy_candidate(2,1)); 
        else
            dy_stancil(dy_candidate(1,2))=1/2*1/abs(dy_candidate(2,1));
            dy_stancil(dy_candidate(2,2))=-1/2*1/abs(dy_candidate(2,1)); 
        end
    elseif length(dy_candidate)==2
        if dy_candidate(1,1)>0
            dy_stancil(k)=-3/2*1/abs(dy_candidate(2,1));
            dy_stancil(dy_candidate(2,2))=2*1/abs(dy_candidate(2,1));
            dy_stancil(dy_candidate(1,2))=-1/2*1/abs(dy_candidate(2,1));
        else
            dy_stancil(k)=3/2*1/abs(dy_candidate(2,1));
            dy_stancil(dy_candidate(1,2))=-2*1/abs(dy_candidate(2,1));
            dy_stancil(dy_candidate(2,2))=1/2*1/abs(dy_candidate(2,1));
        end
    end

elseif order==1
    for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    dy=dx_dy(2);
    r=norm(dx_dy);
    if dy>0
    dy_stancil(k)=dy_stancil(k)+-1/dy;
    dy_stancil(neighbors(i))=dy_stancil(neighbors(i))+1/dy;
    break
    elseif dy<0
    dy_stancil(k)=dy_stancil(k)+-1/dy;
    dy_stancil(neighbors(i))=dy_stancil(neighbors(i))+1/dy;
    break
    end
end
end

end