% code by Lucka Barbeau
function dx_stancil=get_stancil_dx(TR,k,order)
%give dx stancil for the point k in triangulation TR

neighbors=TR.neighbors_v{k,1};
dx_stancil=zeros(1,length(TR.points));
N=length(neighbors);
dx_total=0;
if order==3
    Nb_neighbors=length(neighbors);
    for i=1:Nb_neighbors
        P_neighbors=TR.neighbors_v{neighbors(i),1};
        for j=1:length(P_neighbors)
            if isempty(find(neighbors==P_neighbors(j))) & P_neighbors(j)~=k
                neighbors=[neighbors P_neighbors(j)];
            end
        end
    end
    dx_candidate=[];
    N=length(neighbors);
    for i=1:N
        dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
        if dx_dy(1)~=0 & dx_dy(2)==0
            dx_total=abs(dx_dy(1))+dx_total;
            dx_candidate=[dx_candidate; dx_dy(1) neighbors(i)  ];
        end
    end
    dx_candidate=sortrows(dx_candidate,'descend');

elseif order==2
    for i=1:N
        dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
        if dx_dy(1)~=0 & dx_dy(2)==0
            dx_total=abs(dx_dy(1))+dx_total;
        end
    end
end


if order==2
for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    dx=dx_dy(1);
    r=norm(dx_dy);
    if dx~=0 
    dx_stancil(k)=dx_stancil(k)+-1/dx*abs(dx)/dx_total;
    dx_stancil(neighbors(i))=dx_stancil(neighbors(i))+1/dx*abs(dx)/dx_total;
    end
end

elseif order==3
    if length(dx_candidate)==4
            dx_stancil(dx_candidate(1,2))=-1/12*1/abs(dx_candidate(2,1));
            dx_stancil(dx_candidate(2,2))=2/3*1/abs(dx_candidate(2,1));
            dx_stancil(dx_candidate(3,2))=-2/3*1/abs(dx_candidate(2,1));
            dx_stancil(dx_candidate(4,2))=1/12*1/abs(dx_candidate(2,1));
    elseif length(dx_candidate)==3
        if dx_candidate(1,1)>0 & dx_candidate(2,1)>0
            dx_stancil(dx_candidate(2,2))=1/2*1/abs(dx_candidate(2,1));
            dx_stancil(dx_candidate(3,2))=-1/2*1/abs(dx_candidate(2,1)); 
        else
            dx_stancil(dx_candidate(1,2))=1/2*1/abs(dx_candidate(2,1));
            dx_stancil(dx_candidate(2,2))=-1/2*1/abs(dx_candidate(2,1)); 
        end
    elseif length(dx_candidate)==2
        if dx_candidate(1,1)>0
            dx_stancil(k)=-3/2*1/abs(dx_candidate(2,1));
            dx_stancil(dx_candidate(2,2))=2*1/abs(dx_candidate(2,1));
            dx_stancil(dx_candidate(1,2))=-1/2*1/abs(dx_candidate(2,1));
        else
            dx_stancil(k)=3/2*1/abs(dx_candidate(2,1));
            dx_stancil(dx_candidate(1,2))=-2*1/abs(dx_candidate(2,1));
            dx_stancil(dx_candidate(2,2))=1/2*1/abs(dx_candidate(2,1));
        end
    end

elseif order==1
    for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    dx=dx_dy(1);
    r=norm(dx_dy);
    if dx>0
    dx_stancil(k)=dx_stancil(k)+-1/dx;
    dx_stancil(neighbors(i))=dx_stancil(neighbors(i))+1/dx;
    break
    elseif dx<0
    dx_stancil(k)=dx_stancil(k)+-1/dx;
    dx_stancil(neighbors(i))=dx_stancil(neighbors(i))+1/dx;
    break
    end
end
end
end