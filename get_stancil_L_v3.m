% code by Lucka Barbeau
function L_stancil=get_stancil_L(TR,k,order)
%give laplacien stancil for the point k in triangulation TR

neighbors=TR.neighbors_v{k,1};

L_stancil=zeros(1,length(TR.points));
N=length(neighbors);
dx_total=0;
dy_total=0;
m=1;
l=1;


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
dy_candidate=[];
N=length(neighbors);
for i=1:N
    dx_dy=TR.points(neighbors(i),:)-TR.points(k,:);
    if dx_dy(1)~=0 & dx_dy(2)==0
        dx_total=abs(dx_dy(1))+dx_total;
        dx_candidate=[dx_candidate; dx_dy(1) neighbors(i)  ];
    end
    if dx_dy(2)~=0 & dx_dy(1)==0
        dy_total=abs(dx_dy(2))+dy_total;
        dy_candidate=[dy_candidate; dx_dy(2) neighbors(i)  ];
    end
end
dx_candidate=sortrows(dx_candidate,'descend');
dy_candidate=sortrows(dy_candidate,'descend');
dxx_stancil=zeros(1,length(TR.points));
dyy_stancil=zeros(1,length(TR.points));

if length(dx_candidate)==4
    dxx_stancil(dx_candidate(1,2))=-1/12*1/abs(dx_candidate(2,1))^2;
    dxx_stancil(dx_candidate(2,2))=4/3*1/abs(dx_candidate(2,1))^2;
    dxx_stancil(k)=-5/2*1/abs(dx_candidate(2,1))^2;
    dxx_stancil(dx_candidate(3,2))=4/3*1/abs(dx_candidate(2,1))^2;
    dxx_stancil(dx_candidate(4,2))=-1/12*1/abs(dx_candidate(2,1))^2;
elseif length(dx_candidate)==3
    if dx_candidate(1,1)>0 & dx_candidate(2,1)>0
        dxx_stancil(dx_candidate(2,2))=1*1/abs(dx_candidate(2,1))^2;
        dxx_stancil(k)=-2*1/abs(dx_candidate(2,1))^2;
        dxx_stancil(dx_candidate(3,2))=1*1/abs(dx_candidate(2,1))^2;
    else
        dxx_stancil(dx_candidate(1,2))=1*1/abs(dx_candidate(2,1))^2;
        dxx_stancil(k)=-2*1/abs(dx_candidate(2,1))^2;
        dxx_stancil(dx_candidate(2,2))=1*1/abs(dx_candidate(2,1))^2;
    end
elseif length(dx_candidate)==2
    if dx_candidate(1,1)>0
        dxx_stancil(k)=1*1/abs(dx_candidate(2,1))^2;
        dxx_stancil(dx_candidate(2,2))=-2/abs(dx_candidate(2,1))^2;
        dxx_stancil(dx_candidate(1,2))=1*1/abs(dx_candidate(2,1))^2;
    else
        dxx_stancil(k)=1*1/abs(dx_candidate(2,1))^2;
        dxx_stancil(dx_candidate(1,2))=1*1/abs(dx_candidate(2,1))^2;
        dxx_stancil(dx_candidate(2,2))=-2*1/abs(dx_candidate(2,1))^2;
    end
end

if length(dy_candidate)==4
    dyy_stancil(dy_candidate(1,2))=-1/12*1/abs(dy_candidate(2,1))^2;
    dyy_stancil(dy_candidate(2,2))=4/3*1/abs(dy_candidate(2,1))^2;
    dyy_stancil(k)=-5/2*1/abs(dy_candidate(2,1))^2;
    dyy_stancil(dy_candidate(3,2))=4/3*1/abs(dy_candidate(2,1))^2;
    dyy_stancil(dy_candidate(4,2))=-1/12*1/abs(dy_candidate(2,1))^2;
elseif length(dy_candidate)==3
    if dy_candidate(1,1)>0 & dy_candidate(2,1)>0
        dyy_stancil(dy_candidate(2,2))=1*1/abs(dy_candidate(2,1))^2;
        dyy_stancil(k)=-2*1/abs(dy_candidate(2,1))^2;
        dyy_stancil(dy_candidate(3,2))=1*1/abs(dy_candidate(2,1))^2;
    else
        dyy_stancil(dy_candidate(1,2))=1*1/abs(dy_candidate(2,1))^2;
        dyy_stancil(k)=-2*1/abs(dy_candidate(2,1))^2;
        dyy_stancil(dy_candidate(2,2))=1*1/abs(dy_candidate(2,1))^2;
    end
elseif length(dy_candidate)==2
    if dy_candidate(1,1)>0
        dyy_stancil(k)=1*1/abs(dy_candidate(2,1))^2;
        dyy_stancil(dy_candidate(2,2))=-2/abs(dy_candidate(2,1))^2;
        dyy_stancil(dy_candidate(1,2))=1*1/abs(dy_candidate(2,1))^2;
    else
        dyy_stancil(k)=1*1/abs(dy_candidate(2,1))^2;
        dyy_stancil(dy_candidate(1,2))=1*1/abs(dy_candidate(2,1))^2;
        dyy_stancil(dy_candidate(2,2))=-2*1/abs(dy_candidate(2,1))^2;
    end
end


L_stancil=dxx_stancil+dyy_stancil;
end


