% code by Lucka Barbeau
function cells = raffine_cell(cells,j)
%RAFFINE_CELL and do all the operation necessary for the structure of the
%cells mesh.
Xmin_cell=min([cells.points(cells.connect(j,1),1),cells.points(cells.connect(j,2),1),cells.points(cells.connect(j,3),1),cells.points(cells.connect(j,4),1)]);
Xmax_cell=max([cells.points(cells.connect(j,1),1),cells.points(cells.connect(j,2),1),cells.points(cells.connect(j,3),1),cells.points(cells.connect(j,4),1)]);
Ymin_cell=min([cells.points(cells.connect(j,1),2),cells.points(cells.connect(j,2),2),cells.points(cells.connect(j,3),2),cells.points(cells.connect(j,4),2)]);
Ymax_cell=max([cells.points(cells.connect(j,1),2),cells.points(cells.connect(j,2),2),cells.points(cells.connect(j,3),2),cells.points(cells.connect(j,4),2)]);

point_index_old(1)=find((cells.points(:,1)==Xmin_cell & cells.points(:,2)==Ymin_cell)==1);
point_index_old(2)=find((cells.points(:,1)==Xmax_cell & cells.points(:,2)==Ymin_cell)==1);
point_index_old(3)=find((cells.points(:,1)==Xmax_cell & cells.points(:,2)==Ymax_cell)==1);
point_index_old(4)=find((cells.points(:,1)==Xmin_cell & cells.points(:,2)==Ymax_cell)==1);




% new cells generation
if sum(cells.points(:,1)==(Xmin_cell+Xmax_cell)/2 & cells.points(:,2)==(Ymin_cell+Ymax_cell)/2)==0
cells.points(length(cells.points)+1,:)=[(Xmin_cell+Xmax_cell)/2 (Ymin_cell+Ymax_cell)/2];
point_index(1)=length(cells.points);
cells.hanging(length(cells.hanging)+1)=false;
else
    point_index(1)=find((cells.points(:,1)==(Xmin_cell+Xmax_cell)/2 & cells.points(:,2)==(Ymin_cell+Ymax_cell)/2)==1);
    cells.hanging(point_index(1))=false;
end

if sum(cells.points(:,1)==(Xmin_cell+Xmax_cell)/2 & cells.points(:,2)==Ymin_cell)==0
    cells.points(length(cells.points)+1,:)=[(Xmin_cell+Xmax_cell)/2 Ymin_cell];
    point_index(2)=length(cells.points);
    cells.hanging(length(cells.hanging)+1)=true;
else
    point_index(2)=find((cells.points(:,1)==(Xmin_cell+Xmax_cell)/2 & cells.points(:,2)==Ymin_cell)==1);
    cells.hanging(point_index(2))=false;
end

if sum(cells.points(:,1)==(Xmin_cell+Xmax_cell)/2 & cells.points(:,2)==Ymax_cell)==0
    cells.points(length(cells.points)+1,:)=[(Xmin_cell+Xmax_cell)/2 Ymax_cell];
    point_index(3)=length(cells.points);
    cells.hanging(length(cells.hanging)+1)=true;
else
    point_index(3)=find((cells.points(:,1)==(Xmin_cell+Xmax_cell)/2 & cells.points(:,2)==Ymax_cell)==1);
    cells.hanging(point_index(3))=false;
end

if sum(cells.points(:,1)==Xmin_cell & cells.points(:,2)==(Ymin_cell+Ymax_cell)/2)==0
    cells.points(length(cells.points)+1,:)=[Xmin_cell (Ymin_cell+Ymax_cell)/2];
    point_index(4)=length(cells.points);
    cells.hanging(length(cells.hanging)+1)=true;
else
    point_index(4)=find((cells.points(:,1)==Xmin_cell & cells.points(:,2)==(Ymin_cell+Ymax_cell)/2)==1);
    cells.hanging(point_index(4))=false;
end

if sum(cells.points(:,1)==Xmax_cell & cells.points(:,2)==(Ymin_cell+Ymax_cell)/2)==0
    cells.points(length(cells.points)+1,:)=[Xmax_cell (Ymin_cell+Ymax_cell)/2];
    point_index(5)=length(cells.points);
    cells.hanging(length(cells.hanging)+1)=true;
else
    point_index(5)=find((cells.points(:,1)==Xmax_cell & cells.points(:,2)==(Ymin_cell+Ymax_cell)/2)==1);
    cells.hanging(point_index(5))=false;
end
cells.connect(size(cells.connect,1)+1,:)=[point_index_old(1) point_index(2) point_index(1) point_index(4)];
cells.connect(size(cells.connect,1)+1,:)=[point_index(2) point_index_old(2) point_index(5) point_index(1)];
cells.connect(size(cells.connect,1)+1,:)=[point_index(1) point_index(5) point_index_old(3) point_index(3)];
cells.connect(size(cells.connect,1)+1,:)=[point_index(4) point_index(1) point_index(3) point_index_old(4)];

cells.order(length(cells.order)+1)=cells.order(j)+1;
cells.order(length(cells.order)+1)=cells.order(j)+1;
cells.order(length(cells.order)+1)=cells.order(j)+1;
cells.order(length(cells.order)+1)=cells.order(j)+1;

cells.boundary(size(cells.boundary,1)+1,:)=[Xmin_cell (Xmax_cell+Xmin_cell)/2 Ymin_cell (Ymax_cell+Ymin_cell)/2];
cells.boundary(size(cells.boundary,1)+1,:)=[(Xmax_cell+Xmin_cell)/2 Xmax_cell Ymin_cell (Ymax_cell+Ymin_cell)/2];
cells.boundary(size(cells.boundary,1)+1,:)=[(Xmax_cell+Xmin_cell)/2 Xmax_cell (Ymax_cell+Ymin_cell)/2 Ymax_cell];
cells.boundary(size(cells.boundary,1)+1,:)=[Xmin_cell (Xmax_cell+Xmin_cell)/2 (Ymax_cell+Ymin_cell)/2 Ymax_cell];

cells.center(size(cells.center,1)+1,:)=[(Xmin_cell+(Xmax_cell+Xmin_cell)/2)/2 (Ymin_cell+(Ymax_cell+Ymin_cell)/2)/2];
cells.center(size(cells.center,1)+1,:)=[((Xmax_cell+Xmin_cell)/2+Xmax_cell)/2 (Ymin_cell+(Ymax_cell+Ymin_cell)/2)/2];
cells.center(size(cells.center,1)+1,:)=[((Xmax_cell+Xmin_cell)/2+Xmax_cell)/2 ((Ymax_cell+Ymin_cell)/2+Ymax_cell)/2];
cells.center(size(cells.center,1)+1,:)=[(Xmin_cell+(Xmax_cell+Xmin_cell)/2)/2 ((Ymax_cell+Ymin_cell)/2+Ymax_cell)/2];

cells.area(size(cells.area,1)+1,:)=((Xmax_cell+Xmin_cell)/2-Xmin_cell)*((Ymax_cell+Ymin_cell)/2-Ymin_cell);
cells.area(size(cells.area,1)+1,:)=(Xmax_cell-(Xmax_cell+Xmin_cell)/2)*((Ymax_cell+Ymin_cell)/2-Ymin_cell);
cells.area(size(cells.area,1)+1,:)=(Xmax_cell-(Xmax_cell+Xmin_cell)/2)*(Ymax_cell-(Ymax_cell+Ymin_cell)/2);
cells.area(size(cells.area,1)+1,:)=((Xmax_cell+Xmin_cell)/2-Xmin_cell)*(Ymax_cell-(Ymax_cell+Ymin_cell)/2);


cells.immersed_iteration(size(cells.immersed_iteration,1)+1,1)=0;
cells.immersed_iteration(size(cells.immersed_iteration,1)+1,1)=0;
cells.immersed_iteration(size(cells.immersed_iteration,1)+1,1)=0;
cells.immersed_iteration(size(cells.immersed_iteration,1)+1,1)=0;

cells.active(j)=false;

cells.active(length(cells.active)+1)=true;
cells.active(length(cells.active)+1)=true;
cells.active(length(cells.active)+1)=true;
cells.active(length(cells.active)+1)=true;

cells.parent(length(cells.parent)+1)=j;
cells.parent(length(cells.parent)+1)=j;
cells.parent(length(cells.parent)+1)=j;
cells.parent(length(cells.parent)+1)=j;


end

