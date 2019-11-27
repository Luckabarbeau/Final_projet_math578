% code by Lucka Barbeau
function cells=cells_define(Xmin,Xmax,Ymin,Ymax,initial_refinement)
    
       cells.points(1,:)=[Xmin,Ymin];
       cells.points(2,:)=[Xmax,Ymin];
       cells.points(3,:)=[Xmax,Ymax];
       cells.points(4,:)=[Xmin,Ymax];
       cells.connect(1,:)=[1 2 3 4];
       cells.order(1)=1;
       cells.active(1)=true;
       cells.parent(1)=0;
       
       cells.hanging(1)=false;
       cells.hanging(2)=false;
       cells.hanging(3)=false;
       cells.hanging(4)=false;
       
       cells.boundary(1,:)=[Xmin Xmax Ymin Ymax];
       cells.area(1,:)=(Xmax-Xmin)*(Ymax-Ymin);
       cells.center(1,:)=[(Xmin+Xmax)/2 (Ymin+Ymax)/2];
       cells.immersed_iteration(1)=0;

    for i=1:initial_refinement
        active_cells_index=find(cells.active==1);
        nb_cells=length(active_cells_index);
        disp('initial refinement order')
        disp(i)
        for k=1:nb_cells
            j=active_cells_index(k);
            cells=raffine_cell(cells,j);

        end
    end
    plot(cells.points(:,1),cells.points(:,2),'O')
    active=find(cells.active==1);
    for i=1:length(active)
    cells.connect_active(i,:)=cells.connect(active(i),:);
    end
    cells.active_list=active;
end