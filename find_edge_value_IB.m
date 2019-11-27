function edge_value=find_edge_value_IB(TR,X)
edge_value.edge_index=[];
edge_value.position=[];
X1=X(1:length(X)-1,:);
X2=X(2:length(X),:);
diff_space=X1-X2;
epsilone=mean(sqrt(diff_space(:,1).^2+diff_space(:,2).^2))*10;
for i=1:length(TR.edge)
    if TR.points(TR.edge(i,1),1)-TR.points(TR.edge(i,2),1)==0
        %edge oriente en y
        if TR.points(TR.edge(i,1),2)>TR.points(TR.edge(i,2),2)
        
        X_crossing_candidate=find(X(:,2)>=TR.points(TR.edge(i,2),2)& X(:,2)<=TR.points(TR.edge(i,1),2) & X(:,1)<TR.points(TR.edge(i,1),1)+epsilone & X(:,1)>TR.points(TR.edge(i,1),1)-epsilone);
        X_crossing_candidate=X(X_crossing_candidate,:); 
        else
        X_crossing_candidate=find(X(:,2)>=TR.points(TR.edge(i,1),2)& X(:,2)<=TR.points(TR.edge(i,2),2) & X(:,1)<TR.points(TR.edge(i,1),1)+epsilone & X(:,1)>TR.points(TR.edge(i,1),1)-epsilone);
        X_crossing_candidate=X(X_crossing_candidate,:) ;
        end
        
        for j=1:size(X_crossing_candidate,1)-1
            if X_crossing_candidate(j,1)>=TR.points(TR.edge(i,1),1) & X_crossing_candidate(j+1,1)<= TR.points(TR.edge(i,1),1)
                %Crossing !!!
                edge_value.edge_index=[edge_value.edge_index ; i];
                P=[TR.points(TR.edge(i,1),1) (TR.points(TR.edge(i,1),1)-X_crossing_candidate(j+1,1))/(X_crossing_candidate(j,1)-X_crossing_candidate(j+1,1))*(X_crossing_candidate(j,2)-X_crossing_candidate(j+1,2))+X_crossing_candidate(j,2) ];
                edge_value.position=[edge_value.position ; P ];
                break
            elseif X_crossing_candidate(j,1)<=TR.points(TR.edge(i,1),1) & X_crossing_candidate(j+1,1)>= TR.points(TR.edge(i,1),1)
                %Crossing !!!
                edge_value.edge_index=[edge_value.edge_index ; i];
                P=[TR.points(TR.edge(i,1),1) (TR.points(TR.edge(i,1),1)-X_crossing_candidate(j,1))/(X_crossing_candidate(j+1,1)-X_crossing_candidate(j,1))*(X_crossing_candidate(j+1,2)-X_crossing_candidate(j,2))+X_crossing_candidate(j,2) ];
                edge_value.position=[edge_value.position ; P ];
                break
            end
        end
        
       
    else
        %edge oriente en x
        if TR.points(TR.edge(i,1),1)>TR.points(TR.edge(i,2),1)
        X_crossing_candidate=find(X(:,1)>=TR.points(TR.edge(i,2),1)& X(:,1)<=TR.points(TR.edge(i,1),1) & X(:,2)<TR.points(TR.edge(i,1),2)+epsilone & X(:,2)>TR.points(TR.edge(i,1),2)-epsilone);
        X_crossing_candidate=X(X_crossing_candidate,:) ;
        else
        X_crossing_candidate=find(X(:,1)>=TR.points(TR.edge(i,1),1)& X(:,1)<=TR.points(TR.edge(i,2),1) & X(:,2)<TR.points(TR.edge(i,1),2)+epsilone & X(:,2)>TR.points(TR.edge(i,1),2)-epsilone);
        X_crossing_candidate=X(X_crossing_candidate,:) ;
        end
        
        for j=1:size(X_crossing_candidate,1)-1
            
            if X_crossing_candidate(j,2)>=TR.points(TR.edge(i,1),2) & X_crossing_candidate(j+1,2)<= TR.points(TR.edge(i,1),2)
                %Crossing !!!
                
                edge_value.edge_index=[edge_value.edge_index ; i];
                P=[(TR.points(TR.edge(i,1),2)-X_crossing_candidate(j+1,2))/(X_crossing_candidate(j,2)-X_crossing_candidate(j+1,2))*(X_crossing_candidate(j,1)-X_crossing_candidate(j+1,1))+X_crossing_candidate(j,1) TR.points(TR.edge(i,1),2)];
                edge_value.position=[edge_value.position ; P ];
                break
            elseif X_crossing_candidate(j,2)<=TR.points(TR.edge(i,1),2) & X_crossing_candidate(j+1,2)>= TR.points(TR.edge(i,1),2)
                %Crossing !!!
                edge_value.edge_index=[edge_value.edge_index ; i];
                P=[(TR.points(TR.edge(i,1),2)-X_crossing_candidate(j,2))/(X_crossing_candidate(j+1,2)-X_crossing_candidate(j,2))*(X_crossing_candidate(j+1,1)-X_crossing_candidate(j,1))+X_crossing_candidate(j+1,1) TR.points(TR.edge(i,1),2)];
                edge_value.position=[edge_value.position ; P ];
                break
            end
        end
    end
end



for i =1 :length(edge_value.position)
    edge_value.value(i,:)=immersed_boundary_value_function(edge_value.position(i,:));
end



end