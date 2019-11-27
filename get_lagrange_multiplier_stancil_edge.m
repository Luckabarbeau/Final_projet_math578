function LM_stancil = get_lagrange_multiplier_stancil_edge(TR,edge,index)
%GET_LAGRANGE_MULTIPLIER_STANCIL Summary of this function goes here
% return the stancil for the lagrange multiplier at point X

LM_stancil=zeros(1,length(TR.points));

d=TR.points(TR.edge(edge.edge_index(index),2),:)-TR.points(TR.edge(edge.edge_index(index),1),:);
del=norm(d);

if d(1)>0 
    del_p=norm(edge.position(index,:)-TR.points(TR.edge(edge.edge_index(index),1),:));
    LM_stancil(TR.edge(edge.edge_index(index),2))=1/del*del_p;
    LM_stancil(TR.edge(edge.edge_index(index),1))=1-1/del*del_p;
elseif d(1)<0 
    del_p=norm(edge.position(index,:)-TR.points(TR.edge(edge.edge_index(index),2),:));
    LM_stancil(TR.edge(edge.edge_index(index),1))=1/del*del_p;
    LM_stancil(TR.edge(edge.edge_index(index),2))=1-1/del*del_p;
elseif d(2)>0
    del_p=norm(edge.position(index,:)-TR.points(TR.edge(edge.edge_index(index),1),:));
    LM_stancil(TR.edge(edge.edge_index(index),2))=1/del*del_p;
    LM_stancil(TR.edge(edge.edge_index(index),1))=1-1/del*del_p;
elseif d(2)<0
    del_p=norm(edge.position(index,:)-TR.points(TR.edge(edge.edge_index(index),2),:));
    LM_stancil(TR.edge(edge.edge_index(index),1))=1/del*del_p;
    LM_stancil(TR.edge(edge.edge_index(index),2))=1-1/del*del_p;
end





end