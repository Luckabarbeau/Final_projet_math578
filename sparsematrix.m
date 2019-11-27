% code by Lucka Barbeau
function Sparse_index_V=sparsematrix(A_l,i)
%creat sparse matrix index and value from line adn previous version of the
%sparse matrix
Sparse_index_V=[];
    for j =1 :length(A_l)
        if A_l(j)~=0
        Sparse_index_V=[Sparse_index_V ; i j A_l(j)];
        end
    end
    
end