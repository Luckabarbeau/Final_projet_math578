function New_sparse_index = sparse_matrix_multiply_NS(TR,sparse_index,vector)
%SPARSE_MATRIX_MULTIPLY Summary of this function goes here
%   Multipliy line of sparse matrix by a constante specific for each line
%   FOR NS so not general
New_sparse_index=[];
parfor i=1:length(sparse_index)
    sol_index=mod(sparse_index(i,1),length(TR.points))
    if sol_index==0
        sol_index=length(TR.points);
    end
    temp=[sparse_index(i,1) sparse_index(i,2) sparse_index(i,3)*vector(sol_index)];
    New_sparse_index=[New_sparse_index; temp];
    
end
    
end

