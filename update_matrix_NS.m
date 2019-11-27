function A_Sparse_global = update_matrix_NS(U_x,U_y,mesh,A_Sparse_index_dx,A_Sparse_index_dy,A_Sparse_m_L,Constraint_sparse_pressure,Constraint_sparse_conservation,Constraint_sparse_matrix,nb_zeros,I)
%UPDATE_MATRIX_NS Summary of this function goes here
%   Detailed explanation goes here
A_Sparse_index_dx_iter=sparse_matrix_multiply_NS(mesh,A_Sparse_index_dx,U_x);
A_Sparse_index_dy_iter=sparse_matrix_multiply_NS(mesh,A_Sparse_index_dy,U_y);
A_dx=sparse(A_Sparse_index_dx_iter(:,1),A_Sparse_index_dx_iter(:,2),A_Sparse_index_dx_iter(:,3));
A_dy=sparse(A_Sparse_index_dy_iter(:,1),A_Sparse_index_dy_iter(:,2),A_Sparse_index_dy_iter(:,3));
A_momentum=I-(A_dx+A_dy+A_Sparse_m_L);

A_1=[A_momentum Constraint_sparse_pressure; Constraint_sparse_conservation sparse(length(mesh.points),length(mesh.points))];

A_Sparse_global=[A_1 Constraint_sparse_matrix' ; Constraint_sparse_matrix sparse(nb_zeros,nb_zeros)];


end

