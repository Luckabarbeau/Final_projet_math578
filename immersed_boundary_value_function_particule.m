function B_X = immersed_boundary_value_function_particule(X,type,v,w)
%IMMERSED_BOUNDARY_VALUE_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

% B_X=norm(X);
% B_X=exp(norm(X));
% B_X=1;
% if norm(X-[0.5 0.5])<0.15
%     B_X=1;
% else
%     B_X=0;
% end

if type==1
B_X=v(1);
elseif type==2
B_X=v(2);
end

%  B_X=0;
end
