function B_X = immersed_boundary_value_function(X,type)
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

% if type==1
% B_X=(X(2)-0.25)/0.02*1;
% elseif type==2
% B_X=(X(1)-0.2)/0.02*-1;
% end

 B_X=0;
end

