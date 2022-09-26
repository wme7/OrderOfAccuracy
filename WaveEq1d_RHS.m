function RHS = WE1d_FD_RHS(q,~)
global c r Dx periodic_x idL idR

if not(periodic_x)
    % Set Dirichlet Boundary conditions
    q(idL,1) = 0; % p(x=0)=0
    q(idR,2) = 0; % u(x=L)=0
end

% Get the right-hand-side
RHS = -(Dx*q)*[0,r*c^2;1/r,0];

end