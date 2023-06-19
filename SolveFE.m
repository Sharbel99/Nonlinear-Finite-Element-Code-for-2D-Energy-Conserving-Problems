% Solve Partitioned Finite Element Matrix System
%
% Copyright (C) Arif Masud and Tim Truster
% 7/2009
% UIUC

%Move Constrained DOF to RHS

Fdtilda = zeros(neq,1);
% to accomodate the change of the reference, these modifications were
% necessary
% if step==1
    for i = 1:neq
    
        rhs = 0;
        for j = 1:nieq
            rhs = rhs + Kdf(i,j)*ModelDc(j);
        end
        Fdtilda(i) = Fd(i) - rhs;
    
    end
% else
%     Fdtilda = Fd;
% end

%Solve Kd = F
%Calulate Residual

M_star = (1+alpha)*beta*dt^2*Kdd+Mdd;
R = Fdtilda-((1+alpha)*Fint-alpha*Fint_old)-Mdd*ModelAx;
delAx=M_star\R;

ModelAx = ModelAx + delAx;
ModelVx = ModelVx + dt*gamma*delAx;
ModelDx = ModelDx + dt^2*beta*delAx;
