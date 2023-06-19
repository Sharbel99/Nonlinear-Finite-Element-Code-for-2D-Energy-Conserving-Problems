
format compact


Ro=tol;

stepcount=1;
for  t=dt:dt:tmax

    ModelDx = ModelDx + dt*ModelVx + dt^2/2*(1-2*beta)*ModelAx;
    ModelVx = ModelVx + (1-gamma)*dt*ModelAx;
    
    % Interpret Boundary Conditions and assign Loads; allocate dof's
    assign_bc_load_data
    nneq = neq + nieq;

    NR=1; R=1;
    ModelAx=zeros(4,1);
    
    % Newton Raphson Loop
    while norm(R)>tol*Ro && norm(R) >1e-15 && NR<=NRmax

        %Assemble Stiffness and Internal Force
        isw = 3;
        FormFE 
%pause(inf)
        %Solve Matrix System for FE Solution and Form Residual
        SolveFE
        
        if NR==1
            Ro=norm(R);
        end
        
        residual(NR,stepcount)=norm(R); %Storing Residual
        NR=NR+1;

    end
%     residual(NR:end,stepcount) = 0;
    
    %storing strains, stresses and displacements at int. point 3 (4 in the code)
    time(stepcount) = t;
    STRAIN(stepcount,1:3)=strain(4,1:3); 
    STRESS(stepcount,1:3)=stress(4,1:3);
    d(1:2,stepcount)=ModelDx(3:4);
    v(1:2,stepcount)=ModelVx(3:4);
%     a(1:2,stepcount)=ModelAx(3:4);
    stepcount=stepcount+1;

    Fint_old = Fint;


    %Post procesing

%     Node_U_V = zeros(numnp,ndf);
%     
%     for node = 1:numnp
%         for dir = 1:ndf
%             gDOF = NDOFT(node, dir);
%             if gDOF <= neq
%                 Node_U_V(node, dir) = ModelDx(gDOF,1);
%             else
%                 Node_U_V(node, dir) = ModelDc(gDOF - neq);
%             end
%         end
%     end
%     
%     maxuvw = zeros(ndm,1);
%     maxxyz = zeros(ndm,1);
%     for i = 1:ndm
%         maxuvw(i) = max(abs(Node_U_V(:,i)));
%         maxxyz(i) = max(NodeTable(:,i));
%     end
%     len = sqrt(maxxyz'*maxxyz);
%     perc = 5/100;
%     factor = len/max(maxuvw)*perc;
%     NodeTable2 = NodeTable;
%     for i = 1:ndm
%         NodeTable2(:,i) = NodeTable2(:,i) + Node_U_V(:,i);%*factor;
%     end
%     
%     plotModel(NodeTable2, ix, numel, nen, 1, 1, 1, 'Deformed Configuration', 'y', 'y')
%     plotModel(NodeTable, ix, numel, nen, 1, 1, 1, 'Deformed Configuration', 'n', 'n')


end