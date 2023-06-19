function [ElemK, ElemF, ElemFint,stress,strain,ElemM,S_energy] = Elast2d_Elem(xl,ul,mateprop,nel,ndf,ndm)
%
% Copyright (C) Arif Masud and Tim Truster
%
% Subroutine to compute stiffness matrix and force vector for linear
% 2-dimensional elasticity element. Element currently supports bilinear
% quadrilateral elements with the following node and shape function
% labelling scheme:
%
%  (-1, 1)  4 -------------- 3 ( 1, 1)
%           |       s        |
%           |       ^        |
%           |       |        |
%           |       .-> r    |
%           |                |
%           |                |
%  (-1,-1)  1 -------------- 2 ( 1,-1)
%
% Element local coordinates (r,s) are defined by a coordinate axis with the
% origin at the center of the element; the corners of the element have
% local coordinate values as shown in the figure.
%
% Definitions for input:
%
%   xl:              = local array containing (x,y) coordinates of nodes
%                      forming the element; format is as follows:
%                          Nodes    |        n1  n2  n3  n4
%                          x-coord  |  xl = [x1  x2  x3  x4
%                          y-coord  |        y1  y2  y3  y4];
%
%   mateprop:        = vector of material properties:
%                          mateprop = [E v t]; 
%                                   = [(Young's Modulus) (Poisson's Ratio)
%                                      (thickness)];
%
%   nel:             = number of nodes on current element (4)
%
%   ndf:             = max number of DOF per node (2)
%
%   ndm:             = space dimension of mesh (2)
%
%   PSPS:            = flag for plane stress ('s') or plane strain ('n')
%
% Definitions for output:
%
%   ElemK:           = element stiffness matrix containing stiffness
%                      entries in the following arrangement, where
%                      wij corresponds to weighting function (i), coordinate
%                      direction (j), and ukl corresponds to displacement
%                      function (k), coordinate direction (l):
%                                 u1x  u1y  u2x  u2y  u3x  u3y  u4x  u4y
%                      w1x  ElemK[ .    .    .    .    .    .    .    .
%                      w1y         .    .    .    .    .    .    .    .
%                      w2x         .    .    .    .    .    .    .    .
%                      w2y         .    .    .    .    .    .    .    .
%                      w3x         .    .    .    .    .    .    .    .
%                      w3y         .    .    .    .    .    .    .    .
%                      w4x         .    .    .    .    .    .    .    .
%                      w4y         .    .    .    .    .    .    .    . ];
%
%   ElemF:           = element force vector containing force entries in the
%                      following arrangement:
%                      w1x  ElemF[ . 
%                      w1y         . 
%                      w2x         . 
%                      w2y         . 
%                      w3x         . 
%                      w3y         . 
%                      w4x         . 
%                      w4y         . ];     
%
%   ElemFint:        = element internal force vector containing force 
%                      entries in the following arrangement:
%                      w1x  ElemFint[ . 
%                      w1y            . 
%                      w2x            . 
%                      w2y            . 
%                      w3x            . 
%                      w3y            . 
%                      w4x            . 
%                      w4y            . ];   
%
% Definitions of local constants:
%
%   nst:             = size of element arrays (ndf*nel)
%
%

% Set Material Properties

mu = mateprop(1);
lambda = mateprop(2);
thick = mateprop(3);
rho = mateprop(4);
ul_elem = reshape(ul,ndf*nel,1);

% Initialize Matrix and Vector

nst = nel*ndf;
ElemK = zeros(nst);
ElemF = zeros(nst,1);
ElemFint = zeros(nst,1);
ElemM = zeros(nst);
S_energy = 0;

% Load Guass Integration Points

if nel == 3
    lint = 4;
else
    lint = 4;
end

delta=eye(2);

% Loop over integration points
for l = 1:lint

        if nel == 3
            [Wgt,r,s] =  intpntt(l,lint,0);
        else
            [Wgt,r,s] =  intpntq(l,lint,0);
        end

        % Evaluate local basis functions at integration point
        shp = shpl_2d(r,s,nel);

        % Evaluate first derivatives of basis functions at int. point
        [Qxy_ref, Jdet_ref] = shpg_2d(shp,xl,nel);
        [Qxy, Jdet] = shpg_2d(shp,xl+ul,nel);

        % Form B matrix
        if nel == 3
        Bmat = [Qxy(1,1) 0        Qxy(1,2) 0        Qxy(1,3) 0        
                0        Qxy(2,1) 0        Qxy(2,2) 0        Qxy(2,3)
                Qxy(2,1) Qxy(1,1) Qxy(2,2) Qxy(1,2) Qxy(2,3) Qxy(1,3)];
        else
        Bmat = [Qxy(1,1) 0        Qxy(1,2) 0        Qxy(1,3) 0        Qxy(1,4) 0 
                0        Qxy(2,1) 0        Qxy(2,2) 0        Qxy(2,3) 0        Qxy(2,4)
                Qxy(2,1) Qxy(1,1) Qxy(2,2) Qxy(1,2) Qxy(2,3) Qxy(1,3) Qxy(2,4) Qxy(1,4)];
        Nmat = [shp(3,1) 0        shp(3,2) 0        shp(3,3) 0        shp(3,4) 0        
                0        shp(3,1) 0        shp(3,2) 0        shp(3,3) 0        shp(3,4)];
        end
        
        Gradu=[ Qxy_ref(1,:)*ul(1,:)' Qxy_ref(2,:)*ul(1,:)'
                Qxy_ref(1,:)*ul(2,:)' Qxy_ref(2,:)*ul(2,:)'];
        Fi = eye(2)+Gradu;
        F = Fi;
        J = det(F);
        C = F'*F;
        Bt = F*F';
        trC = C(1,1) + C(2,2);
        
        for i=1:2
            for j=1:2
                for k=1:2
                    for l1=1:2
                        c(i,j,k,l1) = lambda/J*delta(i,j)*delta(k,l1)...
                            + (mu/J-lambda*log(J)/J)*...
                            (delta(i,k)*delta(j,l1)+delta(i,l1)*delta(j,k));
                    end
                end
            end
        end
        
        stress_matrix = mu/J*Bt +(-mu+lambda*log(J))/J*eye(2);
        
%         ElemKi=zeros(ndf*nel); 
%        %this inefficient way can guarantee the accuracy
%        for a=1:nel
%            for b=1:nel
%                for i=1:ndf
%                     for j=1:ndf
%                         for k=1:ndf
%                           for l1=1:ndf
%                               p=ndf*(a-1)+i;
%                               q=ndf*(b-1)+k;
%                               ElemKi(p,q)= ElemKi(p,q)+Qxy(j,a)*(c(i,j,k,l1)...
%                                   +stress_matrix(j,l1)*delta(i,k))*Qxy(l1,b);
%                           end
%                         end
%                     end
%                 end
%             end
%        end  

       %efficiently, tangent can be computed like this

       temp = (mu/J-lambda*log(J)/J);

       Dmat = [lambda/J+2*temp + stress_matrix(1,1), lambda/J, stress_matrix(1,2)
               lambda/J, lambda/J+2*temp + stress_matrix(2,2), 0
               stress_matrix(1,2),0,temp + stress_matrix(2,2)];
    
            %get strain
        strain_matrix = 0.5*(C-eye(2));
        
        strain(l,1) = strain_matrix(1,1);
        strain(l,2) = strain_matrix(2,2);
        strain(l,3) = strain_matrix(1,2);
        
        stress(l,1) = stress_matrix(1,1);
        stress(l,2) = stress_matrix(2,2);
        stress(l,3) = stress_matrix(1,2);%because it is symmetric

        %Get the Lumped Mass part
        ElemM1=(rho)/J*Nmat'*Nmat;
        ElemM2=zeros(length(ElemM1));
        for i=1:length(ElemM1)
            for j=1:length(ElemM1)
                ElemM2(i,i)=ElemM2(i,i)+ElemM1(i,j);
            end
        end


        
        % Update integration weighting factor
        W = Wgt*Jdet*thick;

        % Calculate the strain energy
        S_energy = S_energy + W*(mu/2*(trC-2)-mu*log(J)+lambda/2*log(J)^2); %%%%%%%%% 3 replaced by 2 since W(F) needs to be 0 for F=I
                             
        ElemM = ElemM + W*ElemM1;

%         ElemK = ElemK + W*ElemKi; %inefficient tangent
        ElemK = ElemK + W*Bmat'*Dmat*Bmat; %efficient tangent
        ElemFint = ElemFint + W*Bmat'*stress(l,:)'; %internal force

end %je