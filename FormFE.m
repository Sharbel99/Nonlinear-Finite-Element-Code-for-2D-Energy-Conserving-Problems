% Assemble Quantities from Model Routine
%
% Copyright (C) Arif Masud and Tim Truster
%
% 7/2009
% UIUC

Kdd = zeros(neq,neq);
Kdf = zeros(neq,nieq);
Kfd = zeros(nieq,neq);
Kff = zeros(nieq,nieq);
Mdd = zeros(neq,neq);
Mdf = zeros(neq,nieq);
Mfd = zeros(nieq,neq);
Mff = zeros(nieq,nieq);
Fint = zeros(neq,1); %internal force vector
PSPS='n';

AssemQuant = 'AssemStifForc'; % macro for assembling stiffness matrix and force vectors

for elem = 1:numel
    
    %Determine element size parameters
    if nen == 3
        nel = 3;
    elseif nen == 4
        if ix(elem,nen) == 0
            nel = 3;
        else
            nel = 4;
        end
    elseif nen == 6
        nel = 6;
    else
        if ix(elem,nen) == 0
            nel = 6;
        else
            nel = 9;
        end
    end
    nst = nel*ndf;
    
    %Extract patch nodal coordinates
    xl = zeros(ndm, nel);
    ElemFlag = zeros(nel, 1);
    for k = 1:nel
        node = ix(elem,k);
        ElemFlag(k) = node;
        for l = 1:ndm
            xl(l,k) = NodeTable(node,l);
        end
    end
    
    %Extract patch material properties
    ma = ix(elem,nen+1);
    mateprop = MateT(ma,:);

    %Extract patch nodal displacements
    EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);
    
    %Compute Stress Strain
     ul = zeros(ndm,nel);
    for i = 1:nel*ndf
        ndof_index = EDOFT(i);
        if(ndof_index<=neq)
            ul(i) = ModelDx(ndof_index);
        else
            ul(i) = ModelDc(ndof_index-neq);
        end
    end
    
    %[strain,stress] = CompStrainStress_Elem(xl,ul,mateprop,nel,ndf);
    

    
    %Calculate Element Stiffness and Force
    switch iel
        case 1 %Small-Deformation Isotropic Elastostatics Element
            if ndm == 3
                L_Elem1_3d
            elseif ndm == 2
                [ElemK,ElemF,ElemFint,stress,strain,ElemM,S_energy] = Elast2d_Elem(xl,ul,mateprop,nel,ndf,ndm);
            end
        case 2
            if ndm == 3

            else %ndm == 2
                L_Elem2_2d
            end
        case 3 %Stabilized Mixed Pressure-Displacement Element
            if ndm == 2
                L_Elem3_2d
            else %ndm == 3

            end
        case 4 %Implicit Error Element
            if ndm == 2
                L_Elem4_2d
            else %ndm == 3

            end
        case 5 %Stabilized Mixed Pressure-Displacement Element, Error
            if ndm == 2
                L_Elem5_2d
            else %ndm == 3

            end
    end

  % Plot strain and stress arrays
%     fprintf('\nElement\t%1.0f\n\n',elem);
%     strain
%     stress
    
    %Assemble Element contribution to Model Quantity

    run(AssemQuant)
    
end