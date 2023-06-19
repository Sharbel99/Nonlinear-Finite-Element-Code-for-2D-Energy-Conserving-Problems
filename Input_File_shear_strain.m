% Input File: One Quad Elements Under Uniaxial Load for Nonlinear Small 
% Strain Elasticity
% Copyright (C) Arif Masud, Tim Truster, and Shoaib Goraya.

% This program computes a numerical solution to a finite element model
% using input on the geometry and physical properties of a mesh, and on the
% boundary conditions and applied loads. The routine assembles element
% quantities into the stiffness matrix and force vector to create a system
% of equations which is then solved for the nodal values of the
% displacement field. Boundary conditions are applied to constrain the
% stiffness matrix and to augment the force vector. The output is a list of
% the displacements printed on screen, contour plots of the displacement
% fields, and a plot of the deformed configuration of the mesh.
%
% Mesh input should be uploaded by running an input .m file before
% executing this program.
%
% Format of required input:
%
%   numnp:           = number of nodes in the mesh (length(NodeTable))
%
%   numel:           = number of elements in the mesh
%
%   nen:             = maximum number of nodes per element (4)
%
%   PSPS:            = flag for plane stress ('s') or plane strain ('n')
%
%   NodeTable:       = table of mesh nodal coordinates defining the
%                      geometry of the mesh; format of the table is as
%                      follows:
%                          Nodes  |             x-coord  y-coord
%                          n1     |  NodeTable = [x1     y1
%                          n2     |               x2     y2
%                          ...    |               ..     ..
%                          nnumnp |               xnumnp ynumnp];
%
%   ix:              = table of mesh connectivity information, specifying
%                      how nodes are attached to elements and how materials
%                      are assigned to elements; entries in the first nen
%                      columns correspond to the rows of NodeTable
%                      representing the nodes attached to element e;
%                      entries in the last nen+1 column are rows from MateT
%                      signifying the material properties assigned to
%                      element e; format of the table is as follows:
%                          Elements  |         n1    n2    n3    n4   mat
%                          e1        |  ix = [e1n1  e1n2  e1n3  e1n4 e1mat
%                          e2        |        e2n1  e2n2  e2n3  e2n4 e2mat
%                          ...       |         ..    ..    ..    ..   ..
%                          enumel    |        values for element numel   ];
%
%   MateT:           = table of mesh material properties for each distinct
%                      set of material properties; these sets are
%                      referenced by element e by setting the value of
%                      ix(e,nen+1) to the row number of the desired
%                      material set; format of the table is as follows:
%                          Materials  |           E   v   t
%                          mat1       |  MateT = [E1  v1  t1
%                          mat2       |           E2  v2  t2
%                          ...        |           ..  ..  ..];
%
%   BCLIndex:        = list of the number of boundary conditions and loads
%                      applied to the mesh; first entry is the number of
%                      prescribed displacements at nodes; second entry is
%                      the number of nodal forces
%
%   NodeBC:          = table of prescribed nodal displacement boundary
%                      conditions; it contains lists of nodes, the
%                      direction of the displacement prescribed (x=1, y=2),
%                      and the value of the displacement (set 0 for fixed
%                      boundary); the length of the table must match the
%                      entry in BCLIndex(1), otherwise an error will result
%                      if too few conditions are given or extra BCs will be
%                      ignored in the model input module;  format of the 
%                      table is as follows:
%                          BCs  |            nodes direction value
%                          bc1  |  NodeBC = [bc1n   bc1dir   bc1u
%                          bc2  |            bc2n   bc2dir   bc2u
%                          ...  |             ..     ..       .. ];
%
%   NodeLoad:        = table of prescribed nodal forces; it contains lists 
%                      of nodes, the direction of the force prescribed 
%                      (x=1, y=2), and the value of the force; the length 
%                      of the table must match the entry in BCLIndex(2), 
%                      otherwise an error will result if too few conditions
%                      are given or extra loads will be ignored in the 
%                      model input module; format of the table is as
%                      follows:
%                          Loads  |              nodes direction value
%                          P1     |  NodeLoad = [ P1n    P1dir    P1P
%                          P2     |               P2n    P2dir    P2P
%                          ...    |               ..     ..       .. ];
%
% The following numbering convention is used for 4-node quadrilateral
% elements:
%
%           4 -------------- 3       2
%           |                |       | \
%           |                |       |  \
%           |                |       |   \
%           |                |       |    \
%           |                |       |     \
%           |                |       |      \
%           1 -------------- 2       3-------1
%

% IMPORTANT NOTE:
% This code is for the strain energy functional given below (e=epsilon)
% W(e) = 0.5*a*e_mm*e_nn + 0.25*b*(e_mn*e_mn + e_mn*e_mn) + 1/3 * c *
% (0.75*e_mn*e_np*e_mp + 0.25*e_mn*e_np*e_pm)
% You need to change D matrix and stress tensor for an energy functional
% different than the above. Rest of the code strucure remains same.

clear; clc; close all;

% Arbitrary data for assistance in defining the mesh
L = 1;
H = 1;

% Mesh Nodal Coordinates
NodeTable = [0 0
             L 0
             L H
             0 H];

numnp = length(NodeTable); % total nodes

% Mesh Element Connectivities
ix = [1 2 3 4 1];
nen = 4;
numel = 1;
ndf = 2;
iel = 1;
ndm = 2;

% Mesh Boundary Conditions 
 BCLIndex = [4 0]';
 NodeBC =[1 1 0
           1 2 0
           2 1 0
           2 2 0];
 P = 0;       
 NodeLoad = [3 1 P
              4 1 P];

%Material Parameters
E=100;
nu=0.25;
rho = 1;
thick=1;
mu = E/2/(1+nu);
lambda = E*nu/(1+nu)/(1-2*nu);
MateT = [mu lambda thick, rho];

%Initialization
strain=zeros(4,3);
stress=zeros(4,3);
% ModelDx=zeros(numnp*ndf-BCLIndex(1),1); %unknown displacements
NRmax=5; %max no. of iterations permitted in Newton Raphsion Algorithm
tol=1e-10; %residual tolerance

alpha = -1/3;
gamma = (1-2*alpha)/2; 
beta = (1-alpha)^2/4; 

%Define Time Parameters
dt=0.01; %Time step increment
tmax=10; %End Time

%Initializing
assign_bc_load_data
% ModelDx = [ 0.2; 0; 0.2; 0];
% ModelDx = [0.185984905539714; -0.0887170637599883;0.200020996088810;0.0759712332773637];
ModelDx = [0.200047672616614; -0.0961786891685800; 0.216377545726565; 0.0813547903685381];
ModelVx = zeros(size(ModelDx));
ModelAx = ModelVx;
Strain_energy = 0;
FormFE
M_star = (1+alpha)*beta*dt^2*Kdd+Mdd;
Fint_old = Fint;
Fstar = Fd-((1+alpha)*Fint-alpha*Fint_old)-Mdd*ModelAx; %since Fext does not change
% ModelAx = Mdd\(Fd-Fint);  %% This somehow adds damping

%Run FEA Program
FEA_Program

% Plots

% Plot dx
figure(2)
plot(time,d(1,:))
grid on
title('X Displacement at Node 3','fontsize',15)
axis([0 10 -0.25 0.25])
set(gca,'fontsize',13)
xlabel('Time','fontsize',13)

% Plot dy
figure(3)
plot(time,d(2,:))
grid on
title('Y Displacement at Node 3','fontsize',15)
axis([0 10 -0.2 0.2])
% axis([0 10 -0.15 0.15])
set(gca,'fontsize',13)
xlabel('Time','fontsize',13)

% Plot vx
figure(4)
plot(time,v(1,:))
grid on
title('X Velocity at Node 3','fontsize',15)
axis([0 10 -3.5 3.5])
% axis([0 10 -2 2])
set(gca,'fontsize',13)
xlabel('Time','fontsize',13)

% Plot vy
figure(5)
plot(time,v(2,:))
grid on
title('Y Velocity at Node 3','fontsize',15)
axis([0 10 -3.5 3.5])
% axis([0 10 -1.5 1.5])
set(gca,'fontsize',13)
xlabel('Time','fontsize',13)

%Plot Stress 
figure(6)
plot(time,STRESS(:,1))
grid on
title ('\sigma_{11} at integration point #3','fontsize',15)
axis([0 10 -15 15])
% axis([0 10 -4 4])
set(gca,'fontsize',13)
xlabel('Time','fontsize',13)


figure(7)
plot(time,STRESS(:,2))
grid on
title ('\sigma_{22} at integration point #3','fontsize',15)
axis([0 10 -15 15])
% axis([0 10 -10 10])
set(gca,'fontsize',13)
xlabel('Time','fontsize',13)


figure(8)
plot(time,STRESS(:,3))
grid on
title ('\sigma_{12} at integration point #3','fontsize',15)
axis([0 10 -8 8])
% axis([0 10 -3 3])
set(gca,'fontsize',13)
xlabel('Time','fontsize',13)


