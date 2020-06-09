%% Matlab Initializations

clear;
clc
format shortEng

%% Given data in the problem

thickness = 0.2;
width = 0.5;
rho = 7840;
n_e = 100;                            %Number of elements
n_n = 101;                            %Number of nodes
dof = 2;                              %Dof of elements
E_e = repmat(200e9,n_e,1);                                   %Forming the matrix of Young's Modulus
I_e = repmat(((thickness*(width^3))/12),n_e,1);              %Forming the Moment of Interia Matrix
l_e = repmat(1,n_e,1);
weight = rho*9.81*width*thickness*l_e(1);

%% Initializing and defining the required Global Vectors

Load = zeros(dof*n_n,1);              %Forming the Global Load Vector
U = zeros(dof*n_n,1);                 %Forming the Global displacement Vector
K_g = zeros(dof*n_n);                 %Initializing the Global Stiffness Matrix

%% Defining the Boundary Conditions and Point Loads acting on various nodes

Load(2:2:end) = -10000 - (weight/n_e);

boundary = [1 2 201 202];             %Defining the Boundary Conditions on various nodes

%% Defining the Global Stiffness Matrix

for i = 1:n_e   
    a = 12; b = 6*l_e(i); c = 2*(l_e(i)^2); d = 2*c;
    K_e_loc = ((E_e(i)*I_e(i))/(l_e(i)^3))*[a b -a b; b d -b c; -a -b a -b; b c -b d];    %Calculating the Local Stiffness Matrix
    eldofs = [(dof*(i-1))+1:(dof*(i-1)+2) (dof*((i+1)-1))+1:(dof*((i+1)-1)+2)]
    K_g(eldofs,eldofs) = K_g(eldofs,eldofs) + K_e_loc;                                    %Assembling the Global Stiffness Matrix 
end

%% Caculating the Global Displacement Vector and Global Load Vector using Elimination Approach

K_g_cpy = K_g;                                         %Storing a copy of Global Stiffness Matrix for pos-porcessing requirements
Load_cpy = Load;                                       %Storing a copy of Global Load Vector
K_g_cpy(boundary, :) = [];                             %Eliminating the rows where boundary conditions are defined
K_g_cpy(:, boundary) = [];                             %Eiminating vector excluding boundary conditions
boundary_non = setxor(boundary,1:dof*n_n);             %Eliminating the corresponding columns
Load_cpy(boundary, :) = [];                            %Eliminating the corresponding forces in Global Load Vector
U_non = K_g_cpy\Load_cpy;                              %Calculating the displacementabbing the indices where boundary conditions are not defined
U(boundary_non(1:end)) = U_non(1:end);                  %Calculating the Global Displacement Vector
Reaction = (K_g*U) - Load;                              %Calculating Reactions at the nodes