%% Matlab Initializations

clear
clc
format shortEng

%% Given data in the problem

n_e = 2;                              %Number of elements
n_n = 4;                              %Number of nodes
dof = 2;                              %Dof of elements
E_e = repmat(30e6,n_e,1);             %Forming the matrix of Young's Modulus
t_e = repmat(0.5,n_e,1);
nu_e = repmat(0.25,n_e,1);

%% Initializing and defining the required Global Vectors

Load = zeros(dof*n_n,1);              %Forming the Global Load Vector
U = zeros(dof*n_n,1);                 %Forming the Global displacement Vector
K_g = zeros(dof*n_n);                 %Initializing the Global Stiffness Matrix
B = [];
D = [];
Stress = [];

%% Defining the element connectivity and the global co-ordinates of each node

elems = [1 2 4; 3 4 2];               %Defining the element connectivity
nodes = [3 0; 3 2; 0 2; 0 0];         %Defining the global co-ordinates of each node

%% Defining the Global Stiffness Matrix

for i = 1:n_e
    elnodes = elems(i,:);
    nodexy = nodes(elnodes,:);
    diff_x = [nodexy(3,1)-nodexy(2,1) nodexy(1,1)-nodexy(3,1) nodexy(2,1)-nodexy(1,1)];       
    diff_y = [nodexy(2,2)-nodexy(3,2) nodexy(3,2)-nodexy(1,2) nodexy(1,2)-nodexy(2,2)];
    Det_J = (diff_x(2)*diff_y(1)) - (diff_x(1)*diff_y(2));
    A_e = Det_J/2;                                                                                                                                        %Calculating the Area of Element
    B_e(1,1:2:5) = diff_y; B_e(2,2:2:6) = diff_x; B_e(3,2:2:6) = diff_y; B_e(3,1:2:5) = diff_x;                       
    B_e = B_e/Det_J;                                                                                                                                      %Calculating the B Matrix
    B = [B;B_e];
    D_e = E_e(i)*(1-(nu_e(i)^2))*[1 nu_e(i) 0; nu_e(i) 1 0; 0 0 (1-nu_e(i))/2];                                                                           %Calculating the D Matrix
    D = [D;D_e];
    K_e = t_e(i)*A_e*B_e'*D_e*B_e;                                                                                                                        %Calculating the Local Stiffness Matrix
    eldofs = [(dof*(elnodes(1)-1))+1:(dof*(elnodes(1)-1)+2) (dof*(elnodes(2)-1))+1:(dof*(elnodes(2)-1)+2) (dof*(elnodes(3)-1))+1:(dof*(elnodes(3)-1)+2)];
    K_g(eldofs,eldofs) = K_g(eldofs,eldofs) + K_e;                                                                                                      %Assembling the Global Stiffness Matrix 
end

%{
%% Defining the Boundary Conditions and Point Loads acting on various nodes

% Adding the loads due to Linearly distributed Load in load vector

ldl = [2];
ldl_side = [2 3];
ldl_load = [2e6 3e6];
for i = 1:size(ldl_side,1)
    elnodes = ldl_side(i,:);
    nodexy = nodes(elnodes,:);
    trans = nodexy(2,:) - nodexy(1,:);
    l_e = nom(trans);             
    trans = trans/l_e;
    Tx_1 = ldl_load(i,1)*trans(2); Ty_1 = ldl_load(i,1)*trans(1); Tx_2 = ldl_load(i,2)*trans(2); Ty_2 = ldl_load(i,2)*trans(1);
    T_e = ((t_e(ldl(i))*l_e)/6)*[(2*Tx_1)+Tx_2 (2*Ty_1)+Ty_2 Tx_1+(2*Tx_2) Ty_1+(2*Ty_2)];
    eldofs = [(dof*(elnodes(1)-1))+1:(dof*(elnodes(1)-1)+3) (dof*(elnodes(2)-1))+1:(dof*(elnodes(2)-1)+3)];
    Load(eldofs) = Load(eldofs)+ T_e;
end

% Adding the Point loads on the centre of element in load vector
point_centre = [1];
point_centre_load = [12];
for i = 1:length(point_centre)
    reaction = point_centre_load(i)/2;
    couple = (point_centre_load(i)*l_e(point_centre(i)))/8;
    eldofs = [(dof*(point_centre(i)-1))+1:(dof*(point_centre(i)-1)+3) (dof*((point_centre(i)+1)-1))+1:(dof*((point_centre(i)+1)-1)+3)];
    Load(eldofs) = Load(eldofs)+ [0 -reaction -couple 0 -reaction couple]';
end

%}

Load(4) = -1000;                        %Defining the Point Loads acting on various nodes
boundary = [2 5 6 7 8];                 %Defining the Boundary Conditions on various nodes

%% Caculating the Global Displacement Vector and Global Load Vector using Elimination Approach

K_g_cpy = K_g;                                         %Storing a copy of Global Stiffness Matrix for pos-porcessing requirements
Load_cpy = Load;                                       %Storing a copy of Global Load Vector
K_g_cpy(boundary, :) = [];                             %Eliminating the rows where boundary conditions are defined
K_g_cpy(:, boundary) = [];                             %Eliminating the corresponding columns
Load_cpy(boundary, :) = [];                            %Eliminating the corresponding forces in Global Load Vector
U_non = K_g_cpy\Load_cpy;                              %Calculating the displacement vector excluding boundary conditions
boundary_non = setxor(boundary,1:dof*n_n);             %Grabbing the indices where boundary conditions are not defined
U(boundary_non(1:end)) = U_non(1:end);                 %Calculating the Global Displacement Vector
Reaction = (K_g*U) - Load                             %Calculating Reactions at the nodes

for i = 1:n_e
    elnodes = elems(i,:);
    D_e = [D((3*i)-2,:);D((3*i)-1,:);D((3*i),:)];
    B_e = [B((3*i)-2,:);B((3*i)-1,:);B((3*i),:)];
    eldofs = [(dof*(elnodes(1)-1))+1:(dof*(elnodes(1)-1)+2) (dof*(elnodes(2)-1))+1:(dof*(elnodes(2)-1)+2) (dof*(elnodes(3)-1))+1:(dof*(elnodes(3)-1)+2)];
    Stress_e = D_e*B_e*U(eldofs);
    Stress = [Stress;Stress_e];
end    

disp('Displacment at node 1:')
disp(U(1))
disp('Displacement at node 2:')
disp('u = ')
disp(U(3))
disp('v = ')
disp(U(4))
%disp('Stresses in first element')
%disp(Stress(1:3))
%disp('Stresses in second element')
%disp(Stress(4:6))
