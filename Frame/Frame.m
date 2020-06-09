%% Matlab Initializations

clear
clc
format shortEng

%% Given data in the problem

n_e = 3;                              %Number of elements
n_n = 4;                              %Number of nodes
dof = 3;                              %Dof of elements
E_e = repmat(30e6,n_e,1);             %Forming the matrix of Young's Modulus
I_e = repmat(65,n_e,1);               %Forming the Moment of Interia Matrix
A_e = repmat(6.8,n_e,1);              %Forming the Area matrix
l_e = zeros(n_e,1);

%% Initializing and defining the required Global Vectors

Load = zeros(dof*n_n,1);              %Forming the Global Load Vector
U = zeros(dof*n_n,1);                 %Forming the Global displacement Vector
K_g = zeros(dof*n_n);                 %Initializing the Global Stiffness Matrix

%% Defining the element connectivity and the global co-ordinates of each node

elems = [1 3; 2 4; 3 4];                  %Defining the element connectivity
nodes = [0 0; 144 0; 0 96; 144 96];       %Defining the global co-ordinates of each node

%% Defining the Global Stiffness Matrix

for i = 1:n_e
    elnodes = elems(i,:);
    nodexy = nodes(elnodes,:);
    trans = nodexy(2,:) - nodexy(1,:);
    l_e(i) = norm(trans);                                                                                                                                            %Calculating the length of element
    trans = trans/l_e(i);
    transformation = [trans zeros(1,4); -trans(2) trans(1) zeros(1,4); zeros(1,2) 1 zeros(1,3); zeros(1,3) trans 0; zeros(1,3) -trans(2) trans(1) 0; zeros(1,5) 1];  %Creating the transformational matrix
    a = ((E_e(i)*A_e(i))/l_e(i)); b = 12*((E_e(i)*I_e(i))/(l_e(i)^3)); c = 6*((E_e(i)*I_e(i))/(l_e(i)^2)); d = 4*((E_e(i)*I_e(i))/l_e(i));
    K_e_loc = [a 0 0 -a 0 0; 0 b c 0 -b c; 0 c d 0 -c d/2; -a 0 0 a 0 0; 0 -b -c 0 b -c; 0 c d/2 0 -c d];                                                             %Calculating the Local Stiffness Matrix in local co-ordinates
    K_e_g = transformation'*K_e_loc*transformation;                                                                                                                  %Calculating the Local Stiffness Matrix in global co-ordinates
    eldofs = [(dof*(elnodes(1)-1))+1:(dof*(elnodes(1)-1)+3) (dof*(elnodes(2)-1))+1:(dof*(elnodes(2)-1)+3)];
    K_g(eldofs,eldofs) = K_g(eldofs,eldofs) + K_e_g;                                                                                                                 %Assembling the Global Stiffness Matrix 
end

%% Defining the Boundary Conditions and Point Loads acting on various nodes

% Adding the loads due to UDL in load vector
udl = [3];
udl_load = [500/12];
for i = 1:length(udl)
    elnodes = elems(udl(i),:);
    reaction = (udl_load(i)*l_e(udl(i)))/2;
    couple = (udl_load(i)*(l_e(udl(i))^2))/12;
    eldofs = [(dof*(elnodes(1)-1))+1:(dof*(elnodes(1)-1)+3) (dof*(elnodes(2)-1))+1:(dof*(elnodes(2)-1)+3)];
    Load(eldofs) = Load(eldofs)+ [ 0 -reaction -couple 0 -reaction couple]';
end

%{
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

Load(7) = 3e3;                        %Defining the Point Loads acting on various nodes
boundary = [1 2 3 4 5 6];             %Defining the Boundary Conditions on various nodes

%% Caculating the Global Displacement Vector and Global Load Vector using Elimination Approach

K_g_cpy = K_g;                                         %Storing a copy of Global Stiffness Matrix for pos-porcessing requirements
Load_cpy = Load;                                       %Storing a copy of Global Load Vector
K_g_cpy(boundary, :) = [];                             %Eliminating the rows where boundary conditions are defined
K_g_cpy(:, boundary) = [];                             %Eliminating the corresponding columns
Load_cpy(boundary, :) = [];                            %Eliminating the corresponding forces in Global Load Vector
U_non = K_g_cpy\Load_cpy;                              %Calculating the displacement vector excluding boundary conditions
boundary_non = setxor(boundary,1:dof*n_n);             %Grabbing the indices where boundary conditions are not defined
U(boundary_non(1:end)) = U_non(1:end)                  %Calculating the Global Displacement Vector
Reaction = (K_g*U) - Load                              %Calculating Reactions at the nodes