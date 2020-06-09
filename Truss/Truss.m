%% Matlab Initializations

clear
clc

%% Given data in the problem

n_e = 7;                              %Number of elements
n_n = 5;                              %Number of nodes
dof = 2;                              %Dof of elements
E_e = repmat(2e11,n_e,1);             %Forming the matrix of Young's Modulus
A_e = repmat(1e-5,n_e,1);             %Forming the Area matrix
l_e = zeros(n_e,1);
force = zeros(n_e,1);

%% Initializing and defining the required Global Vectors

F = zeros(dof*n_n,1);                 %Forming the Global Load Vector
U = zeros(dof*n_n,1);                 %Forming the Global displacement Vector
K_g = zeros(dof*n_n);                 %Initializing the Global Stiffness Matrix

%% Defining the Boundary Conditions and Point Loads acting on various nodes

F(3) = 1e4; F(8) = -1e4;               %Defining the Point Loads acting on various nodes
boundary = [1 2 5 9 10];               %Defining the Boundary Conditions on various nodes

%% Defining the element connectivity and the global co-ordinates of each node

elems = [1 2; 1 3; 3 2; 2 4; 2 5; 3 4; 4 5];            %Defining the element connectivity
nodes = [0 4000; 3000 4000; 0 0; 3000 0; 6000 0];       %Defining the global co-ordinates of each node


%% Defining the Global Stiffness Matrix

for i = 1:n_e
    elnodes = elems(i,:);
    nodexy = nodes(elnodes,:);
    trans = nodexy(2,:)-nodexy(1,:);
    l_e(i) = norm(trans);                                                                                    %Calculating the length of element
    trans = trans/l_e(i);
    transformation = [trans zeros(1,2);zeros(1,2) trans];                                                    %Creating the transformational matrix
    K_e_loc = ((A_e(i)*E_e(i))/l_e(i))*[1 -1;-1 1];                                                          %Calculating the Local Stiffness Matrix in local co-ordinates
    K_e_g = transformation'*K_e_loc*transformation;                                                          %Calculating the Local Stiffness Matrix in global co-ordinates
    eldofs = [(dof*(elnodes(1)-1))+1:(dof*(elnodes(1)-1)+2) (dof*(elnodes(2)-1))+1:(dof*(elnodes(2)-1)+2)];
    K_g(eldofs,eldofs) = K_g(eldofs,eldofs) + K_e_g;                                                         %Assembling the Global Stiffness Matrix 
end

%% Caculating the Global Displacement Vector and Global Load Vector using Elimination Approach

K_g_cpy = K_g;                                         %Storing a copy of Global Stiffness Matrix for pos-porcessing requirements
F_cpy = F;                                             %Storing a copy of Global Load Vector
K_g_cpy(boundary, :) = [];                             %Eliminating the rows where boundary conditions are defined
K_g_cpy(:, boundary) = [];                             %Eliminating the corresponding columns
F_cpy(boundary, :) = [];                               %Eliminating the corresponding forces in Global Load Vector
U_non = K_g_cpy\F_cpy;                                 %Calculating the displacement vector excluding boundary conditions
boundary_non = setxor(boundary,1:dof*n_n);             %Grabbing the indices where boundary conditions are not defined
U(boundary_non(1:end)) = U_non(1:end)                  %Calculating the Global Displacement Vector
F = K_g * U

%% Calculating the forces in elements

for i = 1:n_e
    elnodes = elems(i,:);
    nodexy = nodes(elnodes, :);
    trans = nodexy(2,:)-nodexy(1,:);
    trans = trans/l_e(i);
    transformation = [-trans trans];
    force(i) = (E_e(i)*A_e(i)/l_e(i))*transformation*U([elnodes(1)+1 elnodes(1)+2 elnodes(2)+1 elnodes(2)+2]);
end
disp('The force in elements is:')
disp(force)