%% Matlab Initializations

clear
clc

%% Given data in the problem

t = 25;                       %Thickness of the bar
l = 500;                      %Total length of bar
A_start = 150*t;              %Starting Area of bar
A_end = 70*t;                 %Ending Area of bar
n_e = 5;                      %Number of elements
n_n = n_e + 1;
E_e = repmat(200000,n_e,1);   %Forming the matrix of Young's Modulus
L_e = repmat(l/n_e,n_e,1);    %Forming the length matrix

%% Initializing and defining the required Global Vectors

F = zeros(n_n,1);           %Forming the Global Load Vector
U = zeros(n_n,1);           %Forming the Global displacement Vector
K_g = zeros(n_n);           %Initializing the Global Stiffness Matrix
add_e = 0;                

%% Defining the Boundary Conditions and Point Loads acting on various nodes

F(6) = 5000;                  %Defining the Point Loads acting on various nodes
boundary = [1];               %Defining the Boundary Conditions on various nodes

%% Converting the tapered bar into discrete elements with proper areas

element_divide = A_start:-((A_start-A_end)/n_e):A_end; %Dividing the bar into required number of elements
A_e = movmean(element_divide,2);                       %Calculating the area of each element
A_e = A_e(2:end);                                      %Forming the area vector

%% Defining the Global Stiffness Matrix

for i = 1 : n_e
    K_e = ((A_e(i)*E_e(i))/L_e(i))*[1 -1;-1 1];        %Calculating the Local Stiffness Matrix
    K_g(i:i+1,i:i+1) = K_e;                            %Assembling the Global Stiffness Matrix 
    K_g(i,i) = K_g(i,i)+ add_e;
    add_e = K_g(i+1,i+1);
end

%% Caculating the Global Displacement Vector and Global Load Vector using Elimination Approach

K_g_cpy = K_g;                                         %Storing a copy of Global Stiffness Matrix for pos-porcessing requirements
F_cpy = F;                                             %Storing a copy of Global Load Vector
K_g_cpy(boundary, :) = [];                             %Eliminating the rows where boundary conditions are defined
K_g_cpy(:, boundary) = [];                             %Eliminating the corresponding columns
F_cpy(boundary, :) = [];                               %Eliminating the corresponding forces in Global Load Vector
U_non = K_g_cpy\F_cpy;                                 %Calculating the displacement vector excluding boundary conditions
boundary_non = setxor(boundary,1:n_n);                 %Grabbing the indices where boundary conditions are not defined
U(boundary_non(1:end)) = U_non(1:end)                  %Calculating the Global Displacement Vector
F = K_g * U                                            %Calculating the Global Load Vector