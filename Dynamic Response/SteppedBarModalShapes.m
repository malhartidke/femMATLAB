%% Matlab Initializations

clear;close;clc;
format shortEng

%% Given data in the problem

A_e = [1 0.5];                  %Area of bar
n_e = 2;                        %Number of elements
n_n = 3;                        %Number of nodes
dof = 1;
E_e = repmat(3e7,n_e,1);        %Forming the matrix of Young's Modulus
rho_e = repmat(7.324e-4,n_e,1); %Forming the matrix of Density
L_e = [10 5];                   %Forming the length matrix

%% Initializing and defining the required Global Vectors

K_g = zeros(dof*n_n);           %Initializing the Global Stiffness Matrix
M_g = zeros(dof*n_n);           %Initializing the Global Mass Matrix

%% Defining the Boundary Conditions and Point Loads acting on various nodes

boundary = [1];               %Defining the Boundary Conditions on various nodes

%{
%% Converting the tapered bar into discrete elements with proper areas

element_divide = A_start:-((A_start-A_end)/n_e):A_end; %Dividing the bar into required number of elements
A_e = movmean(element_divide,2);                       %Calculating the area of each element
A_e = A_e(2:end);                                      %Forming the area vector
%}

%% Defining the Global Stiffness Matrix

for i = 1 : n_e
    eldofs = [(dof*(i-1))+1:(dof*(i-1)+2)];
    
    K_e = ((A_e(i)*E_e(i))/L_e(i))*[1 -1;-1 1];        %Calculating the Local Stiffness Matrix                            
    K_g(eldofs,eldofs) = K_g(eldofs,eldofs)+ K_e;      %Assembling the Global Stiffness Matrix
    
    M_e = ((rho_e(i)*A_e(i)*L_e(i))/6)*[2 1;1 2];      %Calculating the Local Mass Matrix                    
    M_g(eldofs,eldofs) = M_g(eldofs,eldofs)+ M_e;      %Assembling the Global Mass Matrix
end

%% Caculating the Global Displacement Vector and Global Load Vector using Elimination Approach

K_g_cpy = K_g;                                         %Storing a copy of Global Stiffness Matrix for pos-porcessing requirements
K_g_cpy(boundary, :) = [];                             %Eliminating the rows where boundary conditions are defined
K_g_cpy(:, boundary) = [];                             %Eliminating the corresponding columns

M_g_cpy = M_g;                                         %Storing a copy of Global Mass Matrix for pos-porcessing requirements
M_g_cpy(boundary, :) = [];                             %Eliminating the rows where boundary conditions are defined
M_g_cpy(:, boundary) = [];                             %Eliminating the corresponding columns

[e_vec,e_val] = eig(K_g_cpy,M_g_cpy);
[omg_2,ind] = sort(diag(e_val));
u = e_vec(:,ind);

mode_shape = zeros(dof*n_n,size(u,2));                 %Forming the Global Mode shape Vector
omg = (sqrt(omg_2)/(2*pi))

boundary_non = setxor(boundary,1:(dof*n_n));           %Grabbing the indices where boundary conditions are not defined
mode_shape(boundary_non(1:end),:) = u(1:end,:);        %Calculating the Global Displacement Vector
X = [0 10 15];

figure(1)
plot(X,mode_shape(1:dof:end,1))
hold on
plot(X,mode_shape(1:dof:end,2))
grid on
title('Plot of Modal Shapes')
legend('First Modal Shape','Second Modal Shape')