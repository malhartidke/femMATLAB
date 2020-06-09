%% Matlab Initializations

clear;close;clc;
format shortEng

%% Given data in the problem

n_e = 2;                              %Number of elements
n_n = 3;                              %Number of nodes
dof = 2;                              %Dof of elements
E_e = repmat(200e9,n_e,1);            %Forming the matrix of Young's Modulus
rho_e = repmat(7840,n_e,1);           %Forming the matrix of Density
I_e = repmat(2000e-12,n_e,1);         %Forming the Moment of Interia Matrix
l_e = repmat(0.3,n_e,1);              %Forming the Length Matrix
A_e = repmat(240e-6,n_e,1);           %Forming the Area Matrix


%% Initializing and defining the required Global Vectors

K_g = zeros(dof*n_n);                 %Initializing the Global Stiffness Matrix
M_g = zeros(dof*n_n);                 %Initializing the Global Mass Matrix

%% Defining the Boundary Conditions and Point Loads acting on various nodes

boundary = [1 2];                 %Defining the Boundary Conditions on various nodes

%% Defining the Global Stiffness Matrix

for i = 1:n_e  
    eldofs = [(dof*(i-1))+1:(dof*(i-1)+2) (dof*((i+1)-1))+1:(dof*((i+1)-1)+2)];
    
    ak = 12; bk = 6*l_e(i); ck = 2*(l_e(i)^2); dk = 2*ck;
    K_e = ((E_e(i)*I_e(i))/(l_e(i)^3))*[ak bk -ak bk; bk dk -bk ck; -ak -bk ak -bk; bk ck -bk dk];     %Calculating the Local Stiffness Matrix
    K_g(eldofs,eldofs) = K_g(eldofs,eldofs) + K_e;                                                     %Assembling the Global Stiffness Matrix
    
    am = 22*l_e(i); bm = 13*l_e(i); cm = 3*(l_e(i)^2); dm = 4*(l_e(i)^2);
    M_e = ((rho_e(i)*A_e(i)*l_e(i))/420)*[156 am 54 -bm; am dm bm -cm; 54 bm 156 -am; -bm -cm -am dm]; %Calculating the Local Mass Matrix
    M_g(eldofs,eldofs) = M_g(eldofs,eldofs) + M_e;                                                     %Assembling the Global Mass Matrix
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
X = [0 0.3 0.6];

figure(1)
plot(X,mode_shape(1:dof:end,1))
hold on
plot(X,mode_shape(1:dof:end,2))
grid on
title('Plot of Modal Shapes')
legend('First Modal Shape','Second Modal Shape')
