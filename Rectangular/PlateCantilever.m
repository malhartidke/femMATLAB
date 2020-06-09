%% Matlab Initializations

clear;clc;
format shortEng;

%% Given data in the problem

n_e = 4;                              %Number of elements
n_n = 9;                              %Number of nodes
dof = 2;                              %Dof of elements
E_e = repmat(200e9,n_e,1);            %Forming the matrix of Young's Modulus
t_e = repmat(0.02,n_e,1);
nu_e = repmat(0.3,n_e,1);

%For 2 Point Quadrature
ip = [-(1/sqrt(3)) (1/sqrt(3))];
w = [1 1];

%For 3 Point Quadrature
%ip = [0.7745966692 -0.7745966692 0];
%w = [0.5555555556 0.5555555556 0.8888888889];

%For 4 Point Quadrature
%ip = [0.8611363116 -0.8611363116 0.3399810436 -0.3399810436];
%w = [0.3478548451 0.3478548451 0.6521451549 0.6521451549];

%% Initializing and defining the required Global Vectors

Load = zeros(dof*n_n,1);              %Forming the Global Load Vector
U = zeros(dof*n_n,1);                 %Forming the Global displacement Vector
K_g = zeros(dof*n_n);                 %Initializing the Global Stiffness Matrix

%% Defining the element connectivity and the global co-ordinates of each node

elems = [1 2 5 4; 2 3 6 5; 4 5 8 7; 5 6 9 8];                  %Defining the element connectivity
nodes = [0 0; 1 0; 2 0; 0 1; 1 1; 2 1; 0 2; 1 2; 2 2];         %Defining the global co-ordinates of each node

%% Defining the Boundary Conditions and Point Loads acting on various nodes

Body_F = zeros(dof,n_e);                                     % Initializing the Body force vector
%{
% Adding the loads due to Linearly distributed Load in load vector
ldl_e = [2 4];                                               % Elements on which load is acting
ldl_side = [3 6; 6 9];                                       % Nodes on which load is acting
ldl_load = [2e6; 3e6]; 
for i = 1:size(ldl_side,1)
    elnodes = ldl_side(i,:);
    nodexy = nodes(elnodes,:);
    trans = nodexy(2,:) - nodexy(1,:);
    l_e = norm(trans);             
    trans = trans/l_e;
    Tx = ldl_load(i)*trans(2); Ty = ldl_load(i)*(-trans(1));
    T_e = ((t_e(ldl_e(i))*l_e)/2)*[Tx Ty Tx Ty];
    eldofs = [(dof*(elnodes(1)-1))+1:(dof*(elnodes(1)-1))+2 (dof*(elnodes(2)-1))+1:(dof*(elnodes(2)-1))+2];
    Load(eldofs) = Load(eldofs)+ T_e';
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

Load(18) = -30000;                          %Defining the Point Loads acting on various nodes
boundary = [1 2 7 8 13 14];                 %Defining the Boundary Conditions on various nodes

%% Forming the Global Stiffness Matrix and the Body Force Vector

for n = 1:n_e
    elnodes = elems(n,:);
    nodexy = nodes(elnodes,:);
    phi = zeros(dof*4,dof*4);
    phi_L = zeros(dof*4,2);
    for i = 1:size(ip,2)     
        for j = 1:size(ip,2)
            zeta = ip(i);
            eta = ip(j);
            J1 = [-(1-eta) 1-eta 1+eta -(1+eta)];
            J2 = [-(1-zeta) -(1+zeta) 1+zeta 1-zeta];
            J = (1/4)*[J1*nodexy(:,1) J1*nodexy(:,2); J2*nodexy(:,1) J2*nodexy(:,2)];
            A = (1/det(J))*[J(2,2) -J(1,2) 0 0; 0 0 -J(2,1) J(1,1); -J(2,1) J(1,1) J(2,2) -J(1,2)];
            G1 = reshape([J1' zeros(size(J1,2),1)]',1,[]);
            G2 = reshape([J2' zeros(size(J2,2),1)]',1,[]);
            G = (1/4)*[G1; G2; 0 G1(1:(end-1)); 0 G2(1:(end-1))];
            B = A*G;
            D = (E_e(n)/(1-(nu_e(n)^2)))*[1 nu_e(n) 0; nu_e(n) 1 0; 0 0 (1-nu_e(n))/2];
            phi = (w(i)*w(j)*B'*D*B*det(J))+ phi;
            N1 = (1/4)*(1-zeta)*(1-eta); N2 = (1/4)*(1+zeta)*(1-eta); N3 = (1/4)*(1+zeta)*(1+eta); N4 = (1/4)*(1-zeta)*(1+eta);
            N = [N1 0 N2 0 N3 0 N4 0;0 N1 0 N2 0 N3 0 N4];
            phi_L = (w(i)*w(j)*N'*det(J)) + phi_L;
        end    
    end
    K_e = t_e(n)*phi;
    eldofs = [(dof*(elnodes(1)-1))+1:(dof*(elnodes(1)-1))+2 (dof*(elnodes(2)-1))+1:(dof*(elnodes(2)-1))+2 (dof*(elnodes(3)-1))+1:(dof*(elnodes(3)-1))+2 (dof*(elnodes(4)-1))+1:(dof*(elnodes(4)-1))+2];
    K_g(eldofs,eldofs) = K_g(eldofs,eldofs) + K_e;
    f_e = t_e(n)*phi_L*Body_F(:,n_e);
    Load(eldofs) = Load(eldofs) + f_e;
end

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

disp('Displacement at point A is:')
disp(U(17))
disp(U(18))

disp('Displacement at point B is:')
disp(U(5))
disp(U(6))

disp('Reaction at point C is:')
disp(Reaction(1))
disp(Reaction(2))

disp('Reaction at point D is:')
disp(Reaction(13))
disp(Reaction(14))

%{
for i = 1:n_e
    elnodes = elems(i,:);
    D_e = [D((3*i)-2,:);D((3*i)-1,:);D((3*i),:)];
    B_e = [B((3*i)-2,:);B((3*i)-1,:);B((3*i),:)];
    eldofs = [(dof*(elnodes(1)-1))+1:(dof*(elnodes(1)-1)+2) (dof*(elnodes(2)-1))+1:(dof*(elnodes(2)-1)+2) (dof*(elnodes(3)-1))+1:(dof*(elnodes(3)-1)+2)];
    Stress_e = D_e*B_e*U(eldofs);
    Stress = [Stress;Stress_e];
end   
%}