%% Matlab Initializations

clear
clc
format shortEng

%% Given data in the problem

n_e = 2;                              %Number of elements
n_n = 4;                              %Number of nodes
dof = 2;                              %Dof of elements
E_e = repmat(200e9,n_e,1);            %Forming the matrix of Young's Modulus
nu_e = repmat(0.3,n_e,1);

%% Initializing and defining the required Global Vectors

Load = zeros(dof*n_n,1);              %Forming the Global Load Vector
U = zeros(dof*n_n,1);                 %Forming the Global displacement Vector
K_g = zeros(dof*n_n);                 %Initializing the Global Stiffness Matrix
N = [1/3 0 1/3 0 1/3 0];
B = [];
D = [];
Stress = [];

%% Defining the element connectivity and the global co-ordinates of each node

%Defining the element connectivity
nr = 1;
nc = 1;
for p=1:nr
    for q=1:nc
        cn=q+(1+nc)*(p-1);
        elems(2*nc*(p-1)+q,:)=[cn cn+1 cn+nc+2];
        elems((2*p-1)*nc+q,:)=[cn cn+nc+2 cn+nc+1];
        end
end

nodes = [0.04 0; 0.06 0; 0.04 0.01; 0.06 0.01];         %Defining the global co-ordinates of each node

%% Defining the Global Stiffness Matrix

for i = 1:n_e
    elnodes = elems(i,:);
    noderz = nodes(elnodes,:);
    r_cent = mean(noderz(:,1));
    diff_r = [noderz(3,1)-noderz(2,1) noderz(1,1)-noderz(3,1) noderz(2,1)-noderz(1,1)];       
    diff_z = [noderz(2,2)-noderz(3,2) noderz(3,2)-noderz(1,2) noderz(1,2)-noderz(2,2)];
    J = [diff_r(2) -diff_z(2); -diff_r(1) diff_z(1)];
    A_e = det(J)/2;                                                                                                                                        %Calculating the Area of Element
    B_e(1,1:2:5) = diff_z/det(J); B_e(2,2:2:6) = diff_r/det(J); B_e(3,2:2:6) = diff_z/det(J); B_e(3,1:2:5) = diff_r/det(J);
    B_e(4,:) = (1/r_cent)*N;                                                                                                                                      %Calculating the B Matrix
    B = [B;B_e];
    a = nu_e(i)/(1-nu_e(i));
    D_e = (E_e(i)*(1-nu_e(i))/((1+nu_e(i))*(1-(2*nu_e(i)))))*[1 a 0 a; a 1 0 a; 0 0 (1-(2*nu_e(i)))/(2*(1-nu_e(i))) 0; a a 0 1];                                                                         %Calculating the D Matrix
    D = [D;D_e];
    K_e = 2*pi*r_cent*A_e*B_e'*D_e*B_e;                                                                                                                        %Calculating the Local Stiffness Matrix
    eldofs = [(dof*(elnodes(1)-1))+1:(dof*(elnodes(1)-1)+2) (dof*(elnodes(2)-1))+1:(dof*(elnodes(2)-1)+2) (dof*(elnodes(3)-1))+1:(dof*(elnodes(3)-1)+2)];
    K_g(eldofs,eldofs) = K_g(eldofs,eldofs) + K_e;                                                                                                        %Assembling the Global Stiffness Matrix 
end

%% Defining the Boundary Conditions and Point Loads acting on various nodes

% Adding the loads due to Linearly distributed Load in load vector

ldl_side = [1 3];
ldl_load = [2e6 0];
for i = 1:size(ldl_side,1)
    elnodes = ldl_side(i,:);
    noderz = nodes(elnodes,:);
    trans = noderz(2,:) - noderz(1,:);
    l_e = norm(trans);  
    a = ((2*noderz(1,1))+noderz(2,1))/6 ; b = (noderz(1,1)+(2*noderz(2,1)))/6;
    T_e = (2*pi*l_e)*[a*ldl_load(1) a*ldl_load(2) b*ldl_load(1) b*ldl_load(2)];
    eldofs = [(dof*(elnodes(1)-1))+1:(dof*(elnodes(1)-1)+2) (dof*(elnodes(2)-1))+1:(dof*(elnodes(2)-1)+2)];
    Load(eldofs) = Load(eldofs)+ T_e';
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

%Load(4) = -30000;                        %Defining the Point Loads acting on various nodes
%}
boundary = [2 3 4 6 7 8];                 %Defining the Boundary Conditions on various nodes

%% Caculating the Global Displacement Vector and Global Load Vector using Elimination Approach

K_g_cpy = K_g;                                         %Storing a copy of Global Stiffness Matrix for pos-porcessing requirements
Load_cpy = Load;                                       %Storing a copy of Global Load Vector
K_g_cpy(boundary, :) = [];                             %Eliminating the rows where boundary conditions are defined
K_g_cpy(:, boundary) = [];                             %Eliminating the corresponding columns
Load_cpy(boundary, :) = [];                            %Eliminating the corresponding forces in Global Load Vector
U_non = K_g_cpy\Load_cpy;                              %Calculating the displacement vector excluding boundary conditions
boundary_non = setxor(boundary,1:dof*n_n);             %Grabbing the indices where boundary conditions are not defined
U(boundary_non(1:end)) = U_non(1:end);                 %Calculating the Global Displacement Vector
Reaction = (K_g*U) - Load;                             %Calculating Reactions at the nodes

disp('Displacement at C')
disp(U(1))
disp('Displacement at D')
disp(U(5))
%{
for i = 1:n_e
    elnodes = elems(i,:);
    D_e = [D((4*i)-3,:);D((4*i)-2,:);D((4*i)-1,:);D(4*i,:)];
    B_e = [B((4*i)-3,:);B((4*i)-2,:);B((4*i)-1,:);B(4*i,:)];
    eldofs = [(dof*(elnodes(1)-1))+1:(dof*(elnodes(1)-1)+2) (dof*(elnodes(2)-1))+1:(dof*(elnodes(2)-1)+2) (dof*(elnodes(3)-1))+1:(dof*(elnodes(3)-1)+2)];
    Stress_e = D_e*B_e*U(eldofs);
    Stress = [Stress;Stress_e];
end   

disp('Stresses in first element')
disp(Stress(1:4))
disp('Stresses in second element')
disp(Stress(5:8))
%}