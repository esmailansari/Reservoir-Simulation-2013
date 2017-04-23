%---------------------
% Input/Output Units
%---------------------
% Field units conversion
unit = 'field'; 
%----------------------
% Reservoir description 
%----------------------
Length=200;  %ft 
Width=200;
Height=30;
Nx=20;
Ny=20; % Change Width
phi=0.22;
k=ones(Ny,Nx)*100;   %md
Ct=45e-6;            %psi^-1
%------------------
% Fliud Description
%------------------
% Defining incompressible fluid.
Bw_0=1;
Bo_0=1;
Muo=1;   %cp
Muw=1;   
%------------------
% Rock-Fluid properties
%------------------
SWT=[0.25   0        1       70
0.3    0.0003   0.813   45
0.35   0.0024   0.651   30
0.4    0.008    0.512   16
0.45   0.019    0.394   12
0.5    0.037    0.269    7
0.6    0.1017   0.152    4
0.7    0.216    0.064    2
1      1        0        0];
%-------------------
% Initial condition 
%-------------------
Sw=ones(Ny,Nx)*0.25;  %Connate water
P=ones(Ny,Nx)*3000;   %psi
%----------------------
% Numerical data 
%----------------------
dt=1/4;         %days
Max_days=300;   %days
%----------------------
% Well definition
%----------------------
% Injection well (one well)
x_well_inj=10;
y_well_inj=10;
Qw=200;        %STB/day
% Production wells (well 1)
x_well_prod_1=1;
y_well_prod_1=1; 
P_bh=3000;   
Rw=0.33;     %ft
% Production wells (well 1)
x_well_prod_2=20;
y_well_prod_2=1; 
P_bh=3000;   
Rw=0.33;     %ft
% Production wells (well 3)
x_well_prod_3=1;
y_well_prod_3=20; 
P_bh=3000;   
Rw=0.33;     %ft
% Production wells (well 4)
x_well_prod_4=20;
y_well_prod_4=20; 
P_bh=3000;   
Rw=0.33;     %ft
