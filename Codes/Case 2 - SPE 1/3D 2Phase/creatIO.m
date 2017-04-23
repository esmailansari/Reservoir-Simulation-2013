%---------------------
% Input/Output Units
%---------------------
% Field units conversion
unit = 'field'; 
%----------------------
% Reservoir description (one dimensional for buckley leverett).
%----------------------
Length=200;  %ft 
Width=200;
Height=30;
Nx=7;
Ny=7; 
Nz=3;
phi0=0.22;
k=10;                %md (homogeneous reservoir)
Ct=45e-6;            %psi^-1
%------------------
% Fliud Description
%------------------
% Defining incompressible fluid.
Bw_0=1;
Bo_0=1;
Muo=3;   %cp
Muw=1;   
rho_oil = 45;
rho_water = 60;
%------------------
% Rock-Fluid properties
%------------------
SWT=[0.25   0        1         70
0.3    0.0003   0.813    45
0.35   0.0024   0.651    30
0.4    0.008    0.512    16
0.45   0.019    0.394    12
0.5    0.037    0.269     7
0.6    0.1017   0.152     4
0.7    0.216    0.064     2
1      1        0         0];
%-------------------
% Initial condition 
%-------------------
Sw_init = ones(Nz,Ny,Nx)*0.25; 
P_init = ones(Nz,Ny,Nx)*3000; %psi
%----------------------
% Numerical data 
%----------------------
dt = 1;   %days
Max_days = 400; %days
%----------------------
% Well definition
%----------------------
% Injection well (one well)
x_well_inj=10;
y_well_inj=10;
Qw = 200;  %  STB/day
% Production wells (4 wells) 
pbh = 2500;  %psi  
Rw=0.33;              %ft
