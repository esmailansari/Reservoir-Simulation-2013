%---------------------
% Input/Output Units
%---------------------
% Field units conversion
unit = 'field'; 
%----------------------
% Reservoir description (one dimensional for buckley leverett).
%----------------------
Length=2000;  %ft 
Width=10;
Height=10;
Nx=35;
Ny=1;
Nz=1;
phi=0.22;
k=ones(Nx,1)*1000; %md
%------------------
% Fliud Description
%------------------
% Defining incompressible fluid.
Bw_0=1;
Bo_0=1;
Muo=1;   %cp
Muw=1;   
Pc = 0; % Capillary pressure
%------------------
% Rock-Fluid properties
%------------------
% Corey equation parameters
kro1=1;
krw1=1;
S_wc=0.2;
S_or=0.2;
a=2;b=2;
%-------------------
% Initial condition 
%-------------------
Sw=ones(Nx,1)*S_wc;  %Connate water
P=ones(Nx,1)*3000;  %psi
%----------------------
% Numerical data 
%----------------------
t_d = 0.01;   % for analytical solution. 
Max_days=300; %days
%----------------------
% Well definition
%----------------------
% choosing well location 
x_well_inj=1;
x_well_pro=Nx;
Qw=100; %STB/day
