%------------------
% Pressure Transient Data
%-------------------
Max_days=1.5; %days. total Pressure Transient time
shut_in=1;    %days. Shut in time for Build up test. 
%---------------------
% Input/Output Units
%---------------------
% Field units conversion
unit = 'field'; 
%----------------------
% Reservoir description (one dimensional for buckley leverett).
%----------------------
Length=200;  % ft
Width=200;
Height=30;
Nx=20;
Ny=20;      %Nz=Nx, square griding. Instead Change width.  
Nz=1;       %2D. Instead Change Height. 
phi=0.22;
k=50;   %md  %All grids have the permeability of 50.
%------------------
% Fliud Description
%------------------
% Defining incompressible fluid.
Bo_0=1;
Muo=1;   %cp  
Ct=45e-6;        %psi-1
%------------------
% Rock-Fluid properties
%------------------
% Single phase flow. 
%-------------------
% Initial condition 
%-------------------
P=ones(Ny,Nx)*3000;  %psi
%----------------------
% Numerical data 
%----------------------
dt=1/100;   %days
%----------------------
% Well definition
%----------------------
% choosing well location 
x_well=10;
y_well=10;
Qo=2000; %STB/day
