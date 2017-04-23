% Reading the input data
%using importfile.m function to read the data fron .txt file
Input_cell=importfile('Input.txt');   
FIO = fopen('creatIO.m','w');
for i=1:numel(Input_cell)
    fwrite(FIO,Input_cell{i});
    fprintf(FIO,'\n');
end
creatIO;
close all; clc;
creatIO;
global Sw P dt

%% unit
if strcmp(unit,'field')
    unit_mod1 = 5.615; %STB to ft3
    unit_mod2 = 1.127e-3; %Field unit for darcy equation.
end

%% reservoir  
dx = Length/Nx;
dy = Width/Ny;     %ft
dz = Height/Nz;    

phi=zeros(Ny,Nx,Nz);
for i=1:Nz
    phi(:,:,i) = phi0;
end
phi=permute(phi,[3,1,2]);
k=ones(Nz,Ny,Nx)*k;

% fluid
Bo_0 = ones(Nz,Ny,Nx)*Bo_0;  
Bw_0 = ones(Nz,Ny,Nx)*Bw_0;

%% rock-fluid interpolation
sw_krw = interp1(SWT(:,1),SWT(:,2),'linear','pp');
sw_kro = interp1(SWT(:,1),SWT(:,3),'linear','pp');
sw_pc = interp1(SWT(:,1),SWT(:,4),'linear','pp');
Domain = SWT(1,1):0.005:1;
DpcDsw = diff(ppval(sw_pc,Domain))./diff(Domain);
sw_dpc_dsw = interp1((Domain(1:end-1)),DpcDsw,'linear','pp');

%% inital condition
X_init = zeros(2*Nz*Ny*Nx,1);   
Count=1;
for kk=1:Nz     
    for j=1:Ny
        for i=1:Nx
            X_init(Count)   = P_init(kk,j,i);   
            X_init(Count+1) = Sw_init(kk,j,i);
            Count = Count+2;
        end
    end
end
X0 = X_init; Sw = Sw_init; P = P_init;

%% Well equations
Qw_range = Qw;
Prod_p_range = pbh;
INJ_range = -Qw_range/(dx*dy*dz)* unit_mod1;
Re = sqrt(dy*dx/pi);

WCi = 2*pi*k(1,[1 end],[1 end])* dz/(log(Re/Rw))*unit_mod1*unit_mod2;  %STB/day*5.615=ft3/day   

P_bh_cell = cell(ceil(Max_days/dt)-1,1);
Q_water_cell = cell(ceil(Max_days/dt)-1,1); 
WCi_cell = cell(ceil(Max_days/dt)-1,1);

for i = 1:ceil(Max_days/dt)-1, P_bh_cell{i} = zeros(Nz,Ny,Nx);end
for i = 1:ceil(Max_days/dt)-1, Q_water_cell{i} = zeros(Nz,Ny,Nx);end

for i = 1:ceil(Max_days/dt)-1, WCi_cell{i} = zeros(Nz,Ny,Nx);end
for i = 1:ceil(Max_days/dt)-1,WCi_cell{i}(1,[1 end],[1 end])=WCi;end

% Injection Schedule
for i=1:ceil(Max_days/dt)-1,P_bh_cell{i}(1,1,1)=Prod_p_range;end      
for i=1:ceil(Max_days/dt)-1,P_bh_cell{i}(1,1,end)=Prod_p_range;end  
for i=1:ceil(Max_days/dt)-1,P_bh_cell{i}(1,end,1)=Prod_p_range;end   
for i=1:ceil(Max_days/dt)-1,P_bh_cell{i}(1,end,end)=Prod_p_range;end 

for i = 1:ceil(Max_days/dt)-1,Q_water_cell{i}(1,ceil(Ny/2),ceil(Nx/2)) = INJ_range;end   

varargin = {rho_oil,rho_water,SWT,sw_kro,sw_krw,sw_pc,sw_dpc_dsw,unit_mod1,unit_mod2 ...
          ,Length,dx,dy,Nx,Ny,k,phi,Height,Bo_0,Bw_0,Muo,Muw,Ct...
          ,Q_water_cell{1},Rw,Re,WCi_cell{1},P_bh_cell{1},Nz,dz};    



%% make and solve the equations
Itera = 1;
Sw_cell = cell(1,ceil(Max_days/dt)-1);
P_cell = cell(1,ceil(Max_days/dt)-1);
Landa_cell_O = cell(1,ceil(Max_days/dt)-1);
Landa_cell_W = cell(1,ceil(Max_days/dt)-1);
while Itera*dt<Max_days
 
    varargin{end-6} = Q_water_cell{Itera};
    varargin{end-2} = P_bh_cell{Itera};
    varargin{end-3}=  WCi_cell{Itera};
    Sw_cell{Itera} = Sw;
    P_cell{Itera} = P;
    
    [X,Landa_W,Landa_O] = EquSolve(@EquCalc,X0,1e-6,10,varargin{:});
    
    Count = 1;
    Landa_cell_W{Itera} = Landa_W;
    Landa_cell_O{Itera} = Landa_O;
    for kk=1:Nz                
        for j=1:Ny
            for i=1:Nx
               P(kk,j,i)= X(Count);   
               Sw(kk,j,i) = X(Count+1);  
               Count = Count+2;
            end
        end
    end
    Itera = Itera+1
    X0 = X;

end

 %% Result Visualization
x_well=1;
y_well=1;
z_well=1;

figure('name','Water Cut')
hold on;
for i=1:ceil(Max_days/dt)-1
    water_cut = Landa_cell_W{i}(z_well,y_well,x_well)/(Landa_cell_W{i}(z_well,y_well,x_well)+ ...
                Landa_cell_O{i}(z_well,y_well,x_well));   
    plot(i*dt,water_cut,'b.')

end
    xlabel('Days');ylabel('Water Cut');
    title('Water Cut');



E_Oil_rate=zeros(1,ceil(Max_days/dt)-1);
Cum_oil=zeros(ceil(Max_days/dt),1);

figure('name','Oil Production')
for i=1:ceil(Max_days/dt)-2
    E_Oil_rate(i)=WCi_cell{i}(z_well,y_well,x_well)*(Landa_cell_O{i}(z_well,y_well,x_well)*...
        (P_cell{i}(z_well,y_well,x_well)-P_bh_cell{i}(z_well,y_well,x_well)))/unit_mod1;  
    Cum_oil(i+1)=Cum_oil(i)+E_Oil_rate(i);
    
end
hold on;
plot([2:ceil(Max_days/dt)-1]*dt,E_Oil_rate(2:end),'b.')
xlabel('Days');ylabel('Oil Production Rate (STB/day)');
title('Oil Production Rate vs. Time')


figure('name','water Production')
E_Water_rate=zeros(1,ceil(Max_days/dt)-1);
for i=1:ceil(Max_days/dt)-2
    E_Water_rate(i)=WCi_cell{i}(z_well,y_well,x_well)*(Landa_cell_W{i}(z_well,y_well,x_well)*...
        (P_cell{i}(z_well,y_well,x_well)-P_bh_cell{i}(z_well,y_well,x_well)))/unit_mod1;  
end
hold on;
plot([1:ceil(Max_days/dt)-1]*dt,E_Water_rate,'c.',[1:ceil(Max_days/dt)-1]*dt,E_Water_rate,'-c')
xlabel('Days');ylabel('Water Production Rate (STB/day)');
title('Water Production Rate vs. Time')

% 
% 
figure('name','Cummulative Oil Production')
hold on
for i=1:ceil(Max_days/dt)-2  
    plot(i*dt,Cum_oil(i),'b.')
end
xlabel('Days');ylabel('Cummulative Oil production (STB)');
title('Cummulative Oil Production rate vs. Time')



                 % Water/Oil Production Curve
figure('name','Water-Oil Production Curve')
hold on
for i=1:ceil(Max_days/dt)-2
    water_oil_ratio=Landa_cell_W{i}(z_well,y_well,x_well)/Landa_cell_O{i}(z_well,y_well,x_well);
    plot(i*dt,water_oil_ratio,'b.')

end
xlabel('Days');ylabel('Water/Oil Production Ratio');
title('Water/Oil Production Ratio')

     

%             %Water Injection Downhole Pressure
% figure('name','Water Injection Downholde Pressure')
% B=4;C=4;
% for i=2:ceil(Max_days/dt)-2
%     hold on;
%     water_P=P_cell{i}(1,B,C);  
%     plot(i*dt,water_P-Q_water_cell{i}(1,B,C)*...
%         unit_mod1/(unit_mod2*WCi_cell{i}(1,1,1)*Landa_cell_W{i}(1,B,C)),'b.') 
%        
% end
% xlabel('Days');ylabel('Injection Pressure (Psi)');
% title('Water Injection Downhole Pressure (Psi)')
% legend('Well 1')

ExRecFact = Cum_oil(end-1)*unit_mod1/(phi(z_well,y_well,x_well)*Length^2*(1-SWT(1,1))*Height);

% 3D Visualization
% figure('name','Permeability Distribution')
% [x,y,z] = meshgrid(1:Nx,1:Ny,-1:-1:-Nz);
% xslice = Nx+1; yslice = Ny+1; zslice = -1:-1:-Nz;
% h1=slice(x,y,z,permute(k,[3,2,1]),xslice,yslice,zslice);
% colormap jet 
% set(h1,'facecolor','interp');
% grid on
% title('Permeability Distribution (md)')
% zlabel('Layer')
% colorbar
% 
% figure('name','Porosity Distribution')
% [x,y,z] = meshgrid(1:Nx,1:Ny,-1:-1:-Nz);
% xslice = Nx+1; yslice = Ny+1; zslice = -1:-1:-Nz;
% h1=slice(x,y,z,permute(phi,[3,2,1]),xslice,yslice,zslice);
% colormap jet 
% set(h1,'facecolor','interp');
% title('Porosity Distribution')
% zlabel('Layers')
% grid off
% colorbar

% figure('name','Pressure Propagation')
% for i=2:ceil(Max_days/dt)-1
% %     if i*dt==6 || i*dt==30 || i*dt==60 || i*dt==90 || i*dt==120 || i*dt==200 || i*dt==250
% %     figure('name',['Pressure Propagation after ' num2str(i*dt) ' days'])
%     [x,y,z] = meshgrid(1:Nx,1:Ny,-1:-1:-Nz);
%     xslice = Nx+1; yslice = Ny+1; zslice = -1:-1:-Nz;
%     h1=slice(x,y,z,permute(P_cell{i},[3,2,1]),xslice,yslice,zslice);
%     colormap jet 
%     set(h1,'facecolor','interp');
%     %caxis([min(min(min(P_cell{i}))) P_init(1)])
%     grid off
%     colorbar
%     pause(0.01)
% %     end
% end

% figure('name','Saturation Propagation')
% for i=2:ceil(Max_days/dt)-2
% %     if i*dt==150 || i*dt==180 || i*dt==210 || i*dt==240 || i*dt=270 || i*dt=294
% %     figure('name',['Pressure Propagation after ' num2str(i*dt) ' days'])
%     [x,y,z] = meshgrid(1:Nx,1:Ny,-1:-1:-Nz);
%     xslice = Nx+1; yslice = Ny+1; zslice = -1:-1:-Nz;
%     h1=slice(x,y,z,permute(Sw_cell{i},[3,2,1]),xslice,yslice,zslice);
%     set(h1,'facecolor','interp')
%     colormap hot
% %     caxis([0.2 max(max(max((Sw_cell{i}))))])
%     %set(h1,'facecolor','interp')  
%     colorbar
%     pause(0.01)
% %     end
%     i
% end

