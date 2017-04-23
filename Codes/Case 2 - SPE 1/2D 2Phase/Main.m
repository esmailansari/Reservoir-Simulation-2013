% Reading the input data
%using importfile.m function to read the data fron .txt file
Input_cell=importfile('Input.txt');   
FIO = fopen('creatIO.m','w');
for i=1:numel(Input_cell)
    fwrite(FIO,Input_cell{i});
    fprintf(FIO,'\n');
end
fclose(FIO);
creatIO;
close all; clc;

%% making variables
TOY=zeros(Ny,Nx);TWY=zeros(Ny,Nx);
TOX=zeros(Ny,Nx);TWX=zeros(Ny,Nx);
A=zeros(Nx*Ny);
Q_water=zeros(Ny,Nx);Q_oil=zeros(Ny,Nx);
Sw_t=zeros(Ny,Nx);
b=zeros(Ny*Nx,1);
Pcow=zeros(Ny,Nx);
kro=zeros(Ny,Nx);
krw=zeros(Ny,Nx);
Result=cell(4,floor(Max_days/dt));

%% Defining parameters
dx=Length/Nx;
dy=Width/Ny;
h=Height;  

% linear interpolation for unknown points
S_wc=SWT(1,1);
sw_krw = interp1(SWT(:,1),SWT(:,2),'linear','pp');
sw_kro = interp1(SWT(:,1),SWT(:,3),'linear','pp');
sw_pc = interp1(SWT(:,1),SWT(:,4),'linear','pp');
Domain = S_wc:0.005:1;
%calculating Dpc/Dsw.
DpcDsw = diff(ppval(sw_pc,Domain))./diff(Domain);
sw_dpc_dsw = interp1((Domain(1:end-1)),DpcDsw,'linear','pp');

%unit
if strcmp(unit,'field')
    unit_mod1 = 5.615; %STB to ft3
    unit_mod2 = 1.127e-3; %Field unit for darcy equation.
end

Bo_0=ones(Ny,Nx)*Bo_0;
Bw_0=ones(Ny,Nx)*Bw_0;

%% Well equations
%injector
injection=-Qw/(dx*dy*h)* unit_mod1;  %STB/day/ft3*5.615=1/day  %injection must be minus.

Re=sqrt(dy*dx/pi);
WCi = zeros(Ny,Nx);
WCi(y_well_inj,x_well_inj)=2*pi*k(1,1)*h/(log(Re/Rw))*unit_mod1*unit_mod2;  %STB/day*5.615=ft3/day
WCi(y_well_prod_1,x_well_prod_1)=2*pi*k(1,1)*h/(log(Re/Rw))*unit_mod1*unit_mod2;  %STB/day*5.615=ft3/day   
WCi(y_well_prod_2,x_well_prod_2)=2*pi*k(1,1)*h/(log(Re/Rw))*unit_mod1*unit_mod2;  %STB/day*5.615=ft3/day   
WCi(y_well_prod_3,x_well_prod_3)=2*pi*k(1,1)*h/(log(Re/Rw))*unit_mod1*unit_mod2;  %STB/day*5.615=ft3/day   
WCi(y_well_prod_4,x_well_prod_4)=2*pi*k(1,1)*h/(log(Re/Rw))*unit_mod1*unit_mod2;  %STB/day*5.615=ft3/day   

%% Defining landa and flow rate.
Iteration=1;
P_t=P;
while Iteration<floor(Max_days/dt)
    
    %% Calculating Mobility of each grid
      Bw=Bw_0;%.*(1+C_water.*(P_t-P));  
      Bo=Bo_0;%.*(1+C_oil.*(P_t-P));  
      kro=ppval(sw_kro,Sw);
      krw=ppval(sw_krw,Sw);
      kro(kro>SWT(end,1))=SWT(end,1);kro(kro<0)=0;   
      krw(krw>SWT(end,1))=SWT(end,1);krw(krw<0)=0;  
      Landa_O=kro./(Muo.*Bo);
      Landa_W=krw./(Muw.*Bw);
      
    %% Defining parameters
    C_swo=-phi./(Bo.*dt).*ones(Ny,Nx);
    C_pow=phi.*Sw./dt.*(Ct./Bo);
    dpc_dsw=ppval(sw_dpc_dsw,Sw);
    C_sww=(phi./(Bw.*dt).*ones(Ny,Nx)-dpc_dsw.*C_pow);
    C_poo=phi*(1-Sw)./dt.*(Ct./Bo);
    Pcow=ppval(sw_pc,Sw);
    alpha=-C_swo./C_sww;

    %% injection rate and boundary conditions
    Q_water(y_well_inj,x_well_inj)=injection;
      
    %% calculating transmissiblity for the Y direction.
    for i=1:Nx
         for j=1:Ny-1
            if P(j+1,i)>P(j,i)
                TOY(j+1,i)=unit_mod1*unit_mod2*(2*Landa_O(j+1,i))/(dy*((dy/k(j,i))+(dy/(k(j+1,i)))));
                TWY(j+1,i)=unit_mod1*unit_mod2*(2*Landa_W(j+1,i))/(dy*((dy/k(j,i))+(dy/(k(j+1,i)))));
            else
                TOY(j+1,i)=unit_mod1*unit_mod2*(2*Landa_O(j,i))/(dy*((dy/k(j,i))+(dy/(k(j+1,i)))));
                TWY(j+1,i)=unit_mod1*unit_mod2*(2*Landa_W(j,i))/(dy*((dy/k(j,i))+(dy/(k(j+1,i)))));
            end
         end
    end
    for i=1:Nx
        j=1;
        TOY(j,i)=unit_mod1*unit_mod2*(2*Landa_O(j,i))/(dy*((dy/k(j,i))+(dy/(k(j,i)))));
        TWY(j,i)=unit_mod1*unit_mod2*(2*Landa_W(j,i))/(dy*((dy/k(j,i))+(dy/(k(j,i)))));
    end

    %% calculating transmissibility for the X direction.
    for i=1:Nx-1
         for j=1:Ny
            if P(j,i+1)>P(j,i)
                TOX(j,i+1)=unit_mod1*unit_mod2*(2*Landa_O(j,i+1))/(dx*((dx/k(j,i))+(dx/(k(j,i+1)))));
                TWX(j,i+1)=unit_mod1*unit_mod2*(2*Landa_W(j,i+1))/(dx*((dx/k(j,i))+(dx/(k(j,i+1)))));
            else
                TOX(j,i+1)=unit_mod1*unit_mod2*(2*Landa_O(j,i))/(dx*((dx/k(j,i))+(dx/(k(j,i+1)))));
                TWX(j,i+1)=unit_mod1*unit_mod2*(2*Landa_W(j,i))/(dx*((dx/k(j,i))+(dx/(k(j,i+1)))));
            end
         end
    end 
    for j=1:Ny
        i=1;
        TOX(j,i)=unit_mod1*unit_mod2*(2*Landa_O(j,i))/(dx*((dx/k(j,i))+(dx/(k(j,i)))));
        TWX(j,i)=unit_mod1*unit_mod2*(2*Landa_W(j,i))/(dx*((dx/k(j,i))+(dx/(k(j,i)))));
    end

    %% calculating Coefficient Matrix
    for m=1:1:Nx*Ny
        i=mod(m,Nx);
        if i==0,i=Nx;end;

        j=floor((m-i)/Nx)+1;


        if j==1
            A(m,m+Nx)=TOY(j+1,i)+alpha(j,i)*TWY(j+1,i);   %N(i,j);

            if i==1
                A(m,m+1)=TOX(j,i+1)+alpha(j,i)*TWX(j,i+1);   %E(i,j);
                A(m,m)=-( TOX(j,i+1)+TOY(j+1,i)+C_poo(j,i)  + WCi(j,i)/(dx*dy*h)*Landa_O(j,i) )-alpha(j,i)*(TWX(j,i+1)+TWY(j+1,i)+C_pow(j,i)  + WCi(j,i)/(dx*dy*h)*Landa_W(j,i) );
                
                b(m)=-(C_poo(j,i)+alpha(j,i)*C_pow(j,i))*P(j,i)...
                - WCi(j,i)/(dx*dy*h) * Landa_O(j,i) * P_bh - alpha(j,i)* WCi(j,i)/(dx*dy*h)*Landa_W(j,i)* P_bh ...
                + alpha(j,i)*TWX(j,i+1)*(Pcow(j,i+1)-Pcow(j,i))...
                + alpha(j,i)*TWY(j+1,i)*(Pcow(j+1,i)-Pcow(j,i));
            
            elseif i==Nx
                A(m,m-1)=TOX(j,i)+alpha(j,i)*TWX(j,i);   %W(i,j);
                A(m,m)=-(TOX(j,i)+TOY(j+1,i)+C_poo(j,i)  + WCi(j,i)/(dx*dy*h)*Landa_O(j,i) )-alpha(j,i)*(TWX(j,i)+TWY(j+1,i)+C_pow(j,i)  + WCi(j,i)/(dx*dy*h)*Landa_W(j,i) );
                
                b(m)=-(C_poo(j,i)+alpha(j,i)*C_pow(j,i))*P(j,i)...
                - WCi(j,i)/(dx*dy*h) * Landa_O(j,i) * P_bh - alpha(j,i)* WCi(j,i)/(dx*dy*h)*Landa_W(j,i)* P_bh ...
                - alpha(j,i)*TWX(j,i)*(Pcow(j,i)-Pcow(j,i-1))...
                + alpha(j,i)*TWY(j+1,i)*(Pcow(j+1,i)-Pcow(j,i));
            
            else
                A(m,m)=-(TOX(j,i)+TOX(j,i+1)+TOY(j+1,i)+C_poo(j,i))-alpha(j,i)*(TWX(j,i)+TWX(j,i+1)+TWY(j+1,i)+C_pow(j,i));%M(j,i)
                A(m,m+1)=TOX(j,i+1)+alpha(j,i)*TWX(j,i+1);   %E(i,j);
                A(m,m-1)=TOX(j,i)+alpha(j,i)*TWX(j,i);   %W(i,j);
                
                b(m)=-(C_poo(j,i)+alpha(j,i)*C_pow(j,i))*P(j,i)+alpha(j,i)*Q_water(j,i) ...
                + alpha(j,i)*TWX(j,i+1)*(Pcow(j,i+1)-Pcow(j,i))-alpha(j,i)*TWX(j,i)*(Pcow(j,i)-Pcow(j,i-1))...
                + alpha(j,i)*TWY(j+1,i)*(Pcow(j+1,i)-Pcow(j,i));
            end

        end

        if (j>1) && (j<Ny)
            A(m,m-Nx)=TOY(j,i)+alpha(j,i)*TWY(j,i);     %S(i,j);
            A(m,m+Nx)=TOY(j+1,i)+alpha(j,i)*TWY(j+1,i);     %N(i,j);

            if i==1 
                A(m,m+1)=TOX(j,i+1)+alpha(j,i)*TWX(j,i+1);   %E(i,j);
                A(m,m)=-(TOX(j,i+1)+TOY(j,i)+TOY(j+1,i)+C_poo(j,i))-alpha(j,i)*(TWX(j,i+1)+TWY(j,i)+TWY(j+1,i)+C_pow(j,i));%M(j,i)
                
                b(m)=-(C_poo(j,i)+alpha(j,i)*C_pow(j,i))*P(j,i)+alpha(j,i)*Q_water(j,i) ...
                + alpha(j,i)*TWX(j,i+1)*(Pcow(j,i+1)-Pcow(j,i))...
                + alpha(j,i)*TWY(j+1,i)*(Pcow(j+1,i)-Pcow(j,i))-alpha(j,i)*TWY(j,i)*(Pcow(j,i)-Pcow(j-1,i));
            elseif i==Nx
                A(m,m-1)=TOX(j,i)+alpha(j,i)*TWX(j,i);   %W(i,j);
                A(m,m)=-(TOX(j,i)+TOY(j,i)+TOY(j+1,i)+C_poo(j,i))-alpha(j,i)*(TWX(j,i)+TWY(j,i)+TWY(j+1,i)+C_pow(j,i));%M(j,i)
                
                b(m)=-(C_poo(j,i)+alpha(j,i)*C_pow(j,i))*P(j,i)...
                -alpha(j,i)*TWX(j,i)*(Pcow(j,i)-Pcow(j,i-1))...
                + alpha(j,i)*TWY(j+1,i)*(Pcow(j+1,i)-Pcow(j,i))-alpha(j,i)*TWY(j,i)*(Pcow(j,i)-Pcow(j-1,i));
            else 
                A(m,m+1)=TOX(j,i+1)+alpha(j,i)*TWX(j,i+1);   %E(i,j);
                A(m,m-1)=TOX(j,i)+alpha(j,i)*TWX(j,i);   %W(i,j);
                A(m,m)=-(TOX(j,i)+TOX(j,i+1)+TOY(j,i)+TOY(j+1,i)+C_poo(j,i))-alpha(j,i)*(TWX(j,i)+TWX(j,i+1)+TWY(j,i)+TWY(j+1,i)+C_pow(j,i));%M(j,i)
                
                b(m)=-(C_poo(j,i)+alpha(j,i)*C_pow(j,i))*P(j,i)+alpha(j,i)*Q_water(j,i) ...
                + alpha(j,i)*TWX(j,i+1)*(Pcow(j,i+1)-Pcow(j,i))-alpha(j,i)*TWX(j,i)*(Pcow(j,i)-Pcow(j,i-1))...
                + alpha(j,i)*TWY(j+1,i)*(Pcow(j+1,i)-Pcow(j,i))-alpha(j,i)*TWY(j,i)*(Pcow(j,i)-Pcow(j-1,i));
            end
        end

        if j==Ny
            A(m,m-Nx)=TOY(j,i)+alpha(j,i)*TWY(j,i);     %S(i,j);

            if i==1 
                A(m,m+1)=TOX(j,i+1)+alpha(j,i)*TWX(j,i+1);   %E(i,j);
                A(m,m)=-(TOX(j,i+1)+TOY(j,i)+C_poo(j,i)  + WCi(j,i)/(dx*dy*h)*Landa_O(j,i) )-alpha(j,i)*(TWX(j,i+1)+TWY(j,i)+C_pow(j,i)  + WCi(j,i)/(dx*dy*h)*Landa_W(j,i) );%M(j,i)
                
                b(m)=-(C_poo(j,i)+alpha(j,i)*C_pow(j,i))*P(j,i)...
                - WCi(j,i)/(dx*dy*h) * Landa_O(j,i) * P_bh - alpha(j,i)* WCi(j,i)/(dx*dy*h)*Landa_W(j,i)* P_bh ...
                + alpha(j,i)*TWX(j,i+1)*(Pcow(j,i+1)-Pcow(j,i))...
                -alpha(j,i)*TWY(j,i)*(Pcow(j,i)-Pcow(j-1,i));
            
            elseif i==Nx
                A(m,m-1)=TOX(j,i)+alpha(j,i)*TWX(j,i);   %W(i,j);
                A(m,m)=-(TOX(j,i)+TOY(j,i)+C_poo(j,i)+  WCi(j,i)/(dx*dy*h)*Landa_O(j,i) )-alpha(j,i)*(TWX(j,i)+TWY(j,i)+C_pow(j,i) + WCi(j,i)/(dx*dy*h)*Landa_W(j,i) );
                
                b(m)=-(C_poo(j,i)+alpha(j,i)*C_pow(j,i))*P(j,i)...
                - WCi(j,i)/(dx*dy*h) * Landa_O(j,i) * P_bh - alpha(j,i)* WCi(j,i)/(dx*dy*h)*Landa_W(j,i)* P_bh ...
                -alpha(j,i)*TWX(j,i)*(Pcow(j,i)-Pcow(j,i-1))...
                -alpha(j,i)*TWY(j,i)*(Pcow(j,i)-Pcow(j-1,i));
            
            else 
                 A(m,m+1)=TOX(j,i+1)+alpha(j,i)*TWX(j,i+1);   %E(i,j);
                 A(m,m-1)=TOX(j,i)+alpha(j,i)*TWX(j,i);   %W(i,j);
                 A(m,m)=-(TOX(j,i)+TOX(j,i+1)+TOY(j,i)+C_poo(j,i))-alpha(j,i)*(TWX(j,i)+TWX(j,i+1)+TWY(j,i)+C_pow(j,i));
                 
                 b(m)=-(C_poo(j,i)+alpha(j,i)*C_pow(j,i))*P(j,i)...
                 + alpha(j,i)*TWX(j,i+1)*(Pcow(j,i+1)-Pcow(j,i))-alpha(j,i)*TWX(j,i)*(Pcow(j,i)-Pcow(j,i-1))...
                 -alpha(j,i)*TWY(j,i)*(Pcow(j,i)-Pcow(j-1,i));
             
            end
        end
    end
    P_t=A^-1*b;
    P_t=reshape(P_t,Ny,Nx);
    P_t=P_t';

    %% Explicit saturation calculation
    for j=2:Ny-1
        for i=2:Nx-1
            Sw_t(j,i)=Sw(j,i)+1/C_swo(j,i)* ( TOX(j,i+1)*(P_t(j,i+1)-P_t(j,i))-TOX(j,i)*(P_t(j,i)-P_t(j,i-1))+TOY(j+1,i)*(P_t(j+1,i)-P_t(j,i))+...
                -TOY(j,i)*(P_t(j,i)-P_t(j-1,i))-Q_oil(j,i)-C_poo(j,i)*(P_t(j,i)-P(j,i)));
        end
    end

    for i=2:Nx-1
        j=1;

        Sw_t(j,i)=Sw(j,i)+1/C_swo(j,i)* ( TOX(j,i+1)*(P_t(j,i+1)-P_t(j,i))-TOX(j,i)*(P_t(j,i)-P_t(j,i-1))+TOY(j+1,i)*(P_t(j+1,i)-P_t(j,i))+...
            -Q_oil(j,i)-C_poo(j,i)*(P_t(j,i)-P(j,i)));
    end

    for j=2:Ny-1
        i=1;

        Sw_t(j,i)=Sw(j,i)+1/C_swo(j,i)* ( TOX(j,i+1)*(P_t(j,i+1)-P_t(j,i))+TOY(j+1,i)*(P_t(j+1,i)-P_t(j,i))+...
            -TOY(j,i)*(P_t(j,i)-P_t(j-1,i))-Q_oil(j,i)-C_poo(j,i)*(P_t(j,i)-P(j,i)));
    end

    for j=2:Ny-1
        i=Nx;
        Sw_t(j,i)=Sw(j,i)+1/C_swo(j,i)* ( -TOX(j,i)*(P_t(j,i)-P_t(j,i-1))+TOY(j+1,i)*(P_t(j+1,i)-P_t(j,i))+...
            -TOY(j,i)*(P_t(j,i)-P_t(j-1,i))-Q_oil(j,i)-C_poo(j,i)*(P_t(j,i)-P(j,i)));
    end

    for i=2:Nx-1
        j=Ny;
        Sw_t(j,i)=Sw(j,i)+1/C_swo(j,i)* ( TOX(j,i+1)*(P_t(j,i+1)-P_t(j,i))-TOX(j,i)*(P_t(j,i)-P_t(j,i-1))+...
            -TOY(j,i)*(P_t(j,i)-P_t(j-1,i))-Q_oil(j,i)-C_poo(j,i)*(P_t(j,i)-P(j,i)));
    end

    for j=1
        for i=1
            Sw_t(j,i)=Sw(j,i)+1/C_swo(j,i)* ( TOX(j,i+1)*(P_t(j,i+1)-P_t(j,i))+TOY(j+1,i)*(P_t(j+1,i)-P_t(j,i))...
                -WCi(j,i)/(dx*dy*h)*Landa_O(j,i)*(P_t(j,i)-P_bh)-C_poo(j,i)*(P_t(j,i)-P(j,i)));
        end
    end

    for j=Ny
        for i=Nx
            Sw_t(j,i)=Sw(j,i)+1/C_swo(j,i)* (-TOX(j,i)*(P_t(j,i)-P_t(j,i-1))...
                -TOY(j,i)*(P_t(j,i)-P_t(j-1,i))-WCi(j,i)/(dx*dy*h)*Landa_O(j,i)*(P_t(j,i)-P_bh)-C_poo(j,i)*(P_t(j,i)-P(j,i)));
        end
    end

    for j=Ny
        for i=1
            Sw_t(j,i)=Sw(j,i)+1/C_swo(j,i)* ( TOX(j,i+1)*(P_t(j,i+1)-P_t(j,i))...
                -TOY(j,i)*(P_t(j,i)-P_t(j-1,i))-WCi(j,i)/(dx*dy*h)*Landa_O(j,i)*(P_t(j,i)-P_bh)-C_poo(j,i)*(P_t(j,i)-P(j,i)));
        end
    end

    for j=1
        for i=Nx
            Sw_t(j,i)=Sw(j,i)+1/C_swo(j,i)* ( -TOX(j,i)*(P_t(j,i)-P_t(j,i-1))+TOY(j+1,i)*(P_t(j+1,i)-P_t(j,i))...
                -WCi(j,i)/(dx*dy*h)*Landa_O(j,i)*(P_t(j,i)-P_bh) -C_poo(j,i)*(P_t(j,i)-P(j,i)));
        end
    end 
    
    P=P_t;
    Sw=Sw_t;
  
%% simulating water advance 

%     figure(1);
%     pcolor([1:1:Nx]*dx,[1:1:Ny]*dy,1-Sw); 
%     shading interp
%     colormap jet
%     colorbar
%     xlabel('Horizotal Direction(ft)');
%     ylabel('Vertical Direction(ft)');
%     title('Oil Saturation');
%     pause(0.01);
%     hold on

%% saving attained resutls 
    Result{1,Iteration}=Sw;
    Result{2,Iteration}=P;
    Result{3,Iteration}=Landa_O;
    Result{4,Iteration}=Landa_W;
    
    %% stopping the program at high water production
    x_well=1;y_well=1;
    if Landa_W(x_well,y_well)/(Landa_O(x_well,y_well)+Landa_W(x_well,y_well))>0.95,break;end
    
    Iteration=Iteration+1;
end

save Output Result

%% showing other results
close all
Landa_W=Result(4,:);
Landa_O=Result(3,:);
P=Result(2,:);
Sw=Result(1,:);
n=Iteration-1;
Cummulative_oil=0;
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)]);plot(1,1)
for i=1:n
    subplot(3,2,1)
    water_cut=Landa_W{i}(x_well,y_well)/(Landa_W{i}(x_well,y_well)+Landa_O{i}(x_well,y_well));   
    plot(i*dt,water_cut,'.')
    xlabel('Days');ylabel('Water Cut');
    title('Water Cut');
    hold on;
                    %Oil production rate     
    subplot(3,2,2) 
    Oil_rate=WCi(x_well,y_well)*(Landa_O{i}(x_well,y_well)*(P{i}(x_well,y_well)-P_bh))/unit_mod1;  
    plot(i*dt,Oil_rate,'.')
    xlabel('Days');ylabel('Oil production (STB/day)');
    title('Oil Production Rate vs. Time')
    hold on
                %cummulative produced oil
    subplot(3,2,4)
    Cummulative_oil=Cummulative_oil+Oil_rate;
    plot(i*dt,Cummulative_oil,'.')
    xlabel('Days');ylabel('Cummulative Oil production (STB)');
    title('Cummulative Oil Production rate vs. Time')
    hold on
              %Water/Oil Production Curve
    subplot(3,2,3)
    water_oil_ratio=Landa_W{i}(x_well,y_well)/Landa_O{i}(x_well,y_well);   
    plot(i*dt,water_oil_ratio,'.')
    xlabel('Days');ylabel('Water/Oil Production Ratio');
    title('Water/Oil Production Ratio')
    hold on;
            %Water Injection Downhole Pressure
    subplot(3,2,5)
    water_P=P{i}(y_well_inj,x_well_inj);  
    plot(i*dt,water_P-Q_water(y_well_inj,x_well_inj)*h*dx*dy/(WCi(y_well_inj,x_well_inj)*...
        Landa_W{i}(y_well_inj,x_well_inj)),'.')
    xlabel('Days');ylabel('Water Injection Pressure');
    title('Water Injection Pressure')
    hold on;
    subplot(3,2,6)
    [M_Ny M_Nx]=meshgrid([1:Ny],[1:Nx]);
     mesh(M_Ny,M_Nx,Sw{i})
    title('Sw distribution')
    xlabel('x');ylabel('y');
    zlabel('Sw');
    pause(0.001)
end
Recover_factor = Cummulative_oil*unit_mod1/(phi*Length^2*(1-0.25)*h)