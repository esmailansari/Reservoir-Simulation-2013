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

%% defining parameters
dx =  Length/Nx;
dy = Width/Ny;
h  =  Height/Nz;

%unit
if strcmp(unit,'field')
    unit_mod1 = 5.615; %STB to ft3
    unit_mod2 = 1.127e-3; %Field unit for darcy equation.
end

%Corey equation 
Sw_range = S_wc:0.01:1-S_or;
kro=kro1*((1-Sw_range-S_or)/(1-S_wc-S_or)).^a;
krw=krw1*((Sw_range-S_wc)/(1-S_wc-S_or)).^b;

%   Sw   krw  kro   Pc
%   defining Pc as zero. 
Data=[Sw_range' krw' kro' zeros(numel(Sw_range),1)];

% linear interpolation for unknown points
sw_krw = interp1(Data(:,1),Data(:,2),'linear','pp');
sw_kro = interp1(Data(:,1),Data(:,3),'linear','pp');
sw_pc = interp1(Data(:,1),Data(:,4),'linear','pp');
Domain = S_wc:0.005:1;
%calculating Dpc/Dsw. For buckley leverett case this will be zero. 
DpcDsw = diff(ppval(sw_pc,Domain))./diff(Domain);
sw_dpc_dsw = interp1((Domain(1:end-1)),DpcDsw,'linear','pp');

%% well schedule
injection=-Qw/(dx*dy*h)* unit_mod1;  %STB/day/ft3*5.615=1/day  %injection must be minus.
Q_water=zeros(Nx,1);Q_oil=zeros(Nx,1);
Q_water(x_well_inj)=injection;
Q_oil(x_well_pro)  =-injection;

% These data are not useful but the variables exist in the robust formulation. 
% so putting them to zero. 
P_bh=0;
WCi=zeros(1,Nx); 
Rw=0.33; %ft
Re=sqrt(dy*dx/pi);

%% defining the dt for numerical buckley leverett.
dt= (dy*h)*Length*phi/(Qw*unit_mod1)*t_d;

%% making variables
TOX=zeros(Nx,1);TWX=zeros(Nx,1);
A=zeros(Nx,1);
Sw_t=zeros(Nx,1);
b=zeros(Nx,1);
Pcow=zeros(Nx,1);
kro=zeros(Nx,1);
krw=zeros(Nx,1);
Result=cell(4,floor(Max_days/dt));

%% Defining landa and flow rate.
Iteration=1;
Ct=1e-17; 
P_t=P;
while Iteration<floor(Max_days/dt)
    
    %% Calculating Mobility of each grid
      % incompressible fluid.
      Bw=Bw_0;
      Bo=Bo_0; 
      % interpolation
      kro=ppval(sw_kro,Sw);
      krw=ppval(sw_krw,Sw);
      %putting interpolation endpoints to zero. 
      kro(kro>Data(end,1))=Data(end,1);kro(kro<0)=0;   
      krw(krw>Data(end,1))=Data(end,1);krw(krw<0)=0;  
      % mobility
      Landa_O=kro./(Muo.*Bo);
      Landa_W=krw./(Muw.*Bw);
      
    %% Defining parameters
    C_swo=-phi./(Bo.*dt).*ones(Nx,1);
    C_pow=phi.*Sw./dt.*(Ct./Bo);
    dpc_dsw=ppval(sw_dpc_dsw,Sw);
    C_sww=(phi./(Bw.*dt).*ones(Nx,1)-dpc_dsw.*C_pow); 
    C_poo=phi*(1-Sw)./dt.*(Ct./Bo);
    Pcow=ppval(sw_pc,Sw);
    alpha=-C_swo./C_sww;

    %% calculating transmissibility for the X direction. (upstream weighting).
    for i=1:Nx-1
            if P(i+1)>P(i)
                TOX(i)=unit_mod1*unit_mod2*(2*Landa_O(i+1))/(dx*((dx/k(i))+(dx/(k(i+1)))));
                TWX(i)=unit_mod1*unit_mod2*(2*Landa_W(i+1))/(dx*((dx/k(i))+(dx/(k(i+1)))));
            else
                TOX(i+1)=unit_mod1*unit_mod2*(2*Landa_O(i))/(dx*((dx/k(i))+(dx/(k(i+1)))));
                TWX(i+1)=unit_mod1*unit_mod2*(2*Landa_W(i))/(dx*((dx/k(i))+(dx/(k(i+1)))));
            end
    end
    
    for i=1
                TOX(i)=unit_mod1*unit_mod2*(2*Landa_O(i))/(dx*((dx/k(i))+(dx/(k(i+1)))));
                TWX(i)=unit_mod1*unit_mod2*(2*Landa_W(i))/(dx*((dx/k(i))+(dx/(k(i+1)))));
    end

    %% Implicit pressure calculation
    for m=1:Nx
        i=m;
            if m==1 
                A(m,m+1) =  + TOX(i+1)+alpha(i)*TWX(i+1);   
                A(m,m)   =  -( TOX(i+1)+C_poo(i) )+ ...
                            - alpha(i)*(TWX(i+1)+C_pow(i) ...
                            + WCi(i)/(dx*dy*h)*Landa_W(i))+ ... 
                            + WCi(i)/(dx*dy*h)*Landa_O(i);
                             
                b(m)     =  - ( C_poo(i)+alpha(i)*C_pow(i) ) * P(i)...
                            + alpha(i)*TWX(i+1)*(Pcow(i+1)-Pcow(i))...
                            + Q_oil(i) + alpha(i)*Q_water(i)+ ...
                            - WCi(i)/(dx*dy*h) * Landa_O(i) * P_bh + ...
                            - alpha(i)* WCi(i)/(dx*dy*h)*Landa_W(i)* P_bh;
            
            elseif m==Nx
                A(m,m-1) =    + TOX(i)+alpha(i)*TWX(i);  
                A(m,m)   =    - ( TOX(i)+C_poo(i) ) + ...
                              -alpha(i)*( TWX(i)+C_pow(i) ) + ...
                              + WCi(i)/(dx*dy*h)*Landa_O(i);
                
                b(m)     =     - ( C_poo(i)+alpha(i)*C_pow(i) )*P(i)...
                               - alpha(i)*TWX(i)*( Pcow(i)-Pcow(i-1) )...
                               + Q_oil(i) + alpha(i)*Q_water(i) + ...
                               - WCi(i)/(dx*dy*h) * Landa_O(i) * P_bh ... 
                               - alpha(i)* WCi(i)/(dx*dy*h)*Landa_W(i)* P_bh;
            
            else 
                 A(m,m+1) =    TOX(i+1)+alpha(i)*TWX(i+1);   
                 A(m,m-1) =    TOX(i)+alpha(i)*TWX(i);   
                 A(m,m)   =  - ( TOX(i)+TOX(i+1)+C_poo(i) )+...
                             - alpha(i)*(  TWX(i)+TWX(i+1)+C_pow(i)  );
                 
                 b(m) = - ( C_poo(i)+alpha(i)*C_pow(i))*P(i)...
                        +   alpha(i)*TWX(i+1)*(Pcow(i+1)-Pcow(i) ) ...
                        -   alpha(i)*TWX(i)*(  Pcow(i)-Pcow(i-1)  )...
                        +   Q_oil(i) + alpha(i)*Q_water(i) + ...
                        -   WCi(i)/(dx*dy*h) * Landa_O(i) * P_bh + ...
                        -   alpha(i)* WCi(i)/(dx*dy*h)*Landa_W(i)* P_bh;
             
            end
    end

    P_t=A^-1*b;

    %% Explicit saturation calculation
    for i=1;
        Sw_t(i) = 1-S_or;
    end

    for i=Nx
        Sw_t(i) = S_wc;
    end

    for i=2:Nx-1
        Sw_t(i) = Sw(i)+1/C_swo(i)*(TOX(i+1)*(P_t(i+1)-P_t(i))-TOX(i)*(P_t(i)-P_t(i-1))+ ...
                  -Q_oil(i)-WCi(i)/(dx*dy*h)*Landa_O(i)*(P_t(i)-P_bh)-C_poo(i)*(P_t(i)-P(i)));
    end
    
    P=P_t;
    Sw=Sw_t;
  
%% saving attained resutls 
    Result{1,Iteration}=Sw;
    Result{2,Iteration}=P;
    Result{3,Iteration}=Landa_O;
    Result{4,Iteration}=Landa_W;
    
    Iteration=Iteration+1;
end

save Output Result

%% Solving Buckley leverett Analytically
global swi sor muo 
swi=S_wc;sor=S_or;muw=Muw;muo=Muo;
options = optimset('display','off');
s_wf    = fzero('fw',.5,options,swi,muw);

%% showing results
close all
Landa_W=Result(4,:);
Landa_O=Result(3,:);
P=Result(2,:);
Sw=Result(1,:); 
n=Iteration-1;
Cummulative_oil=0;
scrsz = get(0,'ScreenSize');
figure('Position',[1 1 scrsz(3) scrsz(4)]);plot(1,1)
time=t_d; 
step=5;
for i=1:step:n
    plot(linspace(0,1,Nx),Sw{i+1})
    axis([0 1 0 1])
    hold on
    x_wf=dfw(s_wf,muw)*time;
    plot([x_wf ,x_wf, 1],[s_wf, swi, swi],'r-');
    
    for sw=s_wf:.001:(1-sor)   
        plot(dfw(sw,muw)*time,sw,'r')
    end
    s_wb=fzero('fw',.5,options,0,muw);
    x_wb=dfw(s_wb,muw)*time;
    axis([0 1 0 1]);xlabel('x_w');ylabel('s_w');
    text(0.1*floor((i+1)/step),0.3,['t_D=' num2str(time)])
    pause(0.1)
    time = time + step * t_d; 
   legend('Numerical','Analytical') 
end
