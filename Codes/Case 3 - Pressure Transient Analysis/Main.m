% Reading the input data
% using importfile.m function to read the data fron .txt file
Input_cell=importfile('Input.txt');   
FIO = fopen('creatIO.m','w');
for i=1:numel(Input_cell)
    fwrite(FIO,Input_cell{i});
    fprintf(FIO,'\n');
end
creatIO;
close all; clc;


%% unit conversion
if strcmp(unit,'field')
    unit_mod1 = 5.615; %STB to ft3
    unit_mod2 = 1.127e-3; %Field unit for darcy equation.
end


%% Well production
Production=+Qo/(dx*dy*h)*unit_mod1;

%% processing input data
dx=Length/Nx;
dy=Width/Ny;
h=Height;     %ft

Bo_0=ones(Nx).*Bo_0;
k = ones(Nx)*k;

%% making variables
TOY=zeros(Ny,Nx);
TOX=zeros(Ny,Nx);
A=zeros(Nx*Ny);
Q_oil=zeros(Ny,Nx);
b=zeros(Ny*Nx,1);
kro=zeros(Ny,Nx);
Result=cell(4,floor(Max_days/dt));

%% Defining landa and flow rate.
Iteration=1;
P_t=P;
while Iteration<floor(Max_days/dt)

    
    %% Build up test 
    if Iteration*dt>shut_in
        Production=0;  %STB/day/ft3*5.615=1/day  %injection must be minus
    end
    
    
    %% Calculating Mobility of each grid
      Bo=Bo_0; 
      kro=ones(Ny,Nx);
      Landa_O=kro./(Muo.*Bo);

    %% Defining parameters
    C_poo=phi./dt.*(Ct./Bo);
    
    %% Proudction rate and Boundary conditions 
    Q_oil(y_well,x_well)=Production;
    
    %% calculating transmissiblity for the Y direction.
    for i=1:Nx
         for j=1:Ny-1
            if P(j+1,i)>P(j,i)
                TOY(j+1,i)=unit_mod1*unit_mod2*(2*Landa_O(j+1,i))/(dy*((dy/k(j,i))+(dy/(k(j+1,i)))));
            else
                TOY(j+1,i)=unit_mod1*unit_mod2*(2*Landa_O(j,i))/(dy*((dy/k(j,i))+(dy/(k(j+1,i)))));
            end
         end
    end
    for i=1:Nx
        j=1;
        TOY(j,i)=unit_mod1*unit_mod2*(2*Landa_O(j,i))/(dy*((dy/k(j,i))+(dy/(k(j,i)))));
    end

    %% calculating transmissibility for the X direction.
    for i=1:Nx-1
         for j=1:Ny
            if P(j,i+1)>P(j,i)
                TOX(j,i+1)=unit_mod1*unit_mod2*(2*Landa_O(j,i+1))/(dx*((dx/k(j,i))+(dx/(k(j,i+1)))));
            else
                TOX(j,i+1)=unit_mod1*unit_mod2*(2*Landa_O(j,i))/(dx*((dx/k(j,i))+(dx/(k(j,i+1)))));
            end
         end
    end 
    for j=1:Ny
        i=1;
        TOX(j,i)=unit_mod1*unit_mod2*(2*Landa_O(j,i))/(dx*((dx/k(j,i))+(dx/(k(j,i)))));
    end

    %% Implicit pressure calculation
    for m=1:1:Nx*Ny
        i=mod(m,Nx);
        if i==0,i=Nx;end;

        j=floor((m-i)/Nx)+1;


        if j==1
            A(m,m+Nx)=TOY(j+1,i);

            if i==1
                A(m,m+1)=TOX(j,i+1);                                                %E(i,j);
                A(m,m)=-( TOX(j,i+1)+TOY(j+1,i)+C_poo(j,i) );
                
                b(m)=-(C_poo(j,i))*P(j,i);
            
            elseif i==Nx
                A(m,m-1)=TOX(j,i);                                                  %W(i,j);
                A(m,m)=-(TOX(j,i)+TOY(j+1,i)+C_poo(j,i) );
                
                b(m)=-(C_poo(j,i))*P(j,i);
            else
                A(m,m)=-(TOX(j,i)+TOX(j,i+1)+TOY(j+1,i)+C_poo(j,i));                %M(j,i)
                A(m,m+1)=TOX(j,i+1);                                                %E(i,j);
                A(m,m-1)=TOX(j,i);                                                  %W(i,j);
                
                b(m)=-(C_poo(j,i))*P(j,i)+Q_oil(j,i);
            end

        end

        if (j>1) && (j<Ny)
            A(m,m-Nx)=TOY(j,i);                                                      %S(i,j);
            A(m,m+Nx)=TOY(j+1,i);                                                    %N(i,j);

            if i==1 
                A(m,m+1)=TOX(j,i+1);     %E(i,j);
                A(m,m)=-(TOX(j,i+1)+TOY(j,i)+TOY(j+1,i)+C_poo(j,i));     %M(j,i)
                
                b(m)=-(C_poo(j,i))*P(j,i)+Q_oil(j,i);
            elseif i==Nx
                A(m,m-1)=TOX(j,i);    %W(i,j);
                A(m,m)=-(TOX(j,i)+TOY(j,i)+TOY(j+1,i)+C_poo(j,i)); %M(j,i)
                
                b(m)=-(C_poo(j,i))*P(j,i);
                
            else 
                A(m,m+1)=TOX(j,i+1);    %E(i,j);
                A(m,m-1)=TOX(j,i);    %W(i,j);
                A(m,m)=-(TOX(j,i)+TOX(j,i+1)+TOY(j,i)+TOY(j+1,i)+C_poo(j,i));     %M(j,i)
                
                b(m)=-(C_poo(j,i))*P(j,i)+Q_oil(j,i);
            end
        end

        if j==Ny
            A(m,m-Nx)=TOY(j,i);     %S(i,j);

            if i==1 
                A(m,m+1)=TOX(j,i+1);   %E(i,j);
                A(m,m)=-(TOX(j,i+1)+TOY(j,i)+C_poo(j,i) );   %M(j,i)
                
                b(m)=-(C_poo(j,i))*P(j,i);
                            
            elseif i==Nx
                A(m,m-1)=TOX(j,i);           %W(i,j);
                A(m,m)=-(TOX(j,i)+TOY(j,i)+C_poo(j,i));
                
                b(m)=-(C_poo(j,i))*P(j,i);
                            
            else 
                 A(m,m+1)=TOX(j,i+1);   %E(i,j);
                 A(m,m-1)=TOX(j,i);    %W(i,j);
                 A(m,m)=-(TOX(j,i)+TOX(j,i+1)+TOY(j,i)+C_poo(j,i));                
                 b(m)=-(C_poo(j,i))*P(j,i);         
             
            end
        end
    end
    P_t=A^-1*b;
    P_t=reshape(P_t,Ny,Nx);
    P_t=P_t';
    P=P_t;

    Result{2,Iteration}=P;
    Result{3,Iteration}=Landa_O;
        
    Iteration=Iteration+1;
end

save Output Result

%% Results
close all
Landa_O=Result(3,:);
P=Result(2,:);
n=Iteration-1;

figure(1)
for i=1:n
      t1(i)=i*dt;
      p1(i)=P{i}(floor(Ny/2),floor(Nx/2));
end
plot(t1,p1,'*-')
grid on;
title(['Shuting well after ' num2str(shut_in) ' days flow' ]);
xlabel('Time (days)');ylabel('pressure (psi)');


figure(2)
for i=1:n
    if i*dt<shut_in
        time(i)=i*dt;
        pressure(i)=P{i}(floor(Ny/2),floor(Nx/2));
    end
end
A=polyfit(log10(time(3:10)),pressure(3:10),1);
semilogx(time,pressure,'*b',[1e-2 time(3:10) 1e0],A(1)*log10([1e-2 time(3:10) 1e0])+A(2),'r-');
xlabel('\Deltat  (days)')
ylabel('Pressure (Psi)');
title(' Draw Down Test');
text(1e-2,2000,['P = ' num2str(A(1)) ' log10(\Delta t) + ' num2str(A(2))]);
grid on

figure(3)
for i=1:n
    if i*dt>shut_in
        t3(i)=(i*dt+shut_in)/(i*dt);
        p3(i)=P{i}(floor(Ny/2),floor(Nx/2));
        hold on;
    end
end
t3(t3==0)=[];p3(p3==0)=[];
A=polyfit(log10(t3(end-35:end-30)),p3(end-35:end-30),1);
semilogx(t3,p3,'*b',[1.7 t3(end-35:end-30) 2],A(1)*log10( [1.7 t3(end-35:end-30) 2])+A(2),'r-','linewidth',3);
xlabel('Horner time  (days)')
ylabel('Pressure (Psi)');
title('Build up test');
text(1.7,2000,['P = ' num2str(A(1)) ' log10(t_{horner}) + ' num2str(A(2))]);
grid on