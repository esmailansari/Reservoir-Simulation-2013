function [y,Landa_W,Landa_O] = EquCalc(X,varargin)
global Sw P dt

[rho_oil,rho_water,Data,sw_kro,sw_krw,sw_pc,sw_dpc_dsw,unit_mod1,unit_mod2 ...
 ,Length,dx,dy,Nx,Ny,k,phi,h,Bo_0,Bw_0,Muo,Muw,Ct...
 ,Q_water,Rw,Re,WCi,P_bh,Nz,dz]=varargin{:};
      
%% Distributing X between P, Sw
Sw_t=zeros(Nz,Ny,Nx);
P_t=zeros(Nz,Ny,Nx);
Count=1;
for kk=1:Nz
    for j=1:Ny
        for i=1:Nx

               P_t(kk,j,i)= X(Count);
               Sw_t(kk,j,i) = X(Count+1);
               Count = Count+2;

        end
    end
end

%% making variables
TOZ=zeros(Nz-1,Ny,Nx);TWZ=zeros(Nz-1,Ny,Nx);
TOY=zeros(Nz,Ny-1,Nx);TWY=zeros(Nz,Ny-1,Nx);
TOX=zeros(Nz,Ny,Nx-1);TWX=zeros(Nz,Ny,Nx-1);

%% Calculating Mobility of each grid 
kro=ppval(sw_kro,Sw_t);
krw=ppval(sw_krw,Sw_t);
Bo=Bo_0;
Bw=Bw_0;
Landa_O=kro./(Muo.*Bo);
Landa_W=krw./(Muw.*Bw);

%% Defining parameters
C_swo=-phi./(Bo.*dt).*ones(Nz,Ny,Nx);
C_pow=phi.*Sw_t./dt.*(Ct./Bo);
dpc_dsw=ppval(sw_dpc_dsw,Sw_t);
C_sww=(phi./(Bw.*dt).*ones(Nz,Ny,Nx)-dpc_dsw.*C_pow);  
C_poo=phi.*(1-Sw_t)./dt.*(Ct./Bo);
Pcow=ppval(sw_pc,Sw_t);
       
%% calculating transmissiblity for the Y direction.
for kk=1:Nz
    for i=1:Nx
         for j=1:Ny-1
                if P_t(kk,j+1,i)>P_t(kk,j,i)
                    TOY(kk,j,i)=unit_mod1*unit_mod2*(2*Landa_O(kk,j+1,i))/(dy*((dy/k(kk,j,i))+(dy/(k(kk,j+1,i)))));
                    TWY(kk,j,i)=unit_mod1*unit_mod2*(2*Landa_W(kk,j+1,i))/(dy*((dy/k(kk,j,i))+(dy/(k(kk,j+1,i)))));
                else
                    TOY(kk,j,i)=unit_mod1*unit_mod2*(2*Landa_O(kk,j,i))/(dy*((dy/k(kk,j,i))+(dy/(k(kk,j+1,i)))));
                    TWY(kk,j,i)=unit_mod1*unit_mod2*(2*Landa_W(kk,j,i))/(dy*((dy/k(kk,j,i))+(dy/(k(kk,j+1,i)))));
                end
         end
    end
end

%% calculating transmissibility for the X direction.
for kk=1:Nz
    for i=1:Nx-1
         for j=1:Ny
                if P_t(kk,j,i+1)>P_t(kk,j,i)
                    TOX(kk,j,i)=unit_mod1*unit_mod2*(2*Landa_O(kk,j,i+1))/(dx*((dx/k(kk,j,i))+(dx/(k(kk,j,i+1)))));
                    TWX(kk,j,i)=unit_mod1*unit_mod2*(2*Landa_W(kk,j,i+1))/(dx*((dx/k(kk,j,i))+(dx/(k(kk,j,i+1)))));
                else
                    TOX(kk,j,i)=unit_mod1*unit_mod2*(2*Landa_O(kk,j,i))/(dx*((dx/k(kk,j,i))+(dx/(k(kk,j,i+1)))));
                    TWX(kk,j,i)=unit_mod1*unit_mod2*(2*Landa_W(kk,j,i))/(dx*((dx/k(kk,j,i))+(dx/(k(kk,j,i+1)))));
                end
         end
    end 
end

%% calculating transmissiblity for the Z direction.
for kk=1:Nz-1
    for i=1:Nx
         for j=1:Ny
                if P_t(kk+1,j,i)>P_t(kk,j,i)
                    TOZ(kk,j,i)=unit_mod1*unit_mod2*(2*Landa_O(kk+1,j,i))/(dz*((dz/k(kk,j,i))+(dz/(k(kk+1,j,i)))));
                    TWZ(kk,j,i)=unit_mod1*unit_mod2*(2*Landa_W(kk+1,j,i))/(dz*((dz/k(kk,j,i))+(dz/(k(kk+1,j,i)))));
                else
                    TOZ(kk,j,i)=unit_mod1*unit_mod2*(2*Landa_O(kk,j,i))/(dz*((dz/k(kk,j,i))+(dz/(k(kk+1,j,i)))));
                    TWZ(kk,j,i)=unit_mod1*unit_mod2*(2*Landa_W(kk,j,i))/(dz*((dz/k(kk,j,i))+(dz/(k(kk+1,j,i)))));
                end
         end
    end
end

%% Adding Zeros to Transmissibity and Pressure
[NTX_z NTX_y NTX_x]=size(TOX);
TOX_dummy=zeros(NTX_z,NTX_y,NTX_x+2);
TOX_dummy(:,:,2:end-1)=TOX; TOX=TOX_dummy;
TWX_dummy=zeros(NTX_z,NTX_y,NTX_x+2);
TWX_dummy(:,:,2:end-1)=TWX; TWX=TWX_dummy;

[NTY_z NTY_y NTY_x]=size(TOY);
TOY_dummy=zeros(NTY_z,NTY_y+2,NTY_x);
TOY_dummy(:,2:end-1,:)=TOY; TOY=TOY_dummy;
TWY_dummy=zeros(NTY_z,NTY_y+2,NTY_x);
TWY_dummy(:,2:end-1,:)=TWY; TWY=TWY_dummy;

[NTZ_z NTZ_y NTZ_x]=size(TOZ);
TOZ_dummy=zeros(NTZ_z+2,NTZ_y,NTZ_x);
TOZ_dummy(2:end-1,:,:)=TOZ; TOZ=TOZ_dummy;
TWZ_dummy=zeros(NTZ_z+2,NTZ_y,NTZ_x);
TWZ_dummy(2:end-1,:,:)=TWZ; TWZ=TWZ_dummy;

[NP_z NP_y NP_x]=size(P_t);
P_t_dummy=zeros(NP_z+2,NP_y+2,NP_x+2);
P_t_dummy(2:end-1,2:end-1,2:end-1)=P_t;
P_t_dummy([1 end],2:end-1,2:end-1)=P_t([1 end],:,:);
P_t_dummy(2:end-1,2:end-1,[1 end])=P_t(:,:,[1 end]);
P_t_dummy(2:end-1,[1 end],2:end-1)=P_t(:,[1 end],:);
P_t=P_t_dummy;

[NP_z NP_y NP_x]=size(Pcow);
Pcow_dummy=zeros(NP_z+2,NP_y+2,NP_x+2);
Pcow_dummy(2:end-1,2:end-1,2:end-1)=Pcow;
Pcow_dummy([1 end],2:end-1,2:end-1)=Pcow([1 end],:,:);
Pcow_dummy(2:end-1,2:end-1,[1 end])=Pcow(:,:,[1 end]);
Pcow_dummy(2:end-1,[1 end],2:end-1)=Pcow(:,[1 end],:);
Pcow=Pcow_dummy;

%% Making Eq.s
y=zeros(2*Nz*Nx*Ny,1);
FF=zeros(2*Nz*Nx*Ny,1);
AA=zeros(2*Nz*Nx*Ny,1);
QQ=zeros(2*Nz*Nx*Ny,1);

for kk=1:Nz
    for j=1:Ny
        for i=1:Nx

         m = (kk-1)*Ny*Nx+(j-1)*Nx+i;

          FF(2*m-1)= TOX(kk,j,i+1)* ( P_t(kk+1,j+1,i+2)-P_t(kk+1,j+1,i+1) ) + ...
                    -TOX(kk,j,i) * ( P_t(kk+1,j+1,i+1)-P_t(kk+1,j+1,i) ) + ... 
                    TOY(kk,j+1,i)* ( P_t(kk+1,j+2,i+1)-P_t(kk+1,j+1,i+1) ) + ...
                    -TOY(kk,j,i) * ( P_t(kk+1,j+1,i+1)-P_t(kk+1,j,i+1) ) + ...
                    TOZ(kk+1,j,i)* ( P_t(kk+2,j+1,i+1)-P_t(kk+1,j+1,i+1) - 1/144*rho_oil/Bo(kk,j,i)*(kk-1/2)*dz ) + ...
                    -TOZ(kk,j,i) * ( P_t(kk+1,j+1,i+1)-P_t(kk,j+1,i+1) - 1/144*rho_oil/Bo(kk,j,i)*(kk-3/2)*dz );

          AA(2*m-1)=C_swo(kk,j,i)*( Sw(kk,j,i)-Sw_t(kk,j,i))+...
                     -C_poo(kk,j,i)*( P_t(kk+1,j+1,i+1)-P(kk,j,i)     );

          QQ(2*m-1)= -WCi(kk,j,i)/(dx*dy*dz) * (Landa_O(kk,j,i))*( P_t(kk+1,j+1,i+1)-P_bh(kk,j,i) );


          y(2*m-1) =  FF(2*m-1) + AA(2*m-1) + QQ(2*m-1);

     
          FF(2*m) = TWX(kk,j,i+1)*( (P_t(kk+1,j+1,i+2)-P_t(kk+1,j+1,i+1)) - ...
                                     (Pcow(kk+1,j+1,i+2)-Pcow(kk+1,j+1,i+1) ))+...
                     -TWX(kk,j,i) *( (P_t(kk+1,j+1,i+1)-P_t(kk+1,j+1,i)) - (Pcow(kk+1,j+1,i+1)-Pcow(kk+1,j+1,i) ))+...
                     TWY(kk,j+1,i)*( (P_t(kk+1,j+2,i+1)-P_t(kk+1,j+1,i+1)) - ...
                                     (Pcow(kk+1,j+2,i+1)-Pcow(kk+1,j+1,i+1) ))+...
                     -TWY(kk,j,i) *( (P_t(kk+1,j+1,i+1)-P_t(kk+1,j,i+1)) - (Pcow(kk+1,j+1,i+1)-Pcow(kk+1,j,i+1) ))+...
                     TWZ(kk+1,j,i)*( (P_t(kk+2,j+1,i+1)-P_t(kk+1,j+1,i+1) - 1/144*rho_water/Bw(kk,j,i)*(kk-1/2)*dz  ) - ...
                                      (Pcow(kk+2,j+1,i+1)-Pcow(kk+1,j+1,i+1) ))+...
                     -TWZ(kk,j,i) *( (P_t(kk+1,j+1,i+1)-P_t(kk,j+1,i+1)   -1/144* rho_water/Bw(kk,j,i)*(kk-3/2)*dz) - ...
                                      (Pcow(kk+1,j+1,i+1)-Pcow(kk,j+1,i+1) ));
          
            AA(2*m) = C_sww(kk,j,i)*( Sw(kk,j,i)-Sw_t(kk,j,i)   )+...
                     -C_pow(kk,j,i)*( P_t(kk+1,j+1,i+1)-P(kk,j,i) );

           QQ(2*m) = -Q_water(kk,j,i) -WCi(kk,j,i)/(dx*dy*dz) * (Landa_W(kk,j,i))*( P_t(kk+1,j+1,i+1)-P_bh(kk,j,i) ); ;

           y(2*m) = FF(2*m) + AA(2*m) + QQ(2*m)  ;



        end
    end
end
y=y';