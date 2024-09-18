clear all

M_norm=50e9; %eV

M_DM_GeV=[7e-3, 1e-2, 2e-2, 5e-2, 7e-2, 1e-1, 2e-1, 5e-1, 7e-1, 1, 2, 5, 7, 10, 20, 50, 70, 100, 200, 500, 700, ...
    1000, 2000, 5000, 7000, 1e4]; 
%GeV

M_DM=M_DM_GeV*10^9; %eV

T0=2.35e-4; %eV
z=10;
T10=T0*(z+1); %eV

rho_mod=5.5*10^3; %eV/cm^3
s_mod=2890; %1/cm^3
T_RM=1.2; %eV

for i=1:1:length(M_DM)

    alpha(i)=1/30*(2*0.3*T10^0.3*s_mod*(z+1)^3*M_DM(i)*sqrt(T_RM)/(0.25*rho_mod*1.6e18)*(M_DM(i)/M_norm)^(-0.25))^(10/11);
end

loglog (M_DM_GeV, alpha, 'Color', 'r')
