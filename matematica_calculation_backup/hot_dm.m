clear all




M_norm1=100e9; %eV
M_norm2=50e9; %eV

M_DM_GeV=[1e-4, 1e-3, 2e-3, 5e-3, 7e-3, 1e-2, 2e-2, 5e-2, 7e-2, 1e-1, 2e-1, 5e-1, 7e-1, 1, 2, 5, 7, 10, 20, 50, 70, 100,...
    200, 500, 700, 1000, 2000, 5000, 7000, 1e4, 2e4, 5e4, 7e4, 1e5, 2e5, 5e5, 7e5, 1e6, 1e7]; 
%GeV

alpha_q=[1e-10, 2e-10, 5e-10, 1e-9, 2e-9, 5e-9, 1e-8, 2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, ...
    1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4];

M_DM=M_DM_GeV*10^9;

z=10;
E0_m_al=2.32*10^6*(z+1)^2; %eV

T0=2.35e-4; %eV
z10=10;
z1000=1000;
T10=T0*(z10+1); %eV
T1000=T0*(z1000+1); %eV

rho_mod=5.5*10^3; %eV/cm^3
s_mod=2890; %1/cm^3
T_RM=1.2; %eV
hc=1.97e-5; %eV*fm

ksi=0.448;
mPl=1.22e28; %eV
hMD=3.08e-14;
gs=43/22; %z=10,1000
EV=1.6e-19; %э¬ в дж
c=3e8; %м/с

SM=1.2*(1e13)^2*EV/(1e-3*c^2*hc^2); %с переводом из см^2/г в э¬^-3


v_b=1e6/c; %1000 км/с
sin_theta=pi/4;

for i=1:1:length(M_DM)
one(i)=1;    
    
alpha_max(i)=2*(M_DM(i)/M_norm1)^1.5;

alpha_cq(i)=((15/sqrt(2))*E0_m_al*(M_DM(i))^(-2.5))^(1/4);

alpha_cq0(i)=((15/sqrt(2))*E0_m_al/(z+1)^2*(M_DM(i))^(-2.5))^(1/4);

alpha10(i)=1/30*(0.3*T10^0.3*s_mod*M_DM(i)*sqrt(T_RM)/(0.25*rho_mod*1.6e18)*(M_DM(i)/M_norm2)^(-0.25))^(10/11);

alpha1000(i)=1/30*(0.3*T1000^0.3*s_mod*M_DM(i)*sqrt(T_RM)/(0.25*rho_mod*1.6e18)*(M_DM(i)/M_norm2)^(-0.25))^(10/11);

r0(i)=0.25*rho_mod/(2*M_DM(i)*s_mod);

I(i)=1/2*alpha_cq(i)^2*(M_DM(i)/2);
func_I = @(alpha) 1/2*alpha^2*(M_DM(i)/2);

Tdec(i)=0.2e6*(M_DM(i)/M_norm1)^1.5*0.01/alpha_cq(i);
func_Tdec = @(alpha) 0.2e6*(M_DM(i)/M_norm1)^1.5*0.01/alpha;

term10(i)=ksi*64*sqrt(pi)/(3*sqrt(3))*(alpha_cq(i)/(M_DM(i)/2))^2*sqrt(I(i)*Tdec(i))*sqrt(mPl)/hMD*2*pi^2/45*gs*4*...
   (sqrt(T_RM)*(log(sqrt(I(i)*Tdec(i))/T_RM)+2) - sqrt(T10)*(log(sqrt(I(i)*Tdec(i))/T10)+2));

func_term10 = @(alpha) ksi*64*sqrt(pi)/(3*sqrt(3))*(alpha/(M_DM(i)/2))^2*sqrt(func_I(alpha)*func_Tdec(alpha))...
*sqrt(mPl)/hMD*2*pi^2/45*gs*4*(sqrt(T_RM)*...
(log(sqrt(func_I(alpha)*func_Tdec(alpha))/T_RM)+2) - sqrt(T10)*(log(sqrt(func_I(alpha)*func_Tdec(alpha))/T10)+2));

r10(i)=r0(i)/(1+r0(i)*term10(i));
r10_r0(i)=1/(1+r0(i)*term10(i));

func_r10_r0 = @(alpha) 1/(1+r0(i)*func_term10(alpha));

func = @(alpha) func_r10_r0(alpha) - 0.99;

func_svm = @ (alpha) ksi*64*pi/sqrt(27*pi)*(alpha/(M_DM(i)/2))^2*sqrt(func_I(alpha)/((T0*(1+0.3))^2/func_Tdec(alpha)))*...
    log(func_I(alpha)/((T0*(1+0.3))^2/func_Tdec(alpha)))/M_DM(i);

%alpha_bullet(i)= fsolve(func_svm,1e-20);

if i<15
    alpha_true(i) = fsolve(func,1e-4);
    break_line(i)=1e-1;
else
    alpha_true(i) = alpha_cq(i);
    break_line(i)=1e11;
end

LL=63;
alpha_b(i)=(1/SM*ksi*64*pi/sqrt(27*pi)*2/M_DM(i)^3*sqrt(0.2e6*(M_DM(i)/M_norm1)^1.5*0.01*M_DM(i))/(T0*(0.3+1))*LL)^...
    (-2/5);

x_fr=1/30;

T_fr(i)=M_DM(i)*x_fr;
hRD_fr(i)=5.5*sqrt(G_E(T_fr(i))/11);
H_fr(i)=hRD_fr(i)*T_fr(i)^2/mPl;


for j=1:1:length(alpha_q)
    
    sv_fr(i,j)=pi/2*alpha_q(j)^2/M_DM(i)^2;
    n_fr(i,j)=H_fr(i)/sv_fr(i,j);
    rho_fr(i,j)=2*M_DM(i)*n_fr(i,j)/(hc)^3;
    Omega_fr(i,j)=rho_fr(i,j)*(T0/T_fr(i))^3/rho_mod;
    
    sv_fr_GeV(i,j)=pi/2*alpha_q(j)^2/M_DM_GeV(i)^2;
    Omega_fr(i,j)=3.99e-11/(sv_fr_GeV(i,j)*x_fr);
    
    BULLET(i,j)=func_svm(alpha_q(j));
    
    I_L(i,j)=1/2*alpha_q(j)^2*M_DM(i)/2;
    Tdec_L(i,j)=0.2e6*(M_DM(i)/M_norm1)^1.5*0.01/alpha_q(j);
    L(i,j)=log(I_L(i,j)/(T0^2*1.3^2/Tdec_L(i,j)));
    
end

alpha_q_bul_min(i)=SM*M_DM(i)^3*v_b^4/(4*pi*(1/sin_theta^2-1));

%func_sv_fr_GeV = @(alpha) pi/2*alpha^2/M_DM_GeV(i)^2;
%func_Omega = @(alpha) 3.99e-11/(sv_fr_GeV(alpha)*x_fr) - 0.25;

%alpha_q_true(i) = fsolve(func,3e-7);

end

alpha_q_true=[1e-10, 5.521e-8, 1.1043e-7, 2.7605e-7, 3.865e-7, 5.521e-7, 1.1043e-6, 2.7605e-6, 3.865e-6, 5.521e-6, 1.1043e-5,...
    2.7605e-5, 3.865e-5, 5.521e-5, 1.1043e-4, 2.7607e-4, 3.865e-4];




for j=1:1:length(alpha_q_true)
    
    M_DM_GeV_q_true(j)=M_DM_GeV(j);
    
    sv_fr_GeV_check(j)=pi/2*alpha_q_true(j)^2/M_DM_GeV(j)^2;
    Omega_fr_check(j)=3.99e-11/(sv_fr_GeV_check(j)*x_fr);
    
                    %от параметра хаббла                                %сечение аннигил€ции
    x_fr_1(j)=(log(3*mPl/(4*pi^3)*sqrt(5/(2*G_E(T_fr(j))))*G_S(T_fr(j))*pi/2*alpha_q_true(j)^2/M_DM(j)^2*...
               M_DM(j)*x_fr^(-1/2)))^(-1);
           
   
    T_fr_1(j)=M_DM(j)*x_fr_1(j);
    
    hRD_fr_1(j)=5.5*sqrt(G_E(T_fr_1(j))/11);
    H_fr_1(j)=hRD_fr_1(j)*T_fr_1(j)^2/mPl;
    
end

alpha_q_true_1=[1e-10, 4.1656e-08, 8.499e-08, 2.1788e-7, 3.1e-7, 4.47e-7, 9.095e-7, 2.3325e-6, 3.291e-6, 4.74e-6, 9.628e-6,...
    2.4548e-5, 3.461e-5, 4.981e-5, 1.013e-4, 2.578e-4, 3.682e-4];

%for i=1:1:length(x_fr_1)
 %   for k=1:1:length(alpha_q)
    
  %      sv_fr_GeV(i,k)=pi/2*alpha_q(k)^2/M_DM_GeV(i)^2;
   %     Omega_fr_1(i,k)=3.99e-11/(sv_fr_GeV(i,k)*x_fr_1(i));
    
    %end
%end

for k=1:1:length(alpha_q_true_1)
    
    %Tdecl1(k)=0.2e6*(M_DM(k)/M_norm1)^1.5*0.01/alpha_q_true_1(k);
    %TDM1(k)=T0^2*(10+1)^2/Tdecl1(k);
    %E0(k)=3/2*TDM1(k);
    %ACTION(k)=2*alpha_q_true_1(k)*sqrt(M_DM(k)/E0(k))*(sqrt(3)-atan(sqrt(3)));
    
    
    %прицельный параметр соответствует рассе€нию на пи/2
    
    ACTION(k)=2*alpha_q_true_1(k)/v_b*(sqrt(3)-atan(sqrt(3)));
    s_m(k)=alpha_q_true_1(k)^2/M_DM(k)^3*4*pi/v_b^4*(1/sin_theta^2-1);
    

    
end

%%


subplot(2,1,1);

loglog (M_DM_GeV, alpha_q_bul_min, 'Color', 'black','linewidth',1.3)

h1 = line([M_DM_GeV M_DM_GeV(length(M_DM_GeV)) M_DM_GeV(1)  M_DM_GeV(1) ], ...
   [alpha_max 1e41 1e41 alpha_max(1)],'Color','r','linewidth',2);
hatch(h1,90,'r','-',6,1);

ylim([10^10 1e40])
xlim([1e-3 1e6])
set(gca,'fontsize',14,'FontName','times');
set(gca,'XTick',[])

g1 = subplot(2,1,1) 

p1 = get(g1,'position');
p1(4) = p1(4)*0.5; 

set(g1, 'position', p1)

%FigHandle = figure;
 % set(FigHandle, 'Position', [100, 100, 1049, 895]);


subplot(2,1,2);

loglog (M_DM_GeV_q_true, alpha_q_true_1, 'Color', 'g','linewidth',1.3)

hold on

fill([M_DM_GeV fliplr(M_DM_GeV)], [alpha_cq fliplr(one)], 'white')

loglog (M_DM_GeV, alpha_max, 'Color', 'r')
%h=text(13,10^-5,'COLD DM');
%s = h.FontSize;
%h.FontSize = 18;

h0 = line([M_DM_GeV M_DM_GeV(length(M_DM_GeV)) M_DM_GeV(1)  M_DM_GeV(1) ], ...
    [alpha_max 1e30 1e30 alpha_max(1)],'Color','r','linewidth',2);
hatch(h0,90,'r','-',6,1);

ylim([10^-10 1])
xlim([1e-3 1e6])

%g=text(1e-1,10,'HOT DM', 'BackgroundColor','w');
%s = g.FontSize;
%g.FontSize = 18;

set(gca,'fontsize',14,'FontName','times');
xlabel('$M_{\rm DM}$, GeV','fontsize',18,'interpreter','latex');
ylabel('$\alpha_y$','fontsize',18,'interpreter','latex')


loglog (M_DM_GeV, alpha_true, 'Color', 'm')



h1 = line([M_DM_GeV M_DM_GeV(length(M_DM_GeV)) M_DM_GeV(1)  M_DM_GeV(1) ], ...
    [alpha_true alpha_cq(length(M_DM_GeV)) alpha_cq(1) alpha_true(1)],'Color','m','linewidth',2);
hatch(h1,90,'m','-',6,1);



%h = area(alpha_cq,'LineStyle',':');

loglog (M_DM_GeV, alpha_cq, 'Color', 'b','linewidth',1.3)

loglog (M_DM_GeV, alpha_cq0, 'linewidth',1.3)



g0 = subplot(2,1,2) 

p0 = get(g0,'position');
p0(2) = p0(2)*2; 
set(g0, 'position', p0)



%breakyaxis([1 10]);

%BreakLoglog(M_DM_GeV, break_line ,1,10^10,'RPatch') 


%loglog (M_DM_GeV, alpha10, 'Color', 'g')

%loglog (M_DM_GeV, alpha1000, 'Color', 'black')

%%
%a=20*rand(21,1)+10; 
%figure;hold on; 
%plot(a); 
%breakyaxis([14 21]);

%Y = [1, 5, 3;
 %   3, 2, 7];
%h = area(Y,'LineStyle',':');
%h(1).FaceColor = [0 0 0];
