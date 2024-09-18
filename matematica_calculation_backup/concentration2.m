clear all;
clc 
T0=2.35*10^-4; %eV
r0=4.7e-13;
alpha=1/100;
m=10^11; %eV
m_pl = 1.22*10^(19+9);
Ta = 0.2*10^6;
I=m*alpha^2/2; %eV
kappa=(43/345);
kappa_mod=0.045;

aB=1/(m*alpha); %eV^-1
M=10^12;
 
Trec=1*10^5;
TNY1=Trec*kappa^(-1/3)/sqrt(pi);
TRM=1.2;

i=1;
    T(i)=I;
    x(i)=1;
    x_w_th(i)=1;
    y(i)=1;
    x_real(i)=1;
    %x10(i)=1;
    z(i)=T(i)/T0-1;
    r3(i)=1
 while T(i)>(I/10)*kappa^(-1/3)
    i=i+1;
    T(i)=T(i-1)/2;
    x(i)=1;
    r3(i)=1
    x_w_th(i)=1;
    y(i)=1;
    x_real(i)=1;
    %x10(i)=1;
    z(i)=T(i)/T0-1;
 end

j=i;
T(i)=(I/10)*kappa^(-1/3);
x(i)=1;
y(i)=1;
r3(i)=1
x_w_th(i)=1;
x_real(i)=1;
%x10(1)=1;

while (T(i)>TNY1)
%while (T(i)>T0)
    i=i+1;
    T(i)=T(i-1)*0.9;
    z(i)=T(i)/T0-1;
    x(i)=(1+r0*(3.35/4)*10^14*((T(j))^0.1-(T(i))^0.1))^(-1);
    x_w_th(i)=x(i);
    y(i)=(1+4.6*10^6*r0*(T(j)^0.5*(log(I/(kappa^(1/3)*T(j)))+2)-T(i)^0.5*(log(I/(kappa^(1/3)*T(i)))+2)))^(-1);
    D = (2*pi^2*G_S(T(i))/45)**2 * (4*pi*alpha)^5 * (2*(m)^(1/2)*Ta^(9/2)/M)/(5.5*sqrt(G_E(TRM)/11)/m_pl);
    r3(i)=(1+D*r0*r0*(4/9)*(T(j)^(-9/2)-T(i)^(-9/2)))^(-1/2);
    %x10(i)=y(i);
    if T(i)/T0>1.5*10^6
    x_real(i)=y(i);
    k=i;
    r_rec_real=x_real(i)*r0;
    T_rec_real=T(i);
    else
        x_real(i)=(1+r_rec_real*(3.35/4)*10^14*((T(j))^0.1-(T(i))^0.1))^(-1);
    end
    
end

y_rec=y(i);
r3_t=r3(i)
x_rec=x(i);
x_w_th_rec=x_w_th(i);
r_rec_x=x_rec*r0;
r_rec_x_w_th=x_w_th_rec*r0;
r_rec_y=y_rec*r0;
T_rec=T(i);


while (T(i)>TRM)
    i=i+1;
    T(i)=T(i-1)*0.99;
    z(i)=T(i)/T0-1;
    x(i)=x_rec*(1+((23/20)/0.8)*10^17*r_rec_x*((1/T(i))^0.8-(1/T_rec)^0.8))^(-20/23);
    x_w_th(i)=x_rec*(1+(1/0.8)*10^17*r_rec_x*((1/T(i))^0.8-(1/T_rec)^0.8))^(-1);
    y(i)=y_rec*(1+4.71*10^8*r_rec_y*((log(sqrt(I*TNY1/T(i))))^2-(log(sqrt(I*TNY1/T_rec)))^2))^(-1);
    D = (2*pi^2*G_S(T(i))/45)**2 * (4*pi*alpha)^5 * (2*(m)^(1/2)*Ta^(9/2)/M)/(5.5*sqrt(G_E(TRM)/11)/m_pl);
    r3(i)=r3_t*(1+D*(r3_t*r0)^2*(4/9)*(T(i)^(-9/2)-T_rec^(-9/2)))^(-1/2);
    %x10(i)=y(i);
    if T(i)/T0>1.5*10^6
    x_real(i)=y(i);
    k=i;
    r_rec_real=x_real(i)*r0;
    T_rec_real=T(i);
    else
        x_real(i)=x_real(k)*(1+1.4375*10^17*r_rec_real*((1/T(i))^0.8-(1/T_rec_real)^0.8))^(-20/23);
    end
    
    
end

y_RM=y(i);
x_RM=x(i);
x_w_th_RM=x_w_th(i);
r_RM_x=x_RM*r0;
r_RM_x_w_th=x_w_th_RM*r0;
r_RM_y=y_RM*r0;
T_RM=T(i);
r3_t=r3(i)


while (T(i)>T0)
    i=i+1;
    T(i)=T(i-1)/1.1;
    z(i)=T(i)/T0-1;
    x(i)=x_rec*(1+((23/20)/0.8)*10^17*r_rec_x*((1/T_RM)^0.8-(1/T_rec)^0.8)+...
        (((23/20)/0.3)*10^17/sqrt(TRM))*r_rec_x*((1/T(i))^0.3-(1/T_RM)^0.3))^(-20/23);
    x_w_th(i)=x_rec*(1+(1/0.8)*10^17*r_rec_x*((1/T_RM)^0.8-(1/T_rec)^0.8)+...
        ((1/0.3)*10^17/sqrt(TRM))*r_rec_x*((1/T(i))^0.3-(1/T_RM)^0.3))^(-1);
    y(i)=y_RM*(1+1.7*10^9*y_RM*r0*(T_RM^0.5*(log(sqrt(I*TNY1/T_RM))+2)-T(i)^0.5*(log(sqrt(I*TNY1/T(i)))+2)))^(-1);
    D = (2*pi^2*G_S(T(i))/45)**2 * (4*pi*alpha)^5 * (2*(m)^(1/2)*Ta^(9/2)/M)/(3.08*10^(-14)/m_pl^(1/2));
    r3(i)=r3_t*(1+D*(r3_t*r0)^2*(2/5)*(T(i)^(-5)-T_RM^(-5)))^(-1/2);
    if T(i)/T0>1.5*10^6
    x_real(i)=y(i);
    k=i;
    r_rec_real=x_real(i)*r0;   
    T_rec_real=T(i);
    else
        x_real(i)=x_real(k)*(1+1.4375*10^17*r_rec_real*((1/T_RM)^0.8-(1/T_rec_real)^0.8)+...
        (1.4375*10^17/sqrt(TRM))*r_rec_real*((1/T(i))^0.3-(1/T_RM)^0.3))^(-20/23);
    end
    
    if T(i)/T0>11
    %x10(i)=y(i);
    k10=i;
    r10=x(i)*r0;
    r10_real=x_real(i)*r0;
    %r_rec_10=x10(i)*r0;   
    T_rec_10=T(i);
    w_th=x(i)/x_w_th(i);
    x(i);
    l=i;
    rl=x(i)*r0;
    %else
     %   x10(i)=x10(k10)*(1+(1.4375*10^17/sqrt(TRM))*r_rec_10*((1/T(i))^0.3-(1/T_rec_10)^0.3))^(-20/23);
    end
  
    
end
    %%
%loglog(T, x, T, y, T, x_real)
%axes('Position',[0.1,0.1,0.8,0.8])

%loglog(T, x, '--', 'Color', 'red', 'LineWidth' ,1.3)
%hold on;
#plot(T, y, '-.', 'Color', 'blue', 'LineWidth' ,1.3)
#hold on;
%loglog(T, x_real, 'Color', 'black', 'LineWidth' ,1.3)
%hold on;
%loglog(T, r3, 'Color', 'green', 'LineWidth' ,1.3)
%hold on;
%loglog(T, x_w_th);

#ylim ([0.96 1])
#xlim ([10^-5 10^6])

#set(gca,'fontsize',11,'FontName','times');

#xlabel('$T$, eV','fontsize',18,'interpreter','latex');
#ylabel('$\frac{r}{r_0}$','fontsize',35,'interpreter','latex');
##
##ax = gca;
####
###ax.XTickLabel={10^5, 10^4, 10^3, 10^2, 10^1, 1};
####
###text(9*10^-1,1.5*10^-4,'Q-C','interpreter','latex','fontsize',18, 'Rotation', -15);
##text(50*10^-1,2,'Quantum','interpreter','latex','fontsize',18, 'Color', 'blue');
###text(5*10^3,0.4,'Three-body','interpreter','latex','fontsize',18, 'Color', 'green', 'Rotation', -20);
####
####
####%%
##loglog(T, y, '-.', 'Color', 'blue', 'LineWidth' ,5)
##ylim ([0.95 1])
##xlim ([T0*11 I])
##
##set(gca,'fontsize',20,'FontName','times');
##
##ax = gca;
##
##
##xlabel('$T$, eV','fontsize',28,'interpreter','latex');
##ylabel('$\frac{r}{r_0}$','fontsize',40,'interpreter','latex');
##%%
##for k=1:1:length(T)
##    
##    if T>0.5*10^5
##        E0(k)=1.76*10^(-4)*T(k)/T0; 
##    else
##        E0(k)=8.3*10^(-13)*(T(k)/T0)^2; %eV
##    end;
##    
##Lat(k)=4.6*10^3*M^(1/3)*T0/T(k);
##Lat_real(k)=Lat(k)*(x_real(k))^(-1/3);
##
##RB(k)=(2*sqrt(2)/15)^(2/5)*(alpha/m)*(m/E0(k))^(2/5); %eV^-1
##R10(k)=(1/0.1)^(2/5)*RB(k); %eV^-1
##Lloss(k)=R10(k)-RB(k); %eV^-1
##
##Lloss_Lat(k)=Lloss(k)/Lat(k);
##
##pmax(k)=alpha*((sqrt(2)*pi/4)*(m/E0(k)))^0.2/sqrt(m*E0(k));
##pmaxRB(k)=pmax(k)/RB(k);
##
##pm_Lat(k)=pmax(k)/Lat(k);
##
##r1(k)=(0.5*((alpha/E0(k))^1.5+((alpha/E0(k))^3+(2*sqrt(2)/3)*(alpha^2/(m*E0(k)))^1.5)^0.5))^(2/3);
##r1Lat(k)=r1(k)/Lat(k);
##
##
##r1Lat_real(k)=r1(k)/Lat_real(k);
##
##pm_Lat_real(k)=pmax(k)/Lat_real(k);
##Lloss_Lat_real(k)=Lloss(k)/Lat_real(k);
##
##r2(k)=alpha/E0(k);
##r12(k)=r1(k)/r2(k);
##end
##%%
##loglog (T, Lloss_Lat, '--', 'Color', 'black','linewidth',1.3)
##hold on;
##loglog (T, Lloss_Lat_real, 'Color', 'black','linewidth',1.3);
##
##xlim([T0 I/10])
##ylim([1e-15 0.5e6])
##set(gca,'fontsize',14,'FontName','times');
##xlabel('$T$, eV','fontsize',18,'interpreter','latex');
##
##loglog(T, pm_Lat, '--', 'Color', 'red','linewidth',1.3);
##loglog (T, pm_Lat_real, 'Color', 'red','linewidth',1.3);
##
##loglog (T, r1Lat, '--', 'Color', 'blue','linewidth',1.3);
##loglog(T, r1Lat_real, 'Color', 'blue','linewidth',1.3);
##
##
##ax = gca;
##
##
##text(6*10^-3,0.3*10^-9,'$\frac{L_{\rm{loss}}}{L_{\rm{sp}}}$','interpreter','latex','fontsize',25, 'Rotation', 0);
##text(6*10^-3,5*10^3,'$\frac{r_{*}}{L_{\rm{sp}}}$','interpreter','latex','fontsize',25, 'Color', 'blue');
##text(6*10^-3,2*10^-3,'$\frac{\rho_{\rm{max}}}{L_{\rm{sp}}}$','interpreter','latex','fontsize',25, 'Color', 'red', 'Rotation', 0);
##
##hold off;
##
##%%
loglog(T, y, '--', 'Color', 'red', 'LineWidth' ,1.3)
hold on;
loglog(T, x, '-', 'Color', 'black', 'LineWidth' ,1.3)
hold on;
loglog(T, r3, '-', 'Color', 'blue', 'LineWidth' ,1.3)

%%
##loglog (T, Lloss_Lat, '--', 'Color', 'black','linewidth',1.3)
##hold on;
##loglog (T, Lloss_Lat_real, 'Color', 'black','linewidth',1.3);
##

ylim ([10^(-7) 1])
xlim ([T0 I])

set(gca,'fontsize',18,'FontName','times');

xlabel('$T$, eV','fontsize',18,'interpreter','latex');
ylabel('$\frac{r}{r_0}$','fontsize',35,'interpreter','latex');
##
ax = gca;
##
##
##
##text(2*10^-1,1.5,'Quantum','interpreter','latex','fontsize',18, 'Color', 'red');
##text(5*10^5,5*10^-2,'Classical','interpreter','latex','fontsize',18, 'Color', 'black', 'Rotation', -30);
##hold off;
## %%
## plot (T, r12);
##%%
##loglog (T, Lloss_Lat, '--', 'Color', 'black','linewidth',1.3)
##hold on;
##loglog (T, Lloss_Lat_real, 'Color', 'black','linewidth',1.3);
##
##xlim([T0 I/10])
##ylim([1e-14 1e-8])
##set(gca,'fontsize',14,'FontName','times');
##xlabel('$T$, eV','fontsize',18,'interpreter','latex');
##ylabel('$\frac{L_{\rm loss}}{L_{\rm sp}}$','fontsize',30,'interpreter','latex');
##
##
##
##ax = gca;
##
##
##
##
## 