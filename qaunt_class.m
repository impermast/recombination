clear all;
z=10;
E0_m_al=2.32*10^6*(z+1)^2; %eV
for i=1:1:1000
    M_DM(i)=i; %GeV
    %alpha1(i)=(2*E0/m(i))^0.2
    alpha(i)=((15/(2*sqrt(2)))*E0_m_al*(M_DM(i)*10^9)^(-2.5))^(1/4);
   % m(i)=i;
end

loglog (M_DM, alpha, 'Color', 'b')
h=text(13,10^-3,'Classical approximation');
s = h.FontSize;
h.FontSize = 18;


%hold on;
%for j=1:10:1000
 %   ms(j)=j*10^9; %eV
    %alpha1(i)=(2*E0/m(i))^0.2
  %  alphas(j)=((15/(2*sqrt(2)))*E0/m(j))^0.2;
%end

%stem(ms, alphas, 'Color', 'b', 'Marker', 'none');
%hold on;
h0 = line([M_DM M_DM(length(M_DM))+1 M_DM(1)-0.1 M_DM(1)-0.1], [alpha 10^-10 10^-10 alpha(1)],'Color','b','linewidth',2);
hatch(h0,90,'b','-',6,1);

ylim([10^-6 10^-2])
xlim([M_DM(1) M_DM(length(M_DM))-1])

g=text(1.5,8*10^-6,'Quantum approximation', 'BackgroundColor','w');
s = g.FontSize;
g.FontSize = 18;

set(gca,'fontsize',18,'FontName','times');
xlabel('$\mu$, GeV','fontsize',18,'interpreter','latex');
ylabel('$\alpha_y$','fontsize',18,'interpreter','latex')


