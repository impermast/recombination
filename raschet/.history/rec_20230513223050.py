
from math import gamma as Gamma
import numpy as np
from numpy import sqrt as sqrt
import matplotlib.pyplot as plt
from matplotlib import colors,ticker
from scipy import integrate 
import seaborn as sns
import pandas as pd


#LOGIC -- 1 download data, else generate new
logic = 1
err=10**(-2)
#consts
#E ~ eV
ma = 100*10**9 #ev
mb = 1*10**12 #ev
alpha = 1/10 #none

beta = 14/5 
sigma = 8*np.pi*(alpha/ma)**2/3
r0 = 4.7*10**(-13) #r T_rec
T0 = 2.35*10**(-4) #our time
m_pl = 1.22*10**(19+9) #plank mass ev
h_md = 3.08*10**(-14)  #from artile
t_rm = 1.2             #time of rd md transition
t_mod = 4.4*10**(17+16)/(6.582) #
I = ma*alpha*alpha/2 #ev Ionisation potential
mu = ma*mb/(ma+mb)  #reduced weight
s0 = (4*np.pi)**(2/5) *np.pi* (alpha/mu)**2 #coef in kramers
kappa = 43/345  
B = alpha*mu/2 
Ta = 0.2*10**6*(ma/(100*10**9))**(3/2)*(1/(100*alpha))
Tny = Ta*kappa**(-1/3)/sqrt(np.pi)

def g_s(T):
    return 2+2*8*(T>200*10**6)+ 7/8*(2*2*(T>0.5/3*10**6))+ 2*2*(T>100/3*10**6) + 2*2*(T>2/3*10**9) + 2*1*3*(T>=10**6)+2*1*3*4/11*(T<10**6)+ 2*2*3*2*(T>200*10**6) + 2*2*3*(T>200*10**6) + 2*2*3*(T>1.8/3*10**9) + 2*2*3*(T>4.5/3*10**9)+ 7/8*2*2*3*(T>171/3*10**9)
def g_e(T):
    return 2+2*8*(T>200*10**6)+ 7/8*(2*2*(T>0.5/3*10**6))+ 2*2*(T>100/3*10**6) + 2*2*(T>2/3*10**9) + 2*1*3*(T>=10**6)+ 2*1*3*(4/11)**(4/3)*(T<10**6)+ 2*2*3*2*(T>200*10**6) + 2*2*3*(T>200*10**6) + 2*2*3*(T>1.8/3*10**9) + 2*2*3*(T>4.5/3*10**9)+ 7/8*2*2*3*(T>171/3*10**9)


def H(t):
    # Hubble const on RD and MD from temperature
    if t < t_rm:
        #MD
        return h_md/m_pl**(1/2)
    else:
        return 5.5*np.sqrt(g_e(t_rm)/11)/m_pl
        

def lim_3b(x):
    """
    Function that calculs limit on 3body recombination density

    x(list):    Temperature of O-fotons
    
    return list of density of unrecombined particles
    """
    r_pit = (Ta**2/(x**4)) / ((4*np.pi*alpha)**3 *r0*(2*np.pi**2*g_s(x)/45)*sqrt(ma/mb))
    return r_pit
        
def graph(T,r1,r2,rlim):
    # graphplotter used to plot r(T) graphs
    fig, ax = plt.subplots(1,1)
    ax.plot(T,r1,'b',label ="Трехчастичная рекомбинация")
    ax.plot(T,rlim, "--r",label ="Limit")
    ax.plot(T,r2,'g', label ="Радиационная рекомбинация")
    ax.set_xlabel('T, eV', fontsize=14)
    ax.set_ylabel(r'$r/r_0$', fontsize=14, fontweight='bold')
    ax.tick_params(labelsize=12)
    ax.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.legend(fontsize=12)
    ax.invert_xaxis()
    plt.semilogx()
    plt.semilogy()
    plt.show()

def heatplot(df):
    # heatplot of maximal denity for different model parameters
    sns.set_theme(style="dark")
    sns.set_style("ticks")
    fig, ax = plt.subplots(figsize=(9, 6))
    sns.heatmap(df, annot=True, linewidths=.5, ax=ax, cmap='coolwarm', cbar_kws={'label': 'Значения'})
    ax.set(xlabel=r'$m_a$, eV', ylabel=r'$\alpha_y$', title='Относительная плотность прорекомбинировавших частиц для различных параметров модели')
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f'$10^{{{(x-0.5):.0f}}}$'))
    ylabels = ["{:0.1f}".format(float(item.get_text())) 
                for item in ax.get_yticklabels()] 
    xlabels = ["$10^{%0.1f}$"%(float(item.get_text())) for item in ax.get_xticklabels()]
    ax.set_yticklabels(ylabels)
    ax.set_xticklabels(xlabels)
    plt.xticks(rotation=0)

    # Add contour lines
    # x = np.arange(df.shape[1])
    # y = np.arange(df.shape[0])
    # X, Y = np.meshgrid(x, y)
    # levels = np.linspace(df.values.min(), df.values.max(), 10)
    # CS = ax.contour(X + 0.5, Y + 0.5, df.values, levels, colors='k')
    # ax.clabel(CS, inline=True, fontsize=10, fmt='%1.1f')
    # plt.show()
    
def kriteria(alp_range,m_range):
    """
    Creation of dataframe of maximal density of recombined particles from different alpha and ma (model parameters)

    dens(list):  2D massive of density data
    alp_range(list):    1D list of parameters alpha
    m_range(list):      1D list of electron mass parameter
    krit(list):         1D list of density, that used to find point where limit and 3body recombination density are equal

    create csv file of dataframe
    """
    dens = [0]*len(alp_range)
    dens1 = []
    k=0
    for alpha in alp_range:
        for ma in m_range:
            [T,rk,r3,rc,rlim] = calc(alpha,ma,mb)
            krit=[i for i, j in zip(r3, rlim) if abs(i-j)<err and i<0.1]
            try: dens1.append(float(-np.log(krit[-1])))
            except IndexError: dens1.append(0)
        dens[k]=dens1
        dens1=[]
        print("At ",k+1,"  itteration")
        k+=1
    df = pd.DataFrame(dens, index=alp_range,columns=m_range)
    df.to_csv('data.csv')

def calc(alpha,ma,mb):
    """
    Calculates the temperature-dependent rates for various processes.

    Args:
        I (float): Initial temperature of photons.
        T0 (float): Present-day temperature of photons.
        kappa (float): Ratio of photon number density to baryon number density.
        lim_3b (function): Function that calculates the three-body rate limit.
        t_rm (float): Time of radiation-matter equality.
        alpha (float): Fine-structure constant.
        ma (float): Mass of electron.
        mb (float): Mass of baryon.
        beta (float): Power law index.
        g_s (function): Function that calculates the entropy degrees of freedom.
        s0 (float): Present-day entropy density.
        m_pl (float): Planck mass.
        mu (float): Ratio of proton mass to electron mass.

    Returns:
        list: A list containing the temperature-dependent rates for various processes.
    """
    T = [] # Temperature of photons
    z = [] #redshift
    rk = [] # r for kramerz

    #after decouple of dphoton
    i = 0
    T = [(I)] 
    z = [T[i]/T0 -1]
    rk = [1] #initial r=1 Kramers case
    r3b = [1] #3body case
    rcl = [1] #classical
    rlim =[lim_3b(T[i])] #pitevski lim

    #until O-photons decouple
    while T[i]>(I/10)*kappa**(-1/3):
        i = i+1
        T = T + [T[i-1]/2]
        rk = rk + [1]
        r3b = r3b+[1]
        rcl = rcl+[1]
        rlim = rlim + [lim_3b(T[i])]
        z.append(T[i]/T0 -1)

    # until RD MD transition
    j = i
    while T[i]>t_rm:
        i=i+1
        T = T + [T[i-1]*0.9]
        z.append(T[i]/T0 -1)
        rk = rk + [(1+4.71*10**8*r0*((np.log(np.sqrt(I*Ta/T[i])))**2-(np.log(np.sqrt(I*Ta/T[j])))**2))**(-1)]
        D = (2*np.pi**2*g_s(T[i])/45)**2 * (4*np.pi*alpha)**5 * (2*np.sqrt(ma)*Ta**(9/2)/mb) / H(T[i])   
        r3b = r3b + [(1+ D*2/5*r0*r0*(T[i]**(-5)-I**(-5)))**(-1/2)]
        
        Dcl = Gamma(2-beta/2)*g_s(T[i])*s0*m_pl*(mu*Ta)**((beta-1)/2)/ ((beta-2)*2**((beta-3)/2)*np.sqrt(45*g_s(T[j])))
        rcl = rcl + [(1+Dcl*r0*(T[i]**(2-beta)-T[j]**(2-beta)))**(-1)]
        rlim = rlim + [lim_3b(T[i])]
        
    #until now
    j=i
    rk_rm=rk[j]
    r3b_rm = r3b[j]
    while T[i]>T0:
        i=i+1
        T = T + [T[i-1]*0.99]
        z.append(T[i]/T0 -1)
        rk = rk + [rk_rm*(1+1.7*10**9*rk_rm*r0*(T[j]**0.5*(np.log(np.sqrt(I*Ta/T[j]))+2)-T[i]**0.5*(np.log(np.sqrt(I*Ta/T[i]))+2)))**(-1)]
        D = (2*np.pi**2*g_s(T[i])/45)**2 * (4*np.pi*alpha)**5 * (2*np.sqrt(ma)*Ta**(9/2)/mb) / H(T[i])
        r3b = r3b +[r3b_rm*(1+D*r3b_rm*r3b_rm*r0*r0* 4/9 * (T[i]**(-9/2)-T[j]**(-9/2)))**(-1/2)]
        Dcl = Gamma(2-beta/2)*g_s(T[i])*s0*m_pl*(mu*Ta)**((beta-1)/2)/ ((beta-2.5)*2**((beta-3)/2)*np.sqrt(45*g_s(T[j])))/np.sqrt(t_rm)
        rcl = rcl + [rcl[j]*(1+Dcl*rcl[j]*r0*(T[i]**(2.5-beta)-T[j]**(2.5-beta)))**(-1)]
        rlim = rlim + [lim_3b(T[i])]       

    return [T,rk,r3b,rcl,rlim]
    
def main():
    alp_range = np.logspace(-1,3,10)
    m_range = np.logspace(16,21,10)
    if logic == 1:
        try:
            df = pd.read_csv('data.csv')
        except FileNotFoundError:
            print(u'Saved data not found')
            df = kriteria(alp_range,m_range)
    else:
        kriteria(alp_range,m_range)

    df = pd.read_csv('data.csv', index_col=0)
    heatplot(df)
    print(df)
    
    [T,rk,r3,rcl,rlim] = calc(alpha,ma,mb)
    graph(T,r3,rk,rlim)
main()