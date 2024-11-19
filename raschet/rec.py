#rec.py
from math import gamma as Gamma
from scipy.ndimage import gaussian_filter

import numpy as np
from numpy import sqrt as sqrt
import matplotlib.pyplot as plt
from matplotlib import colors,ticker
import pandas as pd

from tqdm import tqdm


#consts
#E ~ eV
beta = 14/5 
r0 = 4.7*10**(-13) #r T_rec
T0 = 2.35*10**(-4) #our time
m_pl = 1.22*10**(19+9) #plank mass ev
h_md = 3.08*10**(-14)  #from artile
t_rm = 1.2             #time of rd md transition
t_mod = 4.4*10**(17+16)/(6.582) #
kappa = 43/345  


def smooth_step(T, T1):
    """ Smooth transition around T1 with parameter k using arctan,
        returns 0 or 1 for values far from T1 to speed up computation """
    k=1/T1**(0.8)
    if T<T1*0.5:
        return 0 
    elif T>T1*2:
        return 1  
    else:
        return (1 + 2*np.arctan(k * (T - T1)) / np.pi )
def g_s(T):
    return 2+2*8*smooth_step(T,200*10**6)+ 7/8*(2*2*smooth_step(T,0.5/3*10**6))+ 2*2*smooth_step(T,100/3*10**6) + 2*2*smooth_step(T,2/3*10**9) + 2*1*3*smooth_step(T,10**6)+2*1*3*4/11*smooth_step(T,10**6)+ 2*2*3*2*smooth_step(T,200*10**6) + 2*2*3*smooth_step(T,200*10**6) + 2*2*3*smooth_step(T,1.8/3*10**9) + 2*2*3*smooth_step(T,4.5/3*10**9)+ 7/8*2*2*3*smooth_step(T,171/3*10**9)
def g_e(T):
    """ Smooth version of g_e function """
    return (
        2 +
        2 * 8 * smooth_step(T, 200 * 10**6) +
        7 / 8 * (2 * 2 * smooth_step(T, 0.5 / 3 * 10**6)) +
        2 * 2 * smooth_step(T, 100 / 3 * 10**6) +
        2 * 2 * smooth_step(T, 2 / 3 * 10**9) +
        2 * 1 * 3 * smooth_step(T, 10**6) +
        2 * 1 * 3 * 4 / 11 * smooth_step(T, 10**6) +
        2 * 2 * 3 * 2 * smooth_step(T, 200 * 10**6) +
        2 * 2 * 3 * smooth_step(T, 200 * 10**6) +
        2 * 2 * 3 * smooth_step(T, 1.8 / 3 * 10**9) +
        2 * 2 * 3 * smooth_step(T, 4.5 / 3 * 10**9) +
        7 / 8 * 2 * 2 * 3 * smooth_step(T, 171 / 3 * 10**9)
    )

def H(t):
    # Hubble const on RD and MD from temperature
    if t < t_rm:
        #MD
        return h_md/m_pl**(1/2)
    else:
        return 5.5*np.sqrt(g_e(t_rm)/11)/m_pl
        

def lim_3b(x,Ta):
    """
    Function that calculs limit on 3body recombination density

    x(list):    Temperature of O-fotons
    
    return list of density of unrecombined particles
    """
    r_pit = (Ta**2/(x**4)) / ((4*np.pi*alpha)**3 *r0*(2*np.pi**2* g_s(x) / 45)*sqrt(ma/mb))
    return r_pit
        
def exp_constrait(alpha, m):
    
    """
    Observational limit on crosssection: 0.1 < sigma/m < 1 g/sm
    m0 = 100GeV in eV, alpha0 = 0.01, v=200 km/s

    This function returns sigma/m function for retherford cross
    """
    return 1.43* (alpha/0.01)**2 / (m/(100*10**9))**3
    
def kriteria(alp_range,m_range,mb_range):
    """
    Creation of dataframe of maximal density of recombined particles from different alpha and ma (model parameters)

    dens(list):  2D massive of density data
    alp_range(list):    1D list of parameters alpha
    m_range(list):      1D list of electron mass parameter
    krit(list):         1D list of density, that used to find point where limit and 3body recombination density are equal

    create csv file of dataframe
    """
    points = len(alp_range)
    name = "data_nolog_"+str(points)+".scv"
    dens = [0]*len(alp_range)
    for k, a in tqdm(enumerate(alp_range), ncols=100,total=len(alp_range), desc='Calculating'):
        dens1 = []  # Переменная для хранения текущего результата
        for m in m_range:
            [T, rk, r3, rc, rlim] = calc(a, m, m * mb_range)
            krit = r3[np.argmin(np.abs(np.array(r3) - np.array(rlim)))]
            if krit > 0.98:
                dens1.append(0)
            else:
                dens1.append((1 - krit))
        dens[k] = dens1  # Записываем результат

    df = pd.DataFrame(dens, index=alp_range,columns=m_range)
    df.to_csv(name)

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

    I = ma*alpha*alpha/2 #ev Ionisation potential
    mu = ma*mb/(ma+mb)  #reduced weight
    s0 = (4*np.pi)**(2/5) *np.pi* (alpha/mu)**2 #coef in kramers
    B = alpha*mu/2 
    Ta = 0.2*10**6*(ma/(100*10**9))**(3/2)*(1/(100*alpha))
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
    rlim =[lim_3b(T[i],Ta)] #pitevski lim

    #until O-photons decouple
    while T[i]>(I/10)*kappa**(-1/3):
        i = i+1
        T = T + [T[i-1]/2]
        rk = rk + [1]
        r3b = r3b+[1]
        rcl = rcl+[1]
        rlim = rlim + [lim_3b(T[i],Ta)]
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
        rlim = rlim + [lim_3b(T[i],Ta)]
        
    #until now
    j=i
    rk_rm=rk[j]
    rcl_rm=rcl[j]
    r3b_rm = r3b[j]
    while T[i]>T0:
        i=i+1
        T = T + [T[i-1]*0.99]
        z.append(T[i]/T0 -1)
        rk = rk + [rk_rm*(1+1.7*10**9*rk_rm*r0*(T[j]**0.5*(np.log(np.sqrt(I*Ta/T[j]))+2)-T[i]**0.5*(np.log(np.sqrt(I*Ta/T[i]))+2)))**(-1)]
        D = (2*np.pi**2*g_s(T[i])/45)**2 * (4*np.pi*alpha)**5 * (2*np.sqrt(ma)*Ta**(9/2)/mb) / H(T[i])
        r3b = r3b +[r3b_rm*(1+D*r3b_rm*r3b_rm*r0*r0* 4/9 * (T[i]**(-9/2)-T[j]**(-9/2)))**(-1/2)]
        Dcl = Gamma(2-beta/2)*g_s(T[i])*s0*m_pl*(mu*Ta)**((beta-1)/2)/ ((beta-2.5)*2**((beta-3)/2)*np.sqrt(45*g_s(T[j])))/np.sqrt(t_rm)
        rcl = rcl + [rcl_rm*(1+Dcl*rcl_rm*r0*(T[i]**(2.5-beta)-T[j]**(2.5-beta)))**(-1)]
        rlim = rlim + [lim_3b(T[i],Ta)]       

    return [T,rk,r3b,rcl,rlim]
    

def contour_with_constraints(df, lang="eng", PICPATH="/home/kds/sci/threebody/pics/"):

    
    graph_name = "contour_with_constraints.png"

    # Установки для текста
    if lang == "rus":
        lab1 = 'Значения относительной плотности'
        lab2 = r'Относительная плотность прорекомбинировавших частиц: ${1-{r_{lim}}/{r_0}}$'
        title1 = 'Относительная плотность прорекомбинировавших частиц для различных параметров модели'
    elif lang == "eng":
        lab1 = "Relative density"
        lab2 = r'Limiting relative density: ${1-{r_{lim}}/{r_0}}$'
        title1 = "Relative density of recombined particles for different model parameters"

    # Создаем фигуру и оси графика
    fig, ax = plt.subplots(figsize=(12, 8))

    # Данные из DataFrame
    x = np.array(df.columns, dtype=float)
    y = np.array(df.index, dtype=float)
    z = df.values
    z = gaussian_filter(z, sigma=4)  # sigma определяет степень сглаживания


    # Основной контурный график
    colormap = plt.colormaps["coolwarm"].copy()
    cf = ax.contourf(x, y, z, levels=30, cmap=colormap, alpha=0.9)
    cbar = fig.colorbar(cf, ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    cbar.set_label(lab2, fontsize=20)

    ax.set_xlabel(r'$m_a$, eV', fontsize=20)
    ax.set_ylabel(r'$\alpha_y$', fontsize=20)

    # Вычисляем значения exp_constrait(alpha, m)
    X, Y = np.meshgrid(x, y)
    Z_constraint = exp_constrait(Y, X)
    B, A = 0.1, 1 
    ax.contour(X, Y, Z_constraint, levels=[B, A], colors='purple', linewidths=5, linestyles='-')  # Контурные линии
    constraint_mask = (Z_constraint >= B) & (Z_constraint <= A)
    ax.contourf(X, Y, constraint_mask, levels=[0.5, 1], colors='purple', alpha=0.3)


    ax.grid(False)
    ax.tick_params(labelsize=14)
    ax.title.set_size(16)
    cbar.ax.tick_params(labelsize=14)
    plt.tight_layout()
    plt.savefig(PICPATH + graph_name)
    plt.show()


def contourplot(df,lang="eng",PICPATH="/home/kds/sci/threebody/pics/",x0=10**8,y0=8):
    graph_name = "contourplot_lim.png"
    if lang == "rus":
        lab1 = 'Значения относительной плотности'
        lab2 = r'Относительная плотность прорекомбинировавших частиц: ${1-{r_{lim}}/{r_0}}$'
        title1 = 'Относительная плотность прорекомбинировавших частиц для различных параметров модели'
        lab3 = "На другом графике"
    elif lang == "eng":
        lab1 = "Reletive densety"
        lab2 = r'Limiting relative density: ${1-{r_{lim}}/{r_0}}$'
        lab3 = "on another plot"
        title1 = "Reletive density of recombined particles for different model parameters"

    x = np.array(df.columns, dtype=float)
    y = np.array(df.index, dtype=float)
    z = df.values
    z = gaussian_filter(z, sigma=4)

    fig, ax = plt.subplots(figsize=(12, 8))
    colormap = plt.colormaps["coolwarm"].copy()
    # Создаем контурный график с логарифмическим масштабом по обеим осям
    cf = ax.contourf(x, y, z, levels=30, cmap=colormap, alpha=0.9)
    cbar = fig.colorbar(cf, ax=ax)
    ax.set_xscale('log')
    ax.set_yscale('log')
    cbar.set_label(lab2,fontsize = 20)
    # cbar.set_label(r'Значения относительной плотности: $-\log{\frac{r}{r_0}}$',fontsize = 20)
    ax.set_xlabel(r'$m_a$, eV', fontsize=20)
    ax.set_ylabel(r'$\alpha_y$', fontsize=20)


    levels = cf.levels
    levels_done = []
    for i, level in enumerate(levels):
        if i == 0:
            continue
        paths = cf.collections[i-1].get_paths()
        for path in paths:
            x, y = path.vertices[2:, :].mean(0)
            if i%5 ==0 and not(level in levels_done):
                plt.text(x, y, f'{level:.2f}', fontsize=10, color='black')
                levels_done.append(level)
    # Добавляем сглаживание краев контуров
    cf.cmap.set_under('white')
    ax.grid(False) 
    # ax.grid(False, which='both', linestyle='-', color='grey', alpha=0.2)
    # Устанавливаем размеры шрифтов на осях и заголовке
    ax.tick_params(labelsize=14)
    ax.title.set_size(16)
    cbar.ax.tick_params(labelsize=14)
    ax.scatter(x0, y0, s=100, color='black', marker='*')
    # Добавляем подпись для точки
    ax.annotate(lab3, xy=(x0, y0), xytext=(x0+0.5, y0+0.2), fontsize=12)
    # Устанавливаем отступы между подписями и графиком
    plt.tight_layout()
    plt.savefig(PICPATH+graph_name)
    plt.show(block=False)

def lim_graph(T,rk,r3,rcl,rlim,lang="eng",PICPATH="/home/kds/sci/threebody/pics/"):
    graph_name = "lim_graph.png"
    if lang == "eng":
        lab = "Applicability range"
        lab1 = "Radiative recombinaton"
        lab2 = "Three-body recombinaton"
        lab3 = "Classical recombination"
    elif lang=="rus":
        lab = "Предел применимости"
        lab1 = "Радиационная рекомбинация"
        lab2 = "Трехчастичная рекомбинация"
        lab3 = "Классическая рекомбинация"

    krit = np.argmin(np.abs(np.array(r3) - np.array(rlim)))
    # graphplotter used to plot r(T) graphs
    fig, ax = plt.subplots(1,1, figsize = (8,6))

    ax.plot(T,r3,'b', 
            # label ="Трехчастичная рекомбинация")
            label = lab2)
    ax.plot(T,rlim, "-g",
            label =lab)
            # label ="Предел применимости")
    ax.set_xlabel('T, eV', fontsize=16)
    ax.set_ylabel(r'$r/r_0$', fontsize=16)

    ax.text(10**(0.9+np.ceil(np.log10(T[krit]))), 10 ** (-3),
            r'$m_a = 10^{%0.0f}$ eV' % (np.log10(ma)) + "\n" +
            r'$m_b = 10^{%0.0f}$ eV' % (np.log10(mb)) + "\n" +
            r'$\alpha = $' + repr(alpha),
            fontsize=12, bbox={'facecolor': '#FFFFFF', 'alpha': 0.9, 'pad': 10})

    ax.set_xlim([ 10**(1+np.ceil(np.log10(T[krit])))  ,  10**(np.floor(np.log10(T[krit])))])
    ax.set_ylim([10 ** (-4), 10 ** (1)])
    ax.tick_params(labelsize=14)

    ax.grid(color='lightgray', linestyle=':', linewidth=0.5, alpha=0.7)
    ax.legend(fontsize=14)
    ax.fill_between(T, rlim, ax.get_ylim()[1], alpha=0.3, color='green')
    

    plt.xscale('log')
    plt.yscale('log')
    
    plt.tight_layout()  # Оптимизирует расположение элементов
    plt.savefig(PICPATH + graph_name, dpi=300, bbox_inches='tight')  # Сохраняем в высоком разрешении
    plt.show(block=False)

def classkram3body_graph(T,rk,r3,rcl,rlim,lang="eng",PICPATH="/home/kds/sci/threebody/pics/"):
    graph_name = "classkramthreebody.png"
    if lang == "eng":
        lab1 = "Radiative recombinaton"
        lab2 = "Three-body recombinaton"
        lab3 = "Classical recombination"
    elif lang=="rus":
        lab1 = "Радиационная рекомбинация"
        lab2 = "Трехчастичная рекомбинация"
        lab3 = "Классическая рекомбинация"
    # graphplotter used to plot r(T) graphs
    fig, ax = plt.subplots(1,1, figsize = (8,6))
    ax.plot(T,r3,'b',
            label = lab2)
            # label ="Трехчастичная рекомбинация")
    ax.plot(T,rk,'r', 
            # label ="Радиационная рекомбинация")
            label = lab1)
    ax.plot(T,rcl,'g', 
            label = lab3)
            # label ="Классическая рекомбинация")
    ax.set_xlabel('T, eV', fontsize=16)
    ax.set_ylabel(r'$r/r_0$', fontsize=16)

    ax.text(10**(8), 10**(-11),
            r'$m_a = 10^{%0.0f}$ eV'%(np.log10(ma))+"\n"+
            r'$m_b = 10^{%0.0f}$ eV'%(np.log10(mb))+"\n"+
            r'$\alpha = $'+repr(alpha), 
            fontsize=12, bbox={'facecolor': '#FFFFFF', 'alpha': 0.9, 'pad': 10})
    ax.tick_params(labelsize=12)
    ax.grid(color='gray', linestyle='--', linewidth=1, alpha=0.5)

    leg = ax.legend(fontsize=16)
    leg.get_frame().set_linewidth(0.0)
    leg.get_frame().set_facecolor('#FFFFFF')
    leg.get_frame().set_alpha(0.9)
    ax.invert_xaxis()
    plt.semilogx()
    plt.semilogy()
    plt.box(False)
    plt.savefig(PICPATH+graph_name)
    plt.show(block=False)


def kramers_graph(T,rk,r3,rcl,rlim,lang="eng",PICPATH="/home/kds/sci/threebody/pics/"):
    graph_name = "kramers_graph.png"
    if lang == "eng":
        lab1 = "Radiative recombinaton"
    elif lang=="rus":
        lab1 = "Радиационная рекомбинация"

    # graphplotter used to plot r(T) graphs
    fig, ax = plt.subplots(1,1, figsize = (8,6))

    ax.plot(T,rk,'r', 
            label = lab1)
    ax.plot(T,[1]*len(T),'--b')
    ax.set_xlabel('T, eV', fontsize=16)
    ax.set_ylabel(r'$r/r_0$', fontsize=16)
    ax.set_ylim([0.94, 1.005])
    ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:.2f}'))
    ax.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.9)
    ax.legend(fontsize=16)
    ax.invert_xaxis()
    plt.semilogx()
    plt.box(False)
    plt.savefig(PICPATH+graph_name)
    plt.show(block=False)


#LOGIC -- 1 download data, else generate new
logic = 1
lang = "eng"
PICPATH = "/home/kds/sci/threebody/pics/"
err=10**(-3)
points = 100
mb_range = 10 #=mb/ma

def main():
    global ma,alpha,mb,lang,PICPATH
    ma = 100*10**9 #ev
    mb = ma*mb_range #ev
    alpha = 1/10 #none

    ali = 1/((ma/mb)**(3/2)*alpha**3*(2*np.pi**2*g_s(T0)/45))
    print(ali)
    alp_range = np.logspace(-1,1,points)
    m_range = np.logspace(11,15,points)
    if logic == 1:
        try:
            name = "data_nolog_"+str(points)+".scv"
            df = pd.read_csv(name)
        except FileNotFoundError:
            print(u'Saved data not found')
            df = kriteria(alp_range,m_range,mb_range)
    elif logic == 0:
        kriteria(alp_range,m_range,mb_range)
    else:
        pass
        
    [T,rk,r3,rcl,rlim] = calc(alpha,ma,mb)
    # kramers_graph(T,rk,r3,rcl,rlim,lang)
    # classkram3body_graph(T,rk,r3,rcl,rlim,lang)
    # lim_graph(T,rk,r3,rcl,rlim,lang,PICPATH)
    df = pd.read_csv(name, index_col=0)

    alpha = 1
    ma = 10**14
    mb = ma*mb_range
    contour_with_constraints(df,lang,PICPATH)
    contourplot(df,lang,PICPATH,ma,alpha)
    # [T,rk,r3,rcl,rlim] = calc(alpha,ma,mb)
    # lim_graph(T,rk,r3,rcl,rlim,lang,PICPATH)
    # classkram3body_graph(T,rk,r3,rcl,rlim,lang)

    plt.show()
    # [T,rk,r3,rcl,rlim] = calc(alpha,ma,mb)
    

main()