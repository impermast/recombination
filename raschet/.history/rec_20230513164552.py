
from math import gamma as Gamma
import numpy as np
from numpy import sqrt as sqrt
import matplotlib.pyplot as plt
from scipy import integrate 


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


global r_t1,r_t2,r_t3
r_t1, r_t3,r_t2 = r0,r0,r0
print("Starting")

def g_s(T):
    return 2+2*8*(T>200*10**6)+ 7/8*(2*2*(T>0.5/3*10**6))+ 2*2*(T>100/3*10**6) + 2*2*(T>2/3*10**9) + 2*1*3*(T>=10**6)+2*1*3*4/11*(T<10**6)+ 2*2*3*2*(T>200*10**6) + 2*2*3*(T>200*10**6) + 2*2*3*(T>1.8/3*10**9) + 2*2*3*(T>4.5/3*10**9)+ 7/8*2*2*3*(T>171/3*10**9)
def g_e(T):
    return 2+2*8*(T>200*10**6)+ 7/8*(2*2*(T>0.5/3*10**6))+ 2*2*(T>100/3*10**6) + 2*2*(T>2/3*10**9) + 2*1*3*(T>=10**6)+ 2*1*3*(4/11)**(4/3)*(T<10**6)+ 2*2*3*2*(T>200*10**6) + 2*2*3*(T>200*10**6) + 2*2*3*(T>1.8/3*10**9) + 2*2*3*(T>4.5/3*10**9)+ 7/8*2*2*3*(T>171/3*10**9)


#this for three-body rec
def H(t):
    # Hubble const on RD and MD
    if t < t_rm:
        #MD
        return h_md/m_pl**(1/2)
    else:
        return 5.5*np.sqrt(g_e(t_rm)/11)/m_pl
def r(t):
    global r3b_rm
    D = (2*np.pi**2*g_s(t)/45)**2 * (4*np.pi*alpha)**5 * (2*sqrt(ma)*Ta**(9/2)/mb) / H(t)
    if I>t>=t_rm:
        r_t1 = (1+ D*2/5*r0*r0*(t**(-5)-I**(-5)))**(-1/2)
        return r_t1
    elif t_rm>t:
        return  (1+D*r0*r0* 4/9 * (t**(-9/2)-t_rm**(-9/2)))**(-1/2)
    else:
        return 1
        
#limits on 3body
def lim_3b(x):
    r_pit = (Ta**2/(x**4)) / ((4*np.pi*alpha)**3 *r0*(2*np.pi**2*g_s(x)/45)*sqrt(ma/mb))
    return r_pit

#this for classical approach from https://arxiv.org/pdf/1506.03094.pdf
# def r_cl(t):
#     global r_t2
#     D = Gamma(2-beta/2)*g_s(t)*s0*m_pl*(mu*Ta)**((beta-1)/2)/ (2**((beta-3)/2)*np.sqrt(45*g_e(t_rm)))
#     if I>t>=t_rm:
#         #MD era
#         r_t2 = (1+(D/np.sqrt(g_s(t_rm)))*r0*(t**(-beta+1/2)-I**(-beta+1/2))/(beta-1/2))**(-1)
#         return r_t2
#     elif t_rm>t:
#         #RD era
#         print(r_t2)
#         return (1+D*r0*r_t2*(t**(-beta)-t_rm**(-beta))/beta)**(-1)
#     else:
#         return 1

def func(y,t):
    ans = (y**3*(2*np.pi**2*g_s(t)/45)**2/H(t)) * 2*sqrt(ma)*(4*np.pi*alpha)**5*Ta**(9/2)/(t**4*mb)
    print("Ans = ",ans,"   y = ",y, "  t = ",t)
    return ans
        
def graph(T,r1,r2,rlim):
    fig, ax = plt.subplots(1,1)
    ax.plot(T,r1,'b',label ="Трехчастичная рекомбинация")
    ax.plot(T,rlim, "--g",label ="Limit")
    ax.plot(T,r2,'r', label ="Радиационная рекомбинация")
    #ax.plot(r_ode.t,r_ode.y,'y', label ="Ode solved")
    ax.set_xlabel('T,  eV')
    ax.set_ylabel(r'${r}/{r_0}$', fontsize='large', fontweight='bold')
    ax.grid()
    ax.legend()
    ax.invert_xaxis()
    plt.semilogx()
    plt.semilogy()
    plt.show()
    print("graphs")


def main():
    global r3b_rm
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

        


    #solver ode
    #r_ode = integrate.odeint(func,1,T)



    print("Data")
    graph(T,r3b,rcl,rlim)
main()