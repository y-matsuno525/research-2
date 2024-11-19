import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma, factorial
from scipy.stats import linregress

d_eta_list = [1,0.2]
k_c_list = [0.5,1.1]
d_n_list = []

for i,d_eta in enumerate(d_eta_list):
    
    eta = np.linspace(-10,10,100)
    k_c = k_c_list[i]
    
    
    m = 1
    
    
    #Continuum theoryのa
    a_c1 = 1
    a_c2 = 2 
    
    
    #Continuum theoryのエネルギー
    epsilon_c1 = np.sqrt((m*a_c1)**2+k_c**2)
    epsilon_c2 = np.sqrt((m*a_c2)**2+k_c**2)
    
    #Gauss hypergeometric functionのパラメータ
    z_p = 1/2+np.tanh(eta/d_eta)/2 
    z_m = 1/2-np.tanh(eta/d_eta)/2
    alpha_p = 1j/2*(epsilon_c2+epsilon_c1-m*a_c1+m*a_c2)*d_eta
    alpha_m = 1j/2*(epsilon_c2-epsilon_c1-m*a_c1+m*a_c2)*d_eta
    beta_p = 1+1j/2*(epsilon_c2+epsilon_c1+m*a_c1-m*a_c2)*d_eta
    beta_m = 1+1j/2*(epsilon_c2-epsilon_c1+m*a_c1-m*a_c2)*d_eta
    gamma_p = 1+1j*epsilon_c1*d_eta
    gamma_m = 1-1j*epsilon_c1*d_eta
    
    #Continuum theoryの固有ベクトル
    u_c1 = 1/np.sqrt(2*epsilon_c1*(epsilon_c1-m*a_c1))*(-m*a_c1+epsilon_c1)
    v_c1 = 1/np.sqrt(2*epsilon_c1*(epsilon_c1-m*a_c1))*(1j*k_c)
    u_c2 = 1/np.sqrt(2*epsilon_c2*(epsilon_c2-m*a_c2))*(-m*a_c2+epsilon_c2)
    v_c2 = 1/np.sqrt(2*epsilon_c2*(epsilon_c2-m*a_c2))*(1j*k_c)
    
    n_k_qft = np.abs(u_c2*v_c1*(gamma(gamma_p)*gamma(gamma_p-alpha_p-beta_p))/(gamma(gamma_p-alpha_p)*gamma(gamma_p-beta_p))-v_c2*u_c1*(gamma(gamma_m.conjugate())*gamma(-gamma_m.conjugate()+alpha_m.conjugate()+beta_m.conjugate()))/(gamma(alpha_m.conjugate())*gamma(beta_m.conjugate())))**2
    
   
    
    L_list = [64,128,256,512]
    l = 5*2*np.pi
    
    for L in L_list:
        
        epsilon = l/L
        kappa = k_c*epsilon
        
        #Discrete theoryのa
        a_d1 = 1 - (1-np.cos(kappa))/epsilon*m
        a_d2 = 2 - (1-np.cos(kappa))/epsilon*m
            
        #Discrete theoryの固有値に出てくるパラメータ
        #z_1 = 1/epsilon*(1-np.cos(k)-epsilon*m*a_d1)
        #z_2 = 1/epsilon*(1-np.cos(k)-epsilon*m*a_d2)
        #y = 1/epsilon*np.sin(k)
        
        #Discrete theoryのエネルギー
        epsilon_d1 = np.sqrt((m*a_d1)**2+(np.sin(kappa)/epsilon)**2)
        epsilon_d2 = np.sqrt((m*a_d2)**2+(np.sin(kappa)/epsilon)**2)
    
        #Gauss hypergeometric functionのパラメータ
        z_p = 1/2+np.tanh(eta/d_eta)/2 
        z_m = 1/2-np.tanh(eta/d_eta)/2
        alpha_p = 1j/2*(epsilon_d2+epsilon_d1-m*a_d1+m*a_d2)*d_eta
        alpha_m = 1j/2*(epsilon_d2-epsilon_d1-m*a_d1+m*a_d2)*d_eta
        beta_p = 1+1j/2*(epsilon_d2+epsilon_d1+m*a_d1-m*a_d2)*d_eta
        beta_m = 1+1j/2*(epsilon_d2-epsilon_d1+m*a_d1-m*a_d2)*d_eta
        gamma_p = 1+1j*epsilon_d1*d_eta
        gamma_m = 1-1j*epsilon_d1*d_eta
        
        #Discrete theoryの固有ベクトル
        u_d1 = 1/np.sqrt(2*epsilon_d1*(epsilon_d1-m*a_d1))*(-m*a_d1+epsilon_d1)
        v_d1 = 1/np.sqrt(2*epsilon_d1*(epsilon_d1-m*a_d1))*(1j*(np.sin(kappa)/epsilon))
        u_d2 = 1/np.sqrt(2*epsilon_d2*(epsilon_d2-m*a_d2))*(-m*a_d2+epsilon_d2)
        v_d2 = 1/np.sqrt(2*epsilon_d2*(epsilon_d2-m*a_d2))*(1j*(np.sin(kappa)/epsilon))
    
        n_k = abs(u_d2*v_d1*(gamma(gamma_p)*gamma(gamma_p-alpha_p-beta_p))/(gamma(gamma_p-alpha_p)*gamma(gamma_p-beta_p))-v_d2*u_d1*(gamma(gamma_m.conjugate())*gamma(-gamma_m.conjugate()+alpha_m.conjugate()+beta_m.conjugate()))/(gamma(alpha_m.conjugate())*gamma(beta_m.conjugate())))**2
        d_n_list.append(abs(n_k-n_k_qft))

    # 対数スケールでの線形回帰
    log_L = np.log(L_list)
    log_error = np.log(d_n_list)
    slope, intercept, r_value, p_value, std_err = linregress(log_L, log_error)
        
    plt.scatter(L_list, d_n_list,label="d_eta="+str(d_eta)+", slope="+str(slope))
    d_n_list = []

plt.xscale("log",base=2)
plt.yscale("log")
plt.ylim(10**(-8),0.1)
plt.xlim(50,600)
plt.legend()
plt.grid()
plt.show()







import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma, factorial

d_eta_list = [1,0.2]
k_c_list = [7,10]
d_n_list = []

for i,d_eta in enumerate(d_eta_list):
    
    eta = np.linspace(-10,10,100)
    k_c = k_c_list[i]
    
    
    m = 1
    
    
    #Continuum theoryのa
    a_c1 = 1
    a_c2 = 2 
    
    
    #Continuum theoryのエネルギー
    epsilon_c1 = np.sqrt((m*a_c1)**2+k_c**2)
    epsilon_c2 = np.sqrt((m*a_c2)**2+k_c**2)
    
    #Gauss hypergeometric functionのパラメータ
    z_p = 1/2+np.tanh(eta/d_eta)/2 
    z_m = 1/2-np.tanh(eta/d_eta)/2
    alpha_p = 1j/2*(epsilon_c2+epsilon_c1-m*a_c1+m*a_c2)*d_eta
    alpha_m = 1j/2*(epsilon_c2-epsilon_c1-m*a_c1+m*a_c2)*d_eta
    beta_p = 1+1j/2*(epsilon_c2+epsilon_c1+m*a_c1-m*a_c2)*d_eta
    beta_m = 1+1j/2*(epsilon_c2-epsilon_c1+m*a_c1-m*a_c2)*d_eta
    gamma_p = 1+1j*epsilon_c1*d_eta
    gamma_m = 1-1j*epsilon_c1*d_eta
    
    #Continuum theoryの固有ベクトル
    u_c1 = 1/np.sqrt(2*epsilon_c1*(epsilon_c1-m*a_c1))*(-m*a_c1+epsilon_c1)
    v_c1 = 1/np.sqrt(2*epsilon_c1*(epsilon_c1-m*a_c1))*(1j*k_c)
    u_c2 = 1/np.sqrt(2*epsilon_c2*(epsilon_c2-m*a_c2))*(-m*a_c2+epsilon_c2)
    v_c2 = 1/np.sqrt(2*epsilon_c2*(epsilon_c2-m*a_c2))*(1j*k_c)
    
    n_k_qft = np.abs(u_c2*v_c1*(gamma(gamma_p)*gamma(gamma_p-alpha_p-beta_p))/(gamma(gamma_p-alpha_p)*gamma(gamma_p-beta_p))-v_c2*u_c1*(gamma(gamma_m.conjugate())*gamma(-gamma_m.conjugate()+alpha_m.conjugate()+beta_m.conjugate()))/(gamma(alpha_m.conjugate())*gamma(beta_m.conjugate())))**2
    
   
    
    L_list = [64,128,256,512]
    l = 5*2*np.pi
    
    for L in L_list:
        
        epsilon = l/L
        kappa = k_c*epsilon
        
        #Discrete theoryのa
        a_d1 = 1 #- (1-np.cos(kappa))/epsilon*m
        a_d2 = 2 #- (1-np.cos(kappa))/epsilon*m
            
        #Discrete theoryの固有値に出てくるパラメータ
        #z_1 = 1/epsilon*(1-np.cos(k)-epsilon*m*a_d1)
        #z_2 = 1/epsilon*(1-np.cos(k)-epsilon*m*a_d2)
        #y = 1/epsilon*np.sin(k)
        
        #Discrete theoryのエネルギー
        epsilon_d1 = np.sqrt((m*a_d1)**2+(np.sin(kappa)/epsilon)**2)
        epsilon_d2 = np.sqrt((m*a_d2)**2+(np.sin(kappa)/epsilon)**2)
    
        #Gauss hypergeometric functionのパラメータ
        z_p = 1/2+np.tanh(eta/d_eta)/2 
        z_m = 1/2-np.tanh(eta/d_eta)/2
        alpha_p = 1j/2*(epsilon_d2+epsilon_d1-m*a_d1+m*a_d2)*d_eta
        alpha_m = 1j/2*(epsilon_d2-epsilon_d1-m*a_d1+m*a_d2)*d_eta
        beta_p = 1+1j/2*(epsilon_d2+epsilon_d1+m*a_d1-m*a_d2)*d_eta
        beta_m = 1+1j/2*(epsilon_d2-epsilon_d1+m*a_d1-m*a_d2)*d_eta
        gamma_p = 1+1j*epsilon_d1*d_eta
        gamma_m = 1-1j*epsilon_d1*d_eta
        
        #Discrete theoryの固有ベクトル
        u_d1 = 1/np.sqrt(2*epsilon_d1*(epsilon_d1-m*a_d1))*(-m*a_d1+epsilon_d1)
        v_d1 = 1/np.sqrt(2*epsilon_d1*(epsilon_d1-m*a_d1))*(1j*(np.sin(kappa)/epsilon))
        u_d2 = 1/np.sqrt(2*epsilon_d2*(epsilon_d2-m*a_d2))*(-m*a_d2+epsilon_d2)
        v_d2 = 1/np.sqrt(2*epsilon_d2*(epsilon_d2-m*a_d2))*(1j*(np.sin(kappa)/epsilon))
    
        n_k = abs(u_d2*v_d1*(gamma(gamma_p)*gamma(gamma_p-alpha_p-beta_p))/(gamma(gamma_p-alpha_p)*gamma(gamma_p-beta_p))-v_d2*u_d1*(gamma(gamma_m.conjugate())*gamma(-gamma_m.conjugate()+alpha_m.conjugate()+beta_m.conjugate()))/(gamma(alpha_m.conjugate())*gamma(beta_m.conjugate())))**2
        d_n_list.append(abs(n_k-n_k_qft))

    # 対数スケールでの線形回帰
    log_L = np.log(L_list)
    log_error = np.log(d_n_list)
    slope, intercept, r_value, p_value, std_err = linregress(log_L, log_error)
        
    plt.scatter(L_list, d_n_list,label="d_eta="+str(d_eta)+", slope="+str(slope))
    d_n_list = []

plt.xscale("log",base=2)
plt.yscale("log")
plt.ylim(10**(-8),0.1)
plt.xlim(50,600)
plt.legend()
plt.grid()
plt.show()
