import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma, factorial

d_eta_list = [1,0.2]

n_k_list = []

k_m = 90

for d_eta in d_eta_list:
    
    eta = np.linspace(-10,10,100)
    k_c = np.linspace(0,k_m,100)
    
    
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
    
    n_k_qft = np.abs(u_c2*v_c1*(gamma(gamma_p)*gamma(gamma_p-alpha_p-beta_p))/(gamma(gamma_p-alpha_p)*gamma(gamma_p-beta_p))-v_c2*u_c1*(gamma(gamma_m.conj())*gamma(-gamma_m.conj()+alpha_m.conj()+beta_m.conj()))/(gamma(alpha_m.conj())*gamma(beta_m.conj())))**2
    
    plt.plot(k_c,n_k_qft,label="QFT")
    
    L_list = [64,128,256,512]
    
    l = 5*2*np.pi
    
    for L in L_list:
        
        epsilon = l/L
        kappa = np.linspace(0,k_m*epsilon,100)
        
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
    
        n_k = np.abs(u_d2*v_d1*(gamma(gamma_p)*gamma(gamma_p-alpha_p-beta_p))/(gamma(gamma_p-alpha_p)*gamma(gamma_p-beta_p))-v_d2*u_d1*(gamma(gamma_m.conj())*gamma(-gamma_m.conj()+alpha_m.conj()+beta_m.conj()))/(gamma(alpha_m.conj())*gamma(beta_m.conj())))**2
        plt.plot(k_c,n_k,label="L="+str(L))
    
    plt.legend(loc="right")
    plt.grid()
    plt.show()
