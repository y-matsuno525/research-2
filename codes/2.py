import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

epsilon = 1
p = epsilon
alpha = 1
beta = 0
#gamma = 1
#m = 1
zeta = 0
k = np.linspace(-np.pi/epsilon,np.pi/epsilon,100) #波数。分散関係を見るので、これは変えない。
m = np.linspace(0.0001,3,100) #####奥行きのパラメータ
gamma_list = np.linspace(10**(-24),3,100) #####アニメーションで変化させる変数

fig = plt.figure(figsize=(10,10)) #figureインスタンスの作成
ax = fig.add_subplot(projection='3d') #3つ目の軸を追加
K,M = np.meshgrid(k,m) #####格子点を作る。１つ目は波数kで固定。２つ目は奥行きのパラメータ。

def draw_frame(gamma): #####引数はアニメーションで変化させる変数

    ax.cla()

    #軸にラベルをつける
    ax.set_xlabel('k', labelpad=10) #固定
    ax.set_ylabel('m', labelpad=10) #####奥行きのパラメータ
    ax.set_zlabel('Energy', labelpad=10) #固定

    A = -1/(2*epsilon)*(p*np.cos(K)-beta*np.sin(K)-(p-epsilon*M*alpha))
    B = -1/(2*epsilon)*(alpha/gamma*np.cos(zeta)*1j*np.sin(K)+alpha/gamma*np.sin(zeta)*np.sin(K))
    C = -1/(2*epsilon)*(-1*alpha/gamma*np.cos(zeta)*1j*np.sin(K)+alpha/gamma*np.sin(zeta)*np.sin(K))
    D = -1/(2*epsilon)*(-1*p*np.cos(K)-beta*np.sin(K)+(p-epsilon*M*alpha))
    
    #分散関係
    E1 = (A+D)/2+np.sqrt(((A-D)/2)**2+np.abs(B)**2)
    E2 = (A+D)/2-np.sqrt(((A-D)/2)**2+np.abs(B)**2)

    
    ax.plot_surface(K,M,E1, cmap='jet',label="gamma="+str(gamma)) #####２つ目は奥行きのパラメータ,labelはアニメーションで変えるパラメータ
    ax.text2D(0.05, 0.95, f"gamma={gamma}", transform=ax.transAxes) #####
    ax.plot_surface(K,M,E2, cmap='jet') #####２つ目は奥行きのパラメータ
#plt.plot(k,E1,label = "m="+str(m))
#plt.plot(k,E2,label = "m="+str(m))

#plt.legend()
#ax.set_title("ε=1, α=1, γ1=10^-1")
#plt.grid()
ani = FuncAnimation(fig, draw_frame, frames=gamma_list, interval=100) #####3つめにアニメーションで変化させる変数のリストを設定
#ani.save("7.mp4", writer="ffmpeg")
#ax.legend()
plt.show()
