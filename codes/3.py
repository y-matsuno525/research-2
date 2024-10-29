import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

epsilon = 1
p = epsilon
#alpha = 1
beta = 0
gamma = 1
#m = 0.5
zeta = 0
k = np.linspace(-np.pi/epsilon,np.pi/epsilon,100) #波数。分散関係を見るので、これは変えない。
alpha = np.linspace(0.1,3,100) #奥行きのパラメータ
m_list = np.linspace(0.1,2,200) #アニメーションで変化させる変数

fig = plt.figure(figsize=(10,10)) #figureインスタンスの作成
ax = fig.add_subplot(projection='3d') #3つ目の軸を追加
K,Alpha = np.meshgrid(k,alpha) #格子点を作る。１つ目は波数kで固定。２つ目は奥行きのパラメータ。

#軸にラベルをつける
ax.set_xlabel('k', labelpad=10) #固定
ax.set_ylabel('alpha', labelpad=10) #奥行きのパラメータ
ax.set_zlabel('Energy', labelpad=10) #固定


def draw_frame(m): #引数はアニメーションで変化させる変数

    
    ax.cla()

    #軸にラベルをつける
    ax.set_xlabel('k', labelpad=10) #固定
    ax.set_ylabel('alpha', labelpad=10) #奥行きのパラメータ
    ax.set_zlabel('Energy', labelpad=10) #固定

    A = -1/(2*epsilon)*(p*np.cos(K)-beta*np.sin(K)-(p-epsilon*m*Alpha))
    B = -1/(2*epsilon)*(Alpha/gamma*np.cos(zeta)*1j*np.sin(K)+Alpha/gamma*np.sin(zeta)*np.sin(K))
    C = -1/(2*epsilon)*(-1*Alpha/gamma*np.cos(zeta)*1j*np.sin(K)+Alpha/gamma*np.sin(zeta)*np.sin(K))
    D = -1/(2*epsilon)*(-1*p*np.cos(K)-beta*np.sin(K)+(p-epsilon*m*Alpha))
    
    #分散関係
    E1 = (A+D)/2+np.sqrt(((A-D)/2)**2+np.abs(B)**2)
    E2 = (A+D)/2-np.sqrt(((A-D)/2)**2+np.abs(B)**2)

    
    ax.plot_surface(K,Alpha,E1, cmap='jet',label="m="+str(m))
    ax.text2D(0.05, 0.95, f"m={m}", transform=ax.transAxes)
    ax.plot_surface(K,Alpha,E2, cmap='jet')
#plt.plot(k,E1,label = "m="+str(m))
#plt.plot(k,E2,label = "m="+str(m))

#plt.legend()
#ax.set_title("ε=1, α=1, γ1=10^-1")
#plt.grid()
ani = FuncAnimation(fig, draw_frame, frames=m_list, interval=100) #3つめにアニメーションで変化させる変数のリストを設定
ani.save("3.mp4", writer="ffmpeg")
#ax.legend()
plt.show()
