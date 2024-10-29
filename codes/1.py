import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

fig, ax = plt.subplots()
epsilon = 1
p = 1
alpha_list = []
beta = 0
gamma = 1
zeta = 1
k = np.linspace(-np.pi / epsilon, np.pi / epsilon, 100)

# m_listとalpha_listの設定
m_list = np.concatenate([np.linspace(2, 0, 100), np.linspace(0, 2, 150)])
alpha_list = [0] * 100 + [2 * m / (1 + m**2) for m in np.linspace(0, 2, 150)]

# アニメーションフレームの更新関数
def update(i):
    ax.clear()
    ax.set_xlim(-np.pi / epsilon, np.pi / epsilon)
    ax.set_ylim(-2, 2)
    ax.set_xlabel("k")
    ax.set_ylabel("Energy")
    ax.set_title(f"m={m_list[i]:.2f}")

    A = -1/(2*epsilon)*(p*np.cos(k)-beta*np.sin(k)-(p-epsilon*m_list[i]*alpha_list[i]))
    B = -1/(2*epsilon)*(alpha_list[i]/gamma*np.cos(zeta)*1j*np.sin(k)+alpha_list[i]/gamma*np.sin(zeta)*np.sin(k))
    D = -1/(2*epsilon)*(-1*p*np.cos(k)-beta*np.sin(k)+(p-epsilon*m_list[i]*alpha_list[i]))

    # エネルギー分散関係の計算
    E_1 = (A + D) / 2 + np.sqrt(((A - D) / 2)**2 + np.abs(B)**2)
    E_2 = (A + D) / 2 - np.sqrt(((A - D) / 2)**2 + np.abs(B)**2)

    # プロット
    ax.plot(k, E_1.real, color="blue", label="E_1")
    ax.plot(k, E_2.real, color="red", label="E_2")
    ax.legend()

# アニメーションの作成
ani = animation.FuncAnimation(fig, update, frames=len(m_list), interval=10, blit=False)
ani.save("5.mp4", writer="ffmpeg")
plt.show()
