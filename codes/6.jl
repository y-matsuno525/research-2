using ITensors, ITensorMPS
let
    #基本設定
    L = 100 #サイト数
    ε = 1 #格子間距離
    sites = siteinds("Fermion", L) #サイトの定義
    t = -10^24
    m = 1

    #計量に出てくる関数の定義
    a2 = 2
    a1 = 1
    Δt = 1
    
    function α(t, j);
        return (a2+a1)/2+(a2-a1)/2*tanh(t/Δt)
    end

    function β(t, j);
        return 0
    end

    function γ(t, j);
        return (a2+a1)/2+(a2-a1)/2*tanh(t/Δt)
    end

    function p(t, j);
        return 1
    end

    function ζ(t, j);
        return 0
    end

    #空間微分用関数
    function dif(f, t, j);
        return (f(t, (j+1)*ε)-f(t, (j-1)*ε))/(2*ε)
    end

    #Localハミルトニアンの定義
    function hamiltonian(t)
        os += -1/(2*ε)*α(t, j)/γ(t, j)*cos(ζ(t, j)), "C", j+1, "C", j
        os += -1/(2*ε)*α(t, j)/γ(t, j)*cos(ζ(t, j)), "Cdag", j, "Cdag", j+1
        os += -1/(2*ε)*im*α(t, j)/γ(t, j)*sin(ζ(t, j)), "C", j+1, "C", j
        os -= -1/(2*ε)*im*α(t, j)/γ(t, j)*sin(ζ(t, j)), "Cdag", j, "Cdag", j+1
        os += -1/(2*ε)*p(t, j), "Cdag", j, "C", j+1
        os += -1/(2*ε)*p(t, j), "Cdag", j+1, "C", j
        os += -1/(2*ε)*im*β(t, j), "Cdag", j, "C", j+1
        os -= -1/(2*ε)*im*β(t, j), "Cdag", j+1, "C", j
        os += -1/(2*ε)*1/2*(2*p(t, j)-ε*(dif(p, t, j))+2*m*α(t, j)-dif(ζ, t, j)-β(t, j)*dif(ζ, t, j)), "Id", j
        os += -1/(2*ε)*1/2*(2*p(t, j)-ε*(dif(p, t, j))+2*m*α(t, j)-dif(ζ, t, j)-β(t, j)*dif(ζ, t, j))*(-2), "Cdag", j, "C", j
        
    #DMRG
    #Globalハミルトニアンを生成
    os = OpSum()
    for j=1:L-1
        os += hamiltonian(j, t)
    end
    H = MPO(os,sites)
        
    #初期状態をランダムに生成
    psi0 = random_mps(sites;linkdims=10)

    nsweeps = 5
    maxdim = [10,20,100,100,200]
    cutoff = [1E-10]

    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)

    ###時間発展###
    cutoff = 1E-8
    tau = 0.1
    ttotal = 5.0

    #gates作成
    function create_gates(t)
        gates = ITensor[]
        for j in 1:(L-1)
            hj = hamiltonian(t, j)
            Gj = exp(-im * tau / 2 * hj)
            push!(gates, Gj)
        end

    return
end

#生成消滅演算子C,Cdagのフーリエ変換
#量子数nのドメイン生成
n_domain = []
for n=-L/2+1:L/2
    push!(n_domain,n)

#波数kのドメイン生成
k_domain = []
for n in n_domain
    k = (2*π)/L*(n-1/2)
    push!(k_domain,k)
        
#C_kを格納する配列
c_k_list = []

#各波数でフーリエ変換を実行
for k in k_domain
    c_k = OpSum()
    for j=1:L
        c_k += 1/(√L)*exp(-im*k*j), "C", j
    push!(c_k_list, c_k)
#c_nの再定義(独自演算子の定義でやるほうがいい)
c_n_list = []
for j=1:L
    c_n = OpSum()
    for (i,k) in enumerate(k_domain)
        c_n += 1/(√L)*exp(-im*k*j)*c_k_list[i]
#c_n_listには波数表示（つまりフーリエ変換後）のc_nが格納されている