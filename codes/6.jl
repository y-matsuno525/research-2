using ITensors, ITensorMPS
using Plots

let
    ###基本設定#########################################################
    L = 100 #サイト数
    ε = 1 #格子間距離
    sites = siteinds("S=1/2", L) #サイトの定義
    t = 0
    m = 1
    τ = 0.1 #時間ステップの間隔
    ttotal = 5　#終時刻

    ###計量に出てくる関数の定義##########################################
    a2 = 2
    a1 = 1
    Δt = 1
    x_h = 50.5
    
    function α(t, j);
            return 1#(a2+a1)/2+(a2-a1)/2*tanh(t/Δt)
    end

    function β(t, j);
        return √((x_h-2*t)*ε/(j*ε))
    end

    function γ(t, j);
        return -1
    end

    function p(t, j);
        return 1
    end

    function ζ(t, j);
        return 0
    end

    

    ###ハミルトニアンの定義###############################################
    #OpSumを用いたLocalハミルトニアンの定義
    function hamiltonian_opsum(t, j)
        os_tmp = OpSum()
        os_tmp += -1/(4*ε)*(α(t, j)/γ(t, j)*cos(ζ(t, j))+p(t, j)), "Sx", j, "Sx", j+1
        os_tmp -= -1/(4*ε)*(α(t, j)/γ(t, j)*cos(ζ(t, j))-p(t, j)), "Sy", j, "Sy", j+1
        os_tmp -= -1/(4*ε)*(β(t, j)+α(t, j)/γ(t, j)*sin(ζ(t, j))), "Sx", j, "Sy", j+1
        os_tmp += -1/(4*ε)*(β(t, j)-α(t, j)/γ(t, j)*sin(ζ(t, j))), "Sy", j, "Sx", j+1
        os_tmp += -1/(4*ε)*(2*p(t, j)-ε*(dif_x(p, t, j)+2*m*α(t, j)-dif_t(ζ, t, j)-β(t, j)*dif_x(ζ, t, j))), "Sz", j 
        return os_tmp
    end

    #opを用いたLocalハミルトニアンの定義
    function hamiltonian_op(t, j)
        hj_tmp = ITensor()
        hj_tmp += -1/(4*ε)*(α(t, j)/γ(t, j)*cos(ζ(t, j))+p(t, j)) * op("Sx", sites[j]) * op("Sx", sites[j+1])
        hj_tmp -= -1/(4*ε)*(α(t, j)/γ(t, j)*cos(ζ(t, j))-p(t, j)) * op("Sy", sites[j]) * op("Sy", sites[j+1])
        hj_tmp -= -1/(4*ε)*(β(t, j)+α(t, j)/γ(t, j)*sin(ζ(t, j))) * op("Sx", sites[j]) * op("Sy", sites[j+1])
        hj_tmp += -1/(4*ε)*(β(t, j)-α(t, j)/γ(t, j)*sin(ζ(t, j))) * op("Sy", sites[j]) * op("Sx", sites[j+1])
        hj_tmp += -1/(4*ε)*(2*p(t, j)-ε*(dif_x(p, t, j)+2*m*α(t, j)-dif_t(ζ, t, j)-β(t, j)*dif_x(ζ, t, j))) * op("Sz", sites[j]) *op("Id", sites[j+1])
        return hj_tmp
    end
        
    ###その他関数#############################################################3
    #gates作成
    function create_gates(t)
        gates_tmp = ITensor[]
        for j in 1:(L-1)
            hj = hamiltonian_op(t, j)
            Gj = exp(-im * τ / 2 * hj)
            push!(gates_tmp, Gj)
        end
        append!(gates_tmp, reverse(gates_tmp))
        return gates_tmp
    end

    #空間微分用関数
    function dif_x(f, t, j);
        return (f(t, (j+1)*ε)-f(t, (j-1)*ε))/(2*ε)
    end

    #時間微分用関数
    function dif_t(f, t, j);
        return (f(t+τ, j*ε)-f(t-τ, j*ε))/(2*τ)
    end

    ###計算###################################################################
    

    #DMRGでearly timeでの基底状態を求める
    #Globalハミルトニアンを生成
    os = OpSum()
    for j=1:L-1
        os += hamiltonian_opsum(t, j)
    end
    H = MPO(os,sites)

    println("ハミルトニアン生成完了")
        
    #初期状態をランダムに生成
    psi0 = random_mps(sites;linkdims=10)
    
    #DMRGの設定
    nsweeps = 5
    maxdim = [10,20,100,100,200]
    cutoff_DMRG = [1E-10]

    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff=cutoff_DMRG)

    println("DMRG完了")

    #TEBDで時間発展を計算
    sites = s = siteinds(psi) #参考：https://itensor.org/support/3645/the-time-evolution-for-a-h5-file
    Sx_t = []
    cutoff_TEBD = 1E-20
    for (i,t) in enumerate(0.0:τ:ttotal)
        Sx_n = []
        #Sxの期待値をearly timeの基底状態をtまで時間発展させたもので挟む
        for j in 1:L
            Sx = expect(psi, "Sx"; sites=j)
            push!(Sx_n, Sx)
        end
        push!(Sx_t, Sx_n)

        gates = create_gates(t)

        psi = apply(gates, psi; cutoff=cutoff_TEBD)
        normalize!(psi)
        println(string(round(i/(ttotal/τ)*100, digits=1))*"%完了")
    end

    println("時間発展完了")

    ###プロット#################################################################

    anim = Animation()
    sites_array = collect(Float32, 1:L)
    
    for (i, sx_n) in enumerate(Sx_t)
        plt = plot(sites_array, sx_n,ylim=(-0.5, 0.5) , title="t = "*string(round(i*τ, digits=1)), label="")
        frame(anim, plt)
    end

    gif(anim, "plot.gif", fps=20)
end