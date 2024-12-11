#=
計算の流れ
1. DMRGで(3.17)の十分過去での基底状態を計算する

やること
1. (3.17)のハミルトニアンを作る
=#

using ITensors, ITensorMPS

let
    
  N = 100
  sites = siteinds("S=1/2",N)
  ζ = 0
  m = 1
  ε = 1


    
  #空間微分関数の定義（中央差分,最終サイトがN-1）,そもそもこの差分のとり方でいいか？
  function dif(f, n, t = nothing);
    h = 1e-5

    #時間を引数に持つ場合
    if !t == nothing
        if n == 0
            return (f(n + 1,t) - f(n ,t)) / h
        elseif n = N-1
            return (f(n ,t) - f(n - 1,t)) / h
        else
            return (f(n + 1,t) - f(n - 1,t)) / (2 * h)
        end
    #時間を引数に持たない場合
    else
        if n == 0
            return (f(n + 1) - f(n) / h
        elseif n = N-1
            return (f(n) - f(n - 1)) / h
        else
            return (f(n + 1) - f(n - 1)) / (2 * h)
        end
    end
  end


    
  #計量内関数 α,β,γ　の定義
  #α = 0.1 * x + t
　function α(n,t);
        return 0.1 * (n * ε) + t
  end
  #β = 0.05 * x + t
  function β(n,t);
        return 0.05 * (n * ε) + t
  end
  #γ = 0.075 * x + t
  function γ(x,t);
        return 0.075 * (n * ε) + t
  end

  function p(n,t);
        return α(n,t)/γ(n,t)
  end


    
  #ハミルトニアン(ζは定数)
  function H_t(N,t);
    os = OpSum()
    for j=1:N-1
      os += (α(j,t)/γ(j,t)*cos(ζ) + p(j,t)) * "X",j,"X",j+1
      os += -(α(j,t)/γ(j,t)*cos(ζ) - p(j,t)) * "Y",j,"Y",j+1
      os += -(β(j,t) + α(j,t)/γ(j,t)*sin(ζ)) * "X",j,"X",j+1
      os += (β(j,t) - α(j,t)/γ(j,t)*sin(ζ)) * "Y",j,"Y",j+1
      os += (2*p(j,t)-ε(dif(p(j,t),j * ε,t) + 2 * m * α(j,t))) * "Z",j
    end
    H = MPO(os,sites)
    return H
  end


    
  #DMRG実行
  psi0 = random_mps(sites;linkdims=10)

  nsweeps = 5
  maxdim = [10,20,100,100,200]
  cutoff = [1E-10]
　
  energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)

  return
end