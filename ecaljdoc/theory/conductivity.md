# 電気伝導度計算の基本 {#conductivity}

::: info 下書き
**一体問題** として、線形応答（久保公式）から不純物による電気抵抗を導く。
ミソは、不純物の配置について平均を取るのが **グリーン関数の積** だという点である（一本ずつ平均してから掛けるのでは足りない）。
グリーン関数の符号の約束は [Basic I](./basic) と同じで、遅延グリーン関数は $G^R(E)=1/(E-H+i\delta)$。
遅延型の応答関数については [付録](./appendix#causality) も見よ。
:::

## 設定 {#setting}

相互作用の無い電子（質量 $m$、$\varepsilon_{\bf k}=\hbar^2k^2/2m$）が、体積 $\Omega$ の中に乱雑に置かれた不純物のポテンシャル

$$
V({\bf r})=\sum_{j} u({\bf r}-{\bf R}_j)
\tag{C1}
$$

を感じるとする。不純物の数密度を $n_i$、$u$ のフーリエ変換を $u({\bf q})=\int d^3r\, e^{-i{\bf q}\cdot{\bf r}}u({\bf r})$ と書く。
観測される量は、不純物の配置 $\{ {\bf R}_j \}$ についての平均 $\langle\cdots\rangle_{\rm imp}$ である。

**ここで扱うのは一体問題である。** 電子どうしの相互作用は無く、各電子は不純物の作る静的なポテンシャル (C1) の中を独立に動く。
ハミルトニアンは一体の $H=\hat{\bf p}^2/2m+V({\bf r})$ で、以下の $G^R(E)=1/(E-H+i\delta)$ はその一体のグリーン関数（レゾルベント）である。
不純物の配置を一つ決めれば、問題は一体のシュレーディンガー方程式を解くことに尽きる。
以下に出てくる自己エネルギーやバーテックス補正は、多体のダイヤグラムに似ているが、すべて不純物の配置についての平均から来るものである（電子の間の相互作用から来るのではない）。

## 久保公式 {#kubo}

電場 ${\bf E}(t)$ を摂動として加え、電流密度の応答を一次まで取ると、伝導度は電流と電流の遅延相関関数で書ける（久保公式）。
相互作用の無い電子の直流（$\omega\to 0$）、低温の伝導度は、フェルミ準位 $E_F$ での一体の量だけで書ける（Kubo–Greenwood の式）:

$$
\sigma_{xx}=\frac{2\hbar e^2}{\pi\Omega}\,
\Big\langle {\rm Tr}\big[\hat v_x\,{\rm Im}\,G^R(E_F)\,\hat v_x\,{\rm Im}\,G^R(E_F)\big]\Big\rangle_{\rm imp}.
\tag{C2}
$$

$\hat v_x$ は速度の演算子、${\rm Tr}$ は一体の状態についての和、因子 2 はスピンである。
${\rm Im}\,G^R(E)=-\pi\delta(E-H)$ を入れると、固有状態 $n$ についての和 $\sigma_{xx}=(2\pi\hbar e^2/\Omega)\sum_{nm}|\langle n|\hat v_x|m\rangle|^2\delta(E_F-E_n)\delta(E_F-E_m)$ になる。

${\rm Im}\,G^R=(G^R-G^A)/2i$ と書いて展開すると、$G^RG^A$ の項と $G^RG^R$、$G^AG^A$ の項が出る。
不純物が弱いとき主になるのは $G^RG^A$ の項で、

$$
\sigma_{xx}\simeq\frac{\hbar e^2}{\pi\Omega}\,{\rm Re}\,
\Big\langle {\rm Tr}\big[\hat v_x\,G^R(E_F)\,\hat v_x\,G^A(E_F)\big]\Big\rangle_{\rm imp}
\tag{C3}
$$

である（$G^RG^R$ の項は、$1/\tau$ について高次の補正を与えるだけである）。

**ミソ**: 平均を取るのは $G^R$ と $G^A$ の **積** である。

$$
\big\langle \hat v_x G^R \hat v_x G^A\big\rangle_{\rm imp}\neq
\hat v_x\langle G^R\rangle_{\rm imp}\,\hat v_x\langle G^A\rangle_{\rm imp}.
\tag{C4}
$$

$G^R$ と $G^A$ が同じ不純物で散乱される相関（ダイヤグラムでは、二本の線を結ぶ不純物の線）が残る。これがバーテックス補正である。

## 平均した一本のグリーン関数: 寿命 {#averaged-g}

一本の $G^R$ を平均すると並進対称性が戻り、波数 ${\bf k}$ で対角になる。

$$
\langle G^R\rangle_{\rm imp}({\bf k},E)=\frac{1}{E-\varepsilon_{\bf k}-\Sigma^R({\bf k},E)}.
\tag{C5}
$$

不純物について二次（Born 近似）の自己エネルギーは

$$
\Sigma^R({\bf k},E)=n_i\int\frac{d^3k'}{(2\pi)^3}\,
\frac{|u({\bf k}-{\bf k}')|^2}{E-\varepsilon_{\bf k'}+i\delta}
\tag{C6}
$$

で、その虚部が準粒子の寿命 $\tau$ を与える（実部はエネルギーのずれで、$E_F$ に繰り込む）:

$$
\begin{aligned}
&{\rm Im}\,\Sigma^R=-\frac{\hbar}{2\tau},\\
&\frac{1}{\tau}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,\delta(E-\varepsilon_{\bf k'}).
\end{aligned}
\tag{C7}
$$

$1/\tau$ はフェルミの黄金則による散乱の速さそのものである（[Basic I](./basic#td-vs-ti-perturbation) の、黄金則と自己エネルギーの虚部の対応）。

## バーテックス補正を入れない場合: Drude の式 {#drude}

式 (C3) の平均を、一本ずつの平均 (C5) の積で置きかえると（バーテックス補正を落とすと）、${\bf k}$ についての和になる。
$\tau$ が長い（$\hbar/\tau\ll E_F$）とき

$$
\begin{aligned}
\langle G^R\rangle({\bf k},E)\,\langle G^A\rangle({\bf k},E)
&=\frac{1}{(E-\varepsilon_{\bf k})^2+(\hbar/2\tau)^2}\\
&\simeq\frac{2\pi\tau}{\hbar}\,\delta(E-\varepsilon_{\bf k})
\end{aligned}
\tag{C8}
$$

なので、$v_x^2$ の角度平均 $v_F^2/3$ と、スピンあたりの状態密度 $N_1(E_F)$（単位体積あたり）を使って

$$
\sigma_{xx}=\frac{\hbar e^2}{\pi}\cdot\frac{2\pi\tau}{\hbar}\cdot\frac{v_F^2}{3}\,N_1(E_F)
=\frac{2e^2\tau v_F^2N_1(E_F)}{3}=\frac{ne^2\tau}{m}
\tag{C9}
$$

となる。最後は自由電子の $N_1(E_F)=mk_F/(2\pi^2\hbar^2)$、$n=k_F^3/(3\pi^2)$、$v_F=\hbar k_F/m$ を使った。Drude の式である。

## バーテックス補正: 輸送の緩和時間 {#vertex-correction}

積の平均 (C4) を、不純物の線で二本の $G$ を結ぶはしご型のダイヤグラムまで取る（Born 近似の自己エネルギー (C6) と組になる近似）。
電流のバーテックス $\Gamma_x({\bf k})$ は

$$
\begin{aligned}
\Gamma_x({\bf k})={}&v_x({\bf k})+n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\\
&\times\langle G^R\rangle({\bf k}',E_F)\langle G^A\rangle({\bf k}',E_F)\,\Gamma_x({\bf k}')
\end{aligned}
\tag{C10}
$$

を満たす。$\Gamma_x({\bf k})=\gamma\,v_x({\bf k})$ と置き、(C8) を使うと、${\bf k}'$ の向きについての平均で $v_x({\bf k}')$ が $v_x({\bf k})\cos\theta_{\bf kk'}$ になるので、
$\gamma=1+\gamma\,(1-\tau/\tau_{\rm tr})$、すなわち $\gamma=\tau_{\rm tr}/\tau$ を得る。ここで

$$
\frac{1}{\tau_{\rm tr}}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,
\big(1-\cos\theta_{\bf kk'}\big)\,\delta(E_F-\varepsilon_{\bf k'})
\tag{C11}
$$

は輸送の緩和時間である。式 (C9) の片方の $v_x$ を $\Gamma_x$ に替えると

$$
\sigma_{xx}=\frac{ne^2\tau_{\rm tr}}{m}
\tag{C12}
$$

となる。前方への散乱（$\cos\theta\simeq 1$）は電流をほとんど減らさないので、抵抗に効くのは $1-\cos\theta$ の重みを付けた散乱である。
寿命 $\tau$（一本の $G$ の平均で決まる）と、抵抗を決める $\tau_{\rm tr}$（積の平均で決まる）は一般に違う。
$u$ が点状（$u({\bf q})$ が ${\bf q}$ によらない）なら $\cos\theta$ の平均は 0 で、バーテックス補正は消えて $\tau_{\rm tr}=\tau$ になる。
この結果は Boltzmann 方程式で得られるものと同じである。

はしご型でなく、不純物の線が交差するダイヤグラムを足すと、弱局在の補正が出る（ここでは扱わない）。
