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

### 一般の線形応答 {#linear-response}

時間に依存する外場 $F(t)$ が、演算子 $\hat B$ を通して系に加わるとする:

$$
\hat H(t)=\hat H_0-\hat B\,F(t).
\tag{C2}
$$

相互作用表示で状態（または密度行列）を $F$ について一次まで展開すると、別の演算子 $\hat A$ の期待値の変化は

$$
\begin{aligned}
&\delta\langle\hat A\rangle(t)=\int_{-\infty}^{\infty}dt'\,\chi^R_{AB}(t-t')\,F(t'),\\
&\chi^R_{AB}(t)=\frac{i}{\hbar}\,\theta(t)\,\big\langle[\hat A_{\rm H}(t),\hat B]\big\rangle
\end{aligned}
\tag{C3}
$$

となる。$\hat A_{\rm H}(t)=e^{i\hat H_0t/\hbar}\hat A\,e^{-i\hat H_0t/\hbar}$、$\langle\cdots\rangle$ は $\hat H_0$ の熱平衡での期待値である。
$\theta(t)$ があるので $\chi^R$ は遅延型で、応答は原因 $F$ より後にしか起きない（[付録](./appendix#causality)）。

一体の系では、$\hat H_0$ の一体の固有状態 $|n\rangle$（エネルギー $E_n$、占有 $f_n$）を使い、$\hat A$、$\hat B$ を一体の演算子とすると、フーリエ変換は

$$
\chi^R_{AB}(\omega)=\int dt\,e^{i\omega t}\chi^R_{AB}(t)
=\sum_{nm}(f_n-f_m)\,\frac{\langle n|\hat A|m\rangle\langle m|\hat B|n\rangle}{E_m-E_n-\hbar\omega-i\delta}
\tag{C4}
$$

である（$\langle[\hat c^\dagger_n\hat c_m,\hat c^\dagger_m\hat c_n]\rangle=f_n-f_m$ を使った）。

### 電場への応答: 久保公式 {#response-to-field}

一様な電場は、ベクトルポテンシャル $A_x(t)$ を使って $E_x=-\partial A_x/\partial t$ として入れる（周期境界条件の下でも位置の演算子を使わずにすむ）。
電子の電荷を $-e$ とし、運動エネルギー $(\hat{\bf p}+e{\bf A})^2/2m$ を ${\bf A}$ について一次まで展開すると、摂動は

$$
\hat H'(t)=-\hat J_x\,A_x(t),\qquad \hat J_x=-e\sum_i\hat v_{x,i}
\tag{C5}
$$

である。$\hat J_x$ は全電流の演算子（和は電子についてとる）で、式 (C2) の $\hat B=\hat J_x$、$F=A_x$ に当たる。
電流密度は、常磁性の部分（$\hat J_x$ の期待値の変化）と、${\bf A}$ に比例する反磁性の部分の和である（$n$ は電子の数密度）:

$$
j_x(t)=\frac{1}{\Omega}\,\delta\langle\hat J_x\rangle(t)-\frac{ne^2}{m}A_x(t).
\tag{C6}
$$

$e^{-i\omega t}$ の成分では $A_x(\omega)=E_x(\omega)/i\omega$ なので、式 (C3) を使って

$$
\sigma(\omega)=\frac{j_x(\omega)}{E_x(\omega)}=\frac{i}{\omega}\Big[\frac{ne^2}{m}-\frac{1}{\Omega}\,\chi^R_{JJ}(\omega)\Big]
\tag{C7}
$$

を得る。これが久保公式である。散逸（ジュール熱）を決める実部は ${\rm Re}\,\sigma(\omega)={\rm Im}\,\chi^R_{JJ}(\omega)/(\omega\Omega)$ で、反磁性の項は虚部にだけ効く。

### 直流の極限: Kubo–Greenwood の式 {#kubo-greenwood}

式 (C4) を入れ、${\rm Im}\,[1/(x-i\delta)]=\pi\delta(x)$ を使うと

$$
{\rm Re}\,\sigma(\omega)=\frac{\pi}{\omega\Omega}\sum_{nm}(f_n-f_m)\,|\langle n|\hat J_x|m\rangle|^2\,\delta(E_m-E_n-\hbar\omega)
\tag{C8}
$$

となる。$\omega\to 0$ では、デルタ関数により $E_m=E_n+\hbar\omega$ なので $f_n-f_m\simeq -f'(E_n)\,\hbar\omega$ で、$\omega$ が消える。
低温では $-f'(E)=\delta(E-E_F)$ である。$\hat J_x$ の一体の行列要素は $-e\langle n|\hat v_x|m\rangle$ で、スピンの和で因子 2 が付くので（$n$、$m$ は軌道の状態）、
$\sigma_{xx}=(2\pi\hbar e^2/\Omega)\sum_{nm}|\langle n|\hat v_x|m\rangle|^2\delta(E_F-E_n)\delta(E_F-E_m)$ となる。
これは、フェルミ準位 $E_F$ での一体の量だけで決まる。
$\delta(E-H)=-{\rm Im}\,G^R(E)/\pi$ で書き直し、不純物の配置について平均すると、Kubo–Greenwood の式

$$
\sigma_{xx}=\frac{2\hbar e^2}{\pi\Omega}\,
\Big\langle {\rm Tr}\big[\hat v_x\,{\rm Im}\,G^R(E_F)\,\hat v_x\,{\rm Im}\,G^R(E_F)\big]\Big\rangle_{\rm imp}.
\tag{C9}
$$

を得る。$\hat v_x$ は速度の演算子、${\rm Tr}$ は一体の（軌道の）状態についての和、因子 2 はスピンである。

### $G^RG^A$ の項 {#grga}


${\rm Im}\,G^R=(G^R-G^A)/2i$ と書いて展開すると、$G^RG^A$ の項と $G^RG^R$、$G^AG^A$ の項が出る。
不純物が弱いとき主になるのは $G^RG^A$ の項で、

$$
\sigma_{xx}\simeq\frac{\hbar e^2}{\pi\Omega}\,{\rm Re}\,
\Big\langle {\rm Tr}\big[\hat v_x\,G^R(E_F)\,\hat v_x\,G^A(E_F)\big]\Big\rangle_{\rm imp}
\tag{C10}
$$

である（$G^RG^R$ の項は、$1/\tau$ について高次の補正を与えるだけである）。

**ミソ**: 平均を取るのは $G^R$ と $G^A$ の **積** である。

$$
\big\langle \hat v_x G^R \hat v_x G^A\big\rangle_{\rm imp}\neq
\hat v_x\langle G^R\rangle_{\rm imp}\,\hat v_x\langle G^A\rangle_{\rm imp}.
\tag{C11}
$$

$G^R$ と $G^A$ が同じ不純物で散乱される相関（ダイヤグラムでは、二本の線を結ぶ不純物の線）が残る。これがバーテックス補正である。

## 平均した一本のグリーン関数: 寿命 {#averaged-g}

一本の $G^R$ を平均すると並進対称性が戻り、波数 ${\bf k}$ で対角になる。

$$
\langle G^R\rangle_{\rm imp}({\bf k},E)=\frac{1}{E-\varepsilon_{\bf k}-\Sigma^R({\bf k},E)}.
\tag{C12}
$$

不純物について二次（Born 近似）の自己エネルギーは

$$
\Sigma^R({\bf k},E)=n_i\int\frac{d^3k'}{(2\pi)^3}\,
\frac{|u({\bf k}-{\bf k}')|^2}{E-\varepsilon_{\bf k'}+i\delta}
\tag{C13}
$$

で、その虚部が準粒子の寿命 $\tau$ を与える（実部はエネルギーのずれで、$E_F$ に繰り込む）:

$$
\begin{aligned}
&{\rm Im}\,\Sigma^R=-\frac{\hbar}{2\tau},\\
&\frac{1}{\tau}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,\delta(E-\varepsilon_{\bf k'}).
\end{aligned}
\tag{C14}
$$

$1/\tau$ はフェルミの黄金則による散乱の速さそのものである（[Basic I](./basic#td-vs-ti-perturbation) の、黄金則と自己エネルギーの虚部の対応）。

## バーテックス補正を入れない場合: Drude の式 {#drude}

式 (C10) の平均を、一本ずつの平均 (C12) の積で置きかえると（バーテックス補正を落とすと）、${\bf k}$ についての和になる。
$\tau$ が長い（$\hbar/\tau\ll E_F$）とき

$$
\begin{aligned}
\langle G^R\rangle({\bf k},E)\,\langle G^A\rangle({\bf k},E)
&=\frac{1}{(E-\varepsilon_{\bf k})^2+(\hbar/2\tau)^2}\\
&\simeq\frac{2\pi\tau}{\hbar}\,\delta(E-\varepsilon_{\bf k})
\end{aligned}
\tag{C15}
$$

なので、$v_x^2$ の角度平均 $v_F^2/3$ と、スピンあたりの状態密度 $N_1(E_F)$（単位体積あたり）を使って

$$
\sigma_{xx}=\frac{\hbar e^2}{\pi}\cdot\frac{2\pi\tau}{\hbar}\cdot\frac{v_F^2}{3}\,N_1(E_F)
=\frac{2e^2\tau v_F^2N_1(E_F)}{3}=\frac{ne^2\tau}{m}
\tag{C16}
$$

となる。最後は自由電子の $N_1(E_F)=mk_F/(2\pi^2\hbar^2)$、$n=k_F^3/(3\pi^2)$、$v_F=\hbar k_F/m$ を使った。Drude の式である。

## バーテックス補正: 輸送の緩和時間 {#vertex-correction}

積の平均 (C11) を、不純物の線で二本の $G$ を結ぶはしご型のダイヤグラムまで取る（Born 近似の自己エネルギー (C13) と組になる近似）。
電流のバーテックス $\Gamma_x({\bf k})$ は

$$
\begin{aligned}
\Gamma_x({\bf k})={}&v_x({\bf k})+n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\\
&\times\langle G^R\rangle({\bf k}',E_F)\langle G^A\rangle({\bf k}',E_F)\,\Gamma_x({\bf k}')
\end{aligned}
\tag{C17}
$$

を満たす。$\Gamma_x({\bf k})=\gamma\,v_x({\bf k})$ と置き、(C15) を使うと、${\bf k}'$ の向きについての平均で $v_x({\bf k}')$ が $v_x({\bf k})\cos\theta_{\bf kk'}$ になるので、
$\gamma=1+\gamma\,(1-\tau/\tau_{\rm tr})$、すなわち $\gamma=\tau_{\rm tr}/\tau$ を得る。ここで

$$
\frac{1}{\tau_{\rm tr}}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,
\big(1-\cos\theta_{\bf kk'}\big)\,\delta(E_F-\varepsilon_{\bf k'})
\tag{C18}
$$

は輸送の緩和時間である。式 (C16) の片方の $v_x$ を $\Gamma_x$ に替えると

$$
\sigma_{xx}=\frac{ne^2\tau_{\rm tr}}{m}
\tag{C19}
$$

となる。前方への散乱（$\cos\theta\simeq 1$）は電流をほとんど減らさないので、抵抗に効くのは $1-\cos\theta$ の重みを付けた散乱である。
寿命 $\tau$（一本の $G$ の平均で決まる）と、抵抗を決める $\tau_{\rm tr}$（積の平均で決まる）は一般に違う。
$u$ が点状（$u({\bf q})$ が ${\bf q}$ によらない）なら $\cos\theta$ の平均は 0 で、バーテックス補正は消えて $\tau_{\rm tr}=\tau$ になる。
この結果は Boltzmann 方程式で得られるものと同じである。

はしご型でなく、不純物の線が交差するダイヤグラムを足すと、弱局在の補正が出る（ここでは扱わない）。
