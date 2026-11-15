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

応答を観測する演算子 $\hat A$ は、それ自身が $F$ を含んでいてもよい（後で使う電流の演算子がそうである）。
期待値 $\langle\hat A\rangle(t)={\rm Tr}[\rho(t)\hat A(t)]$ を $F$ について一次まで取ると、$\hat A$ 自身の変化 $(\partial\hat A/\partial F)F$ を平衡の密度行列 $\rho_0$ で平均した項と、密度行列の変化 $\delta\rho$ による項の和になる。
後者は、相互作用表示で $\delta\rho_{\rm I}(t)=(i/\hbar)\int_{-\infty}^{t}dt'\,[\hat B_{\rm I}(t'),\rho_0]\,F(t')$ と一次まで展開して得られる。結果は

$$
\begin{aligned}
&\delta\langle\hat A\rangle(t)=\Big\langle\frac{\partial\hat A}{\partial F}\Big\rangle F(t)+\int dt'\,\chi^R_{AB}(t-t')\,F(t'),\\
&\chi^R_{AB}(t)=\frac{i}{\hbar}\,\theta(t)\,\big\langle\big[\big(\hat A\big|_{F=0}\big)_{\rm H}(t),\hat B\big]\big\rangle
\end{aligned}
\tag{C3}
$$

である。$\hat A_{\rm H}(t)=e^{i\hat H_0t/\hbar}\hat A\,e^{-i\hat H_0t/\hbar}$、$\langle\cdots\rangle$ は $\hat H_0$ の熱平衡での期待値である。第 2 項の交換子には、$F$ について 0 次の $\hat A\big|_{F=0}$ が入る。
$\hat A$ が $F$ を含まなければ第 1 項は無い。
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
電子の電荷を $-e$ とすると、ハミルトニアンの運動エネルギーは $\sum_i(\hat{\bf p}_i+e{\bf A})^2/2m$ で、${\bf A}$ について展開すると

$$
\begin{aligned}
&\hat H(t)=\hat H_0-\hat J^0_x\,A_x(t)+\frac{e^2}{2m}\hat N A_x(t)^2,\qquad \hat J^0_x=-e\sum_i\frac{\hat p_{x,i}}{m},\\
&\hat J_x\equiv-\frac{\partial\hat H}{\partial A_x}=\hat J^0_x-\frac{e^2}{m}\hat N A_x
\end{aligned}
\tag{C5}
$$

である（$\hat N$ は電子の数の演算子、和は電子についてとる）。
$\hat J_x$ が全電流の演算子（電荷 $-e$ と速度 $(\hat{\bf p}_i+e{\bf A})/m$ の積の和）で、$\hat J^0_x$ は ${\bf A}=0$ での $\hat J_x$ である。

式 (C2)、(C3) と見比べると、$\hat A=\hat J_x$、$\hat B=\hat J^0_x$、$F=A_x$ である。
これを式 (C3) に代入すると

$$
\begin{aligned}
&\delta\langle\hat J_x\rangle(t)=\Big\langle\frac{\partial\hat J_x}{\partial A_x}\Big\rangle A_x(t)+\int dt'\,\chi^R_{JJ}(t-t')\,A_x(t'),\\
&\chi^R_{JJ}(t)=\frac{i}{\hbar}\,\theta(t)\,\big\langle\big[\big(\hat J_x\big|_{A_x=0}\big)_{\rm H}(t),\hat J^0_x\big]\big\rangle
=\frac{i}{\hbar}\,\theta(t)\,\big\langle[\hat J^0_{x,{\rm H}}(t),\hat J^0_x]\big\rangle
\end{aligned}
\tag{C6}
$$

である（最後の等号で式 (C5) の $\hat J_x\big|_{A_x=0}=\hat J^0_x$ を使った）。
式 (C5) から $\partial\hat J_x/\partial A_x=-(e^2/m)\hat N$ で、その期待値は $-(e^2/m)\langle\hat N\rangle=-(e^2/m)n\Omega$ である（系の体積を $\Omega$、電子の数密度を $n$ とする）。したがって

$$
\delta\langle\hat J_x\rangle(t)=-\frac{e^2}{m}n\Omega\,A_x(t)+\int dt'\,\chi^R_{JJ}(t-t')\,A_x(t')
\tag{C7}
$$

となる。第 1 項が反磁性の部分、第 2 項が常磁性の部分で、どちらも線形応答の式 (C3) から出る。

$e^{-i\omega t}$ の成分では、第 2 項のたたみこみは積 $\chi^R_{JJ}(\omega)A_x(\omega)$ になる。
電流密度は $\delta\langle\hat J_x\rangle/\Omega$ で、$A_x(\omega)=E_x(\omega)/i\omega$ なので

$$
\sigma(\omega)=\frac{\delta\langle\hat J_x\rangle(\omega)}{\Omega\,E_x(\omega)}=\frac{i}{\omega}\Big[\frac{ne^2}{m}-\frac{1}{\Omega}\,\chi^R_{JJ}(\omega)\Big]
\tag{C8}
$$

を得る。これが久保公式の結果である。散逸（ジュール熱）を決めるのは $\sigma$ の実部である。反磁性の項は、電場と電流の位相が直交しており、散逸を与えない（$\sigma$ の虚部にだけ効く。虚数 $i$ を掛けることは、位相が 90 度遅れることである。これはインダクタンス効果である。ただし磁場によるものではなく、電子の慣性による運動インダクタンスで、エネルギーは電子の運動エネルギーとして蓄えられる）。

### 直流の極限: Kubo–Greenwood の式 {#kubo-greenwood}

式 (C8) の実部を取ると、反磁性の項は落ちて

$$
\begin{aligned}
{\rm Re}\,\sigma(\omega)&=\frac{1}{\omega\Omega}\,{\rm Im}\,\chi^R_{JJ}(\omega)\\
&=\frac{1}{\omega\Omega}\,{\rm Im}\sum_{nm}(f_n-f_m)\,\frac{|\langle n|\hat J^0_x|m\rangle|^2}{E_m-E_n-\hbar\omega-i\delta}\\
&=\frac{\pi}{\omega\Omega}\sum_{nm}(f_n-f_m)\,|\langle n|\hat J^0_x|m\rangle|^2\,\delta(E_m-E_n-\hbar\omega)
\end{aligned}
\tag{C9}
$$

となる。2 行目では式 (C4) を $\hat A=\hat B=\hat J^0_x$ として入れ、3 行目では ${\rm Im}\,[1/(x-i\delta)]=\pi\delta(x)$ を使った。
$\omega\to 0$ では、デルタ関数により $E_m=E_n+\hbar\omega$ なので $f_n-f_m\simeq -f'(E_n)\,\hbar\omega$ で、$\omega$ が消える。
低温では $-f'(E)=\delta(E-E_F)$ である。$\hat J^0_x$ の一体の行列要素は $-e\langle n|\hat v_x|m\rangle$ で、スピンの因子を考えて 2 倍することにすると（$n$、$m$ は軌道の状態）、

$$
\sigma_{xx}=\frac{2\pi\hbar e^2}{\Omega}\sum_{nm}|\langle n|\hat v_x|m\rangle|^2\,\delta(E_F-E_n)\,\delta(E_F-E_m)
\tag{C10}
$$

となる。
これは、フェルミ準位 $E_F$ での一体の量だけで決まる。
$\delta(E-H)=-{\rm Im}\,G^R(E)/\pi$ で書き直すと、Kubo–Greenwood の式

$$
\sigma_{xx}=\frac{2\hbar e^2}{\pi\Omega}\,
{\rm Tr}\big[\hat v_x\,{\rm Im}\,G^R(E_F)\,\hat v_x\,{\rm Im}\,G^R(E_F)\big]
\tag{C11}
$$

を得る。$\hat v_x$ は速度の演算子、${\rm Tr}$ は一体の（軌道の）状態についての和、因子 2 はスピンである。

一方、一体の系では、リングダイアグラム（$G$ 二つの積で、端点に $J$ を置いて閉じたもの）が厳密な結果を与える。
その線形応答関数は、時間順序のグリーン関数 $G^T(E)=\sum_m|m\rangle\langle m|/\big(E-E_m+i\delta\,{\rm sgn}(E_m-\mu)\big)$ を使って

$$
\chi^T_{JJ}(\omega)=\frac{i}{2\pi}\int dE\,{\rm Tr}\big[J\,G^T(E)\,J\,G^T(E+\hbar\omega)\big]
\tag{C12}
$$

と書ける（絶対零度。$J=-e\hat v_x$ は一体の演算子）。$E$ の積分で留数を取ると、占有の $n$ と非占有の $m$ の組だけが残る。極の位置は式 (C4) と同じで、違うのは $\omega\lt 0$ の側の $i\delta$ の向きだけである。
この $\chi^T$ を遅延型にするには、$\omega$ の極をすべて下半面に移せばよい（$\omega\to\omega+i\delta$。遅延型の応答関数は上半面で解析的である）。
占有を $f(E)$ で書けば、任意の温度で次のようになる:

$$
\begin{aligned}
\chi^R_{JJ}(\omega)&=\frac{1}{2\pi i}\int dE\,f(E)\\
&\quad\times{\rm Tr}\big[J\,\big(G^R(E)-G^A(E)\big)\,J\,G^R(E+\hbar\omega)\\
&\qquad+J\,G^A(E-\hbar\omega)\,J\,\big(G^R(E)-G^A(E)\big)\big]
\end{aligned}
\tag{C13}
$$

となる。$G^R={\rm Re}\,G+i\,{\rm Im}\,G$、$G^A={\rm Re}\,G-i\,{\rm Im}\,G$（${\rm Re}\,G$、${\rm Im}\,G$ はどちらもエルミート）と分けると、虚部は ${\rm Re}\,G\cdot{\rm Im}\,G$ と ${\rm Im}\,G\cdot{\rm Im}\,G$ の組み合わせになる。
${\rm Re}\,G$ を含む項のトレースは実数なので虚部には効かず、${\rm Im}\,G$ どうしの積（極の寄与）だけが残る。第 2 項で $E\to E+\hbar\omega$ とずらすと

$$
\begin{aligned}
{\rm Re}\,\sigma(\omega)&=\frac{1}{\omega\Omega}\,{\rm Im}\,\chi^R_{JJ}(\omega)\\
&=\frac{1}{\pi\omega\Omega}\int dE\,\big[f(E)-f(E+\hbar\omega)\big]\\
&\quad\times{\rm Tr}\big[J\,{\rm Im}\,G^R(E)\,J\,{\rm Im}\,G^R(E+\hbar\omega)\big]
\end{aligned}
\tag{C14}
$$

となる。これは任意の $\omega$ と温度で成り立つ（$\omega\neq 0$ では光学伝導度を与える）。$\omega\to 0$、低温の極限では、式 (C10) と一致する。

式 (C10) のデルタ関数の積の和（積分）は分かりにくいと思うので、以下で説明しておく。
$i\delta$ を有限に残すと、式 (C10) のデルタ関数は幅 $\delta$ のローレンツ関数 $\delta_\delta(x)=(1/\pi)\,\delta/(x^2+\delta^2)$ になる。
$\sum_n\delta_\delta(E-E_n)$ は状態密度 $D(E)$（系全体の、エネルギーあたりの状態の数）である。
もし行列要素 $|\langle n|\hat v_x|m\rangle|^2$ が $E_F$ の近くで一様（平均 $\overline{|v_x|^2}$）なら

$$
\begin{aligned}
\sigma_{xx}&\simeq\frac{2\pi\hbar e^2}{\Omega}\,\overline{|v_x|^2}\,\Big[\sum_n\delta_\delta(E_F-E_n)\Big]^2\\
&=\frac{2\pi\hbar e^2}{\Omega}\,\overline{|v_x|^2}\,D(E_F)^2
\end{aligned}
\tag{C15}
$$

で、$\delta$ によらない。一方、対角成分だけの場合は $n=m$ の項だけが残り、$\int dx\,\delta_\delta(x)^2=1/(2\pi\delta)$ から

$$
\begin{aligned}
\sigma_{xx}&\simeq\frac{2\pi\hbar e^2}{\Omega}\sum_n|\langle n|\hat v_x|n\rangle|^2\,\delta_\delta(E_F-E_n)^2\\
&=\frac{\hbar e^2}{\Omega}\,\overline{v_x^2}\,\frac{D(E_F)}{\delta}
\end{aligned}
\tag{C16}
$$

で、$\delta\to 0$ で発散する（散乱が無ければ電流は減衰しない）。
不純物があると平均した $G$ の幅は $\hbar/2\tau$ になり（[平均した一本のグリーン関数](#averaged-g)）、$\delta=\hbar/2\tau$ とすると式 (C16) は $ne^2\tau/m$、すなわち [Drude の式](#drude) になる。
二つの場合は、完全系から $\sum_m|\langle n|\hat v_x|m\rangle|^2=\langle n|\hat v_x^2|n\rangle$ が有限であることでつながる。$n$ とつながる $m$ がエネルギーの幅 $W$ に広がっていれば $|\langle n|\hat v_x|m\rangle|^2\approx\overline{v_x^2}/(D W)$ で、式 (C15) は式 (C16) の $\delta$ を $W$ に替えた形になる。

### $G^RG^A$ の項 {#grga}

式 (C11) は、不純物の配置を一つ決めたときの式である。観測される伝導度は、これを不純物の配置について平均したもの $\langle\cdots\rangle_{\rm imp}$ である。
${\rm Im}\,G^R=(G^R-G^A)/2i$ と書いて展開すると、$G^RG^A$ の項と $G^RG^R$、$G^AG^A$ の項が出る。
不純物が弱いとき主になるのは $G^RG^A$ の項で、

$$
\sigma_{xx}\simeq\frac{\hbar e^2}{\pi\Omega}\,{\rm Re}\,
\Big\langle {\rm Tr}\big[\hat v_x\,G^R(E_F)\,\hat v_x\,G^A(E_F)\big]\Big\rangle_{\rm imp}
\tag{C17}
$$

である（$G^RG^R$ の項は、$1/\tau$ について高次の補正を与えるだけである）。

**ミソ**: 平均を取るのは $G^R$ と $G^A$ の **積** である。

$$
\big\langle \hat v_x G^R \hat v_x G^A\big\rangle_{\rm imp}\neq
\hat v_x\langle G^R\rangle_{\rm imp}\,\hat v_x\langle G^A\rangle_{\rm imp}.
\tag{C18}
$$

$G^R$ と $G^A$ が同じ不純物で散乱される相関（ダイヤグラムでは、二本の線を結ぶ不純物の線）が残る。これがバーテックス補正である。

## 平均した一本のグリーン関数: 寿命 {#averaged-g}

一本の $G^R$ を平均すると並進対称性が戻り、波数 ${\bf k}$ で対角になる。

$$
\langle G^R\rangle_{\rm imp}({\bf k},E)=\frac{1}{E-\varepsilon_{\bf k}-\Sigma^R({\bf k},E)}.
\tag{C19}
$$

不純物について二次（Born 近似）の自己エネルギーは

$$
\Sigma^R({\bf k},E)=n_i\int\frac{d^3k'}{(2\pi)^3}\,
\frac{|u({\bf k}-{\bf k}')|^2}{E-\varepsilon_{\bf k'}+i\delta}
\tag{C20}
$$

で、その虚部が準粒子の寿命 $\tau$ を与える（実部はエネルギーのずれで、$E_F$ に繰り込む）:

$$
\begin{aligned}
&{\rm Im}\,\Sigma^R=-\frac{\hbar}{2\tau},\\
&\frac{1}{\tau}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,\delta(E-\varepsilon_{\bf k'}).
\end{aligned}
\tag{C21}
$$

$1/\tau$ はフェルミの黄金則による散乱の速さそのものである（[Basic I](./basic#td-vs-ti-perturbation) の、黄金則と自己エネルギーの虚部の対応）。

## バーテックス補正を入れない場合: Drude の式 {#drude}

式 (C17) の平均を、一本ずつの平均 (C19) の積で置きかえると（バーテックス補正を落とすと）、${\bf k}$ についての和になる。
$\tau$ が長い（$\hbar/\tau\ll E_F$）とき

$$
\begin{aligned}
\langle G^R\rangle({\bf k},E)\,\langle G^A\rangle({\bf k},E)
&=\frac{1}{(E-\varepsilon_{\bf k})^2+(\hbar/2\tau)^2}\\
&\simeq\frac{2\pi\tau}{\hbar}\,\delta(E-\varepsilon_{\bf k})
\end{aligned}
\tag{C22}
$$

なので、$v_x^2$ の角度平均 $v_F^2/3$ と、スピンあたりの状態密度 $N_1(E_F)$（単位体積あたり）を使って

$$
\sigma_{xx}=\frac{\hbar e^2}{\pi}\cdot\frac{2\pi\tau}{\hbar}\cdot\frac{v_F^2}{3}\,N_1(E_F)
=\frac{2e^2\tau v_F^2N_1(E_F)}{3}=\frac{ne^2\tau}{m}
\tag{C23}
$$

となる。最後は自由電子の $N_1(E_F)=mk_F/(2\pi^2\hbar^2)$、$n=k_F^3/(3\pi^2)$、$v_F=\hbar k_F/m$ を使った。Drude の式である。

## バーテックス補正: 輸送の緩和時間 {#vertex-correction}

積の平均 (C18) を、不純物の線で二本の $G$ を結ぶはしご型のダイヤグラムまで取る（Born 近似の自己エネルギー (C20) と組になる近似）。
電流のバーテックス $\Gamma_x({\bf k})$ は

$$
\begin{aligned}
\Gamma_x({\bf k})={}&v_x({\bf k})+n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\\
&\times\langle G^R\rangle({\bf k}',E_F)\langle G^A\rangle({\bf k}',E_F)\,\Gamma_x({\bf k}')
\end{aligned}
\tag{C24}
$$

を満たす。$\Gamma_x({\bf k})=\gamma\,v_x({\bf k})$ と置き、(C22) を使うと、${\bf k}'$ の向きについての平均で $v_x({\bf k}')$ が $v_x({\bf k})\cos\theta_{\bf kk'}$ になるので、
$\gamma=1+\gamma\,(1-\tau/\tau_{\rm tr})$、すなわち $\gamma=\tau_{\rm tr}/\tau$ を得る。ここで

$$
\frac{1}{\tau_{\rm tr}}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,
\big(1-\cos\theta_{\bf kk'}\big)\,\delta(E_F-\varepsilon_{\bf k'})
\tag{C25}
$$

は輸送の緩和時間である。式 (C23) の片方の $v_x$ を $\Gamma_x$ に替えると

$$
\sigma_{xx}=\frac{ne^2\tau_{\rm tr}}{m}
\tag{C26}
$$

となる。前方への散乱（$\cos\theta\simeq 1$）は電流をほとんど減らさないので、抵抗に効くのは $1-\cos\theta$ の重みを付けた散乱である。
寿命 $\tau$（一本の $G$ の平均で決まる）と、抵抗を決める $\tau_{\rm tr}$（積の平均で決まる）は一般に違う。
$u$ が点状（$u({\bf q})$ が ${\bf q}$ によらない）なら $\cos\theta$ の平均は 0 で、バーテックス補正は消えて $\tau_{\rm tr}=\tau$ になる。
この結果は Boltzmann 方程式で得られるものと同じである。

はしご型でなく、不純物の線が交差するダイヤグラムを足すと、弱局在の補正が出る（ここでは扱わない）。
