# 電気伝導度計算の基本 {#conductivity}

::: info 下書き
**一体問題** として、線形応答（久保公式）から不純物による電気抵抗を導く。
ミソは、不純物の配置について平均を取るのが **グリーン関数の積** だという点である（一本ずつ平均してから掛けるのでは足りない）。
グリーン関数の符号の約束は [Basic I](./basic) と同じで、遅延グリーン関数は $G^R(E)=1/(E-H+i\delta)$。
遅延型の応答関数については [付録](./appendix#causality) も見よ。
:::

## 設定 {#setting}

ここでは、真空の中に原子（不純物）がランダムに並んでいる系を考える。
体積 $\Omega$ の中に $N_i$ 個の同じ原子があり、その位置 ${\bf R}_j$ は互いに独立に、一様に分布する。
相互作用の無い電子（質量 $m$、$\varepsilon_{\bf k}=\hbar^2k^2/2m$）が、これらの原子のポテンシャル

$$
V({\bf r})=\sum_{j} u({\bf r}-{\bf R}_j)
\tag{C1}
$$

を感じるとする（$u$ は原子一つのポテンシャル）。不純物の数密度を $n_i=N_i/\Omega$、$u$ のフーリエ変換を $u({\bf q})=\int d^3r\, e^{-i{\bf q}\cdot{\bf r}}u({\bf r})$ と書く。
観測される量は、不純物の配置 $\{ {\bf R}_j \}$ についての平均 $\langle\cdots\rangle_{\rm imp}$ である。

ここでの議論は、固体中の不純物の問題にそのまま直せる。平面波 $|{\bf k}\rangle$ と $\varepsilon_{\bf k}$ を母体の結晶のブロッホ状態とバンドのエネルギーに、$u$ を不純物が置き換わったことによるポテンシャルの変化に替えればよい。

**ここで扱うのは一体問題である。** 電子どうしの相互作用は無く、各電子は不純物の作る静的なポテンシャル (C1) の中を独立に動く。
ハミルトニアンは一体の $H=\hat{\bf p}^2/2m+V({\bf r})$ で、以下の $G^R(E)=1/(E-H+i\delta)$ はその一体のグリーン関数である。
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

> 式 (C3) などに交換子が出る理由に注意しておく。これは $\hat A_{\rm H}(t)$ を $t$ で微分すると $d\hat A_{\rm H}/dt=(i/\hbar)[\hat H_0,\hat A]_{\rm H}$ と交換子が出ることから理解できる。
> 時間発展はヒルベルト空間のユニタリーな回転で、微小な回転 $e^{i\epsilon\hat G}$ による演算子の変化は $e^{i\epsilon\hat G}\hat A\,e^{-i\epsilon\hat G}\simeq\hat A+i\epsilon[\hat G,\hat A]$ である。
> 摂動 $-\hat B F$ が加わると、時間発展に $\hat B$ による微小な回転が加わるので、$\hat A$ の変化は $\hat B$ との交換子になる。
> これは、$3\times 3$ のエルミート（実対称）行列 $M$ が微小な空間回転 $R=e^{i\epsilon L}\simeq 1+i\epsilon L$（$L$ はエルミート）で $RMR^{\dagger}\simeq M+i\epsilon[L,M]$ と変換を受けるのと同じである。

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

となる。$G^R={\rm Re}\,G+i\,{\rm Im}\,G$、$G^A={\rm Re}\,G-i\,{\rm Im}\,G$（${\rm Re}\,G$、${\rm Im}\,G$ はどちらもエルミート。${\rm Im}\,G=-\pi\delta(E-H)$ は極だけの部分である）と分けると、虚部は ${\rm Re}\,G\cdot{\rm Im}\,G$ と ${\rm Im}\,G\cdot{\rm Im}\,G$ の組み合わせになる。
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

すこし考察してみる。もし行列要素 $|\langle n|\hat v_x|m\rangle|^2$ が $E_F$ の近くで一様（平均 $\overline{|v_x|^2}$）なら

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

で、$\delta\to 0$ で発散する（散乱が無ければ電流は減衰しない）。大雑把には、$\delta$ はフォノンなどによる準粒子の寿命 $\tau$ から来る幅 $\hbar/2\tau$ と考えてもよい。
実際には、行列要素は対角成分だけではないにせよ、$n$ とつながる $m$ の数は（完全系から $\sum_m|\langle n|\hat v_x|m\rangle|^2=\langle n|\hat v_x^2|n\rangle$ が有限なので）有限であり、式 (C16) が正しいと考えてよい。

式 (C16) の「フェルミ面の電子が電気伝導に効く」は、Drude の「すべての電子が $\tau$ の間加速される」と、次のようにつながる。
式 (C16) の右辺の次の因子は、一つのバンド、スピンあたりで、$v_x\,\delta(E_F-\varepsilon_{\bf k})=-\hbar^{-1}\partial_{k_x}\theta(E_F-\varepsilon_{\bf k})$（$\theta$ は占有の階段関数、$v_x=\hbar^{-1}\partial\varepsilon_{\bf k}/\partial k_x$）を使って $k_x$ について部分積分すると（境界の項はブリルアン帯の周期性で消える）

$$
\begin{aligned}
\frac{\overline{v_x^2}\,D(E_F)}{\Omega}&=\int\frac{d^3k}{(2\pi)^3}\,v_x^2\,\delta(E_F-\varepsilon_{\bf k})\\
&=\frac{1}{\hbar}\int\frac{d^3k}{(2\pi)^3}\,\theta(E_F-\varepsilon_{\bf k})\,\frac{\partial v_x}{\partial k_x}\\
&=\int_{\rm occ}\frac{d^3k}{(2\pi)^3}\,\frac{1}{\hbar^2}\frac{\partial^2\varepsilon_{\bf k}}{\partial k_x^2}
\end{aligned}
\tag{C17}
$$

となる。$v_x^2$ の一方の $v_x$ はデルタ関数と組になって $\theta$ の微分になり（$\theta$ を $k_x$ で微分すると $\partial\varepsilon_{\bf k}/\partial k_x=\hbar v_x$ が出る）、部分積分で微分がもう一方の $v_x$ に移って $\partial v_x/\partial k_x=\hbar^{-1}\partial^2\varepsilon_{\bf k}/\partial k_x^2$（逆有効質量）になる。
左辺はフェルミ面の量、右辺は占有された全状態の逆有効質量の和である。
自由電子では $\partial^2\varepsilon_{\bf k}/\partial k_x^2=\hbar^2/m$ が一定なので、右辺は $n/2m$ となる（$n$ は両スピンを含む電子の数密度で、$n/2$ がスピンあたり）。
したがって、式 (C16) で $\delta=\hbar/2\tau$ とし、式 (C17) を使うと

$$
\sigma_{xx}=2e^2\tau\,\frac{\overline{v_x^2}\,D(E_F)}{\Omega}=\frac{ne^2\tau}{m}
\tag{C18}
$$

となり、[Drude の式](#drude) が出る。
このフェルミ面での有効質量を一定とする近似は、Na、Cu、Al などではよく成り立つはずである。

> 占有された領域の積分とフェルミ面の上の積分が部分積分で結ばれるという式 (C17) の構造は、ベリー曲率を使った異常ホール伝導度の計算（占有状態でのベリー曲率の積分を、フェルミ面の上のベリー位相の積分に書き直す）と同じである。


### $G^RG^A$ の項 {#grga}

式 (C11) は、不純物の配置を一つ決めたときの式である。観測される伝導度は、これを不純物の配置について平均したもの $\langle\cdots\rangle_{\rm imp}$ である。
${\rm Im}\,G^R=(G^R-G^A)/2i$ と書いて展開すると、$G^RG^A$ の項と $G^RG^R$、$G^AG^A$ の項が出る。
不純物が弱いとき主になるのは $G^RG^A$ の項で、

$$
\sigma_{xx}\simeq\frac{\hbar e^2}{\pi\Omega}\,{\rm Re}\,
\Big\langle {\rm Tr}\big[\hat v_x\,G^R(E_F)\,\hat v_x\,G^A(E_F)\big]\Big\rangle_{\rm imp}
\tag{C19}
$$

である（$G^RG^R$ の項は、$1/\tau$ について高次の補正を与えるだけである）。

**ミソ**: 平均を取るのは $G^R$ と $G^A$ の **積** である。

$$
\big\langle \hat v_x G^R \hat v_x G^A\big\rangle_{\rm imp}\neq
\hat v_x\langle G^R\rangle_{\rm imp}\,\hat v_x\langle G^A\rangle_{\rm imp}.
\tag{C20}
$$

$G^R$ と $G^A$ が同じ不純物で散乱される相関（ダイヤグラムでは、二本の線を結ぶ不純物の線）が残る。これがバーテックス補正である。

もし式 (C20) で等号としてしまうと（一本ずつの平均の積にすると）、式 (C16) の形になる。このままでは、一体問題において不純物は抵抗を与えない（$\delta\to 0$ で発散する）という結論になってしまう。

## 不純物の配置についての平均 {#impurity-average}

式 (C20) の積の平均を計算するには、そもそもの一体問題に戻る。$G^R(E)=(E-H_0-V+i\delta)^{-1}$ は、不純物の配置を一つ決めた一体のハミルトニアン $H=H_0+V$（$H_0$ は運動エネルギー、$V$ は式 (C1)）のグリーン関数である。これを $V$ について展開し、積の平均を、$V$ の対の平均の積に分解する。
平面波 $|{\bf k}\rangle$ での行列要素は $\langle{\bf k}|V|{\bf k}'\rangle=\Omega^{-1}u({\bf k}-{\bf k}')\sum_j e^{-i({\bf k}-{\bf k}')\cdot{\bf R}_j}$ である。
不純物の位置は乱雑なので、位相 $e^{-i{\bf q}\cdot{\bf R}_j}$ は平均すると打ち消し合い（ランダム位相）、残るのは同じ不純物 $j$ について位相が消える組だけである:

$$
\big\langle \langle{\bf k}|V|{\bf k}'\rangle\langle{\bf k}'|V|{\bf k}\rangle\big\rangle_{\rm imp}=\frac{n_i}{\Omega}\,|u({\bf k}-{\bf k}')|^2
\tag{C21}
$$

（${\bf q}=0$ の成分は一様なエネルギーのずれで、$E_F$ に繰り込む。同じ不純物に三回以上かかる項は、不純物が弱いとして落とす。）
これを使って式 (C19) を展開する。平面波 $|{\bf k}\rangle$ では $\hat v_x$ は対角（$v_x({\bf k})=\hbar k_x/m$）なので、トレースは

$$
\begin{aligned}
{\rm Tr}\big[\hat v_x G^R\hat v_x G^A\big]
&=\sum_{ {\bf k}{\bf k}'}v_x({\bf k})\,G^R_{ {\bf k}{\bf k}'}\,v_x({\bf k}')\,G^A_{ {\bf k}'{\bf k}}
\end{aligned}
\tag{C22}
$$

である（$G_{ {\bf k}{\bf k}'}=\langle{\bf k}|G|{\bf k}'\rangle$。$E=E_F$ は書かない）。$G^R=(E-H_0-V+i\delta)^{-1}$ を $V$ について展開すると

$$
\begin{aligned}
G^R_{ {\bf k}{\bf k}'}={}&G^R_0({\bf k})\,\delta_{ {\bf k}{\bf k}'}+G^R_0({\bf k})\,V_{ {\bf k}{\bf k}'}\,G^R_0({\bf k}')\\
&+\sum_{ {\bf k}_1}G^R_0({\bf k})\,V_{ {\bf k}{\bf k}_1}\,G^R_0({\bf k}_1)\,V_{ {\bf k}_1{\bf k}'}\,G^R_0({\bf k}')+\cdots
\end{aligned}
\tag{C23}
$$

である。ここで $G^R_0({\bf k})=1/(E-\varepsilon_{\bf k}+i\delta)$、$\varepsilon_{\bf k}=\hbar^2k^2/2m$、$V_{ {\bf k}{\bf k}'}=\langle{\bf k}|V|{\bf k}'\rangle$ で、$G^A$ も同様（$G^A_0({\bf k})=1/(E-\varepsilon_{\bf k}-i\delta)$）である。
これを式 (C22) に入れて平均する。$V$ について 1 次の項は ${\bf q}=0$ の成分だけが残り、$E_F$ に繰り込まれている。2 次の項は、二つの $V$ がどこにあるかで三つに分かれ、式 (C21) を使うと

$$
\begin{aligned}
&\Big\langle{\rm Tr}\big[\hat v_x G^R\hat v_x G^A\big]\Big\rangle_{\rm imp}\\
&\quad=\sum_{\bf k}v_x({\bf k})^2\,G^R_0({\bf k})\,G^A_0({\bf k})\\
&\qquad+\sum_{\bf k}v_x({\bf k})^2\,G^R_0({\bf k})\,\Sigma^R_2({\bf k})\,G^R_0({\bf k})\,G^A_0({\bf k})\\
&\qquad+\sum_{\bf k}v_x({\bf k})^2\,G^R_0({\bf k})\,G^A_0({\bf k})\,\Sigma^A_2({\bf k})\,G^A_0({\bf k})\\
&\qquad+\sum_{ {\bf k}{\bf k}'}v_x({\bf k})\,G^R_0({\bf k})\,G^A_0({\bf k})\\
&\qquad\quad\times\frac{n_i}{\Omega}\,|u({\bf k}-{\bf k}')|^2\,G^R_0({\bf k}')\,G^A_0({\bf k}')\,v_x({\bf k}')\\
&\qquad+O(V^4)
\end{aligned}
\tag{C24}
$$

となる。ここで $\Sigma^{R}_2({\bf k})=\sum_{ {\bf k}_1}(n_i/\Omega)\,|u({\bf k}-{\bf k}_1)|^2\,G^R_0({\bf k}_1)$（$\Sigma^A_2$ は $G^R_0$ を $G^A_0$ にしたもの）である。
右辺の第 1 項は不純物の無い項、第 2 項と第 3 項は一本の $G$ の中の二つの $V$ を結んだ項（自己エネルギー $\Sigma_2$ が一つ挟まる）、第 4 項は $G^R$ の $V$ と $G^A$ の $V$ を結んだ項である。
一本の $G$ の中の平均 $\langle V_{ {\bf k}{\bf k}_1}V_{ {\bf k}_1{\bf k}'}\rangle_{\rm imp}$ は、ランダム位相のため ${\bf k}={\bf k}'$ のときだけ残る。平均した後の $G$ が ${\bf k}$ について対角（並進対称）になるのはこのためである。

$V$ について高次の項も同じで、$V$ の積の平均は対の平均の積の和になる。ダイヤグラムでは、組になった二つの $V$ を不純物の線で結ぶ。線には次の三通りがある。

- 一本の $G$ の中で結ぶ線: 自己エネルギーになる（次節）。
- $G^R$ と $G^A$ を結ぶ線: バーテックス補正になる（[バーテックス補正](#vertex-correction) の節）。
- 交差する線: 弱局在の補正になる（ここでは扱わない）。

式 (C20) で等号とすることは、$G^R$ と $G^A$ を結ぶ線（式 (C24) の第 4 項とその高次）を落とすことにあたる。

## 自己エネルギーを含むグリーン関数: 寿命 {#averaged-g}

ランダム位相の平均の結果、一本の $G$ の中で結ぶ線はすべて自己エネルギーにまとまる。つまり $G^R$ を平均する代わりに、次の $\bar G^R$（並進対称で、波数 ${\bf k}$ について対角）を考え、その積を取ればよい。

$$
\bar G^R({\bf k},E)=\frac{1}{E-\varepsilon_{\bf k}-\Sigma^R({\bf k},E)}.
\tag{C25}
$$

一本の $G$ の中で結ぶ線を、互いに交差しないもの（Born 近似）まで取ると、式 (C21) から自己エネルギーは

$$
\Sigma^R({\bf k},E)=n_i\int\frac{d^3k'}{(2\pi)^3}\,
\frac{|u({\bf k}-{\bf k}')|^2}{E-\varepsilon_{\bf k'}+i\delta}
\tag{C26}
$$

で（式 (C24) の $\Sigma^R_2$ の和を積分にしたもの）、その虚部が準粒子の寿命 $\tau$ を与える（実部はエネルギーのずれで、$E_F$ に繰り込む）:

$$
\begin{aligned}
&{\rm Im}\,\Sigma^R=-\frac{\hbar}{2\tau},\\
&\frac{1}{\tau}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,\delta(E-\varepsilon_{\bf k'}).
\end{aligned}
\tag{C27}
$$

$1/\tau$ はフェルミの黄金則による散乱の速さそのものである（[Basic I](./basic#td-vs-ti-perturbation) の、黄金則と自己エネルギーの虚部の対応）。

## バーテックス補正を入れない場合: Drude の式 {#drude}

式 (C20) で等号とした近似（$G^R$ と $G^A$ を結ぶ線を落とす近似）でも、$\bar G$ には自己エネルギーによる幅 $\hbar/2\tau$ が入っている。
$\tau$ が長い（$\hbar/\tau\ll E_F$）とき

$$
\begin{aligned}
\bar G^R({\bf k},E)\,\bar G^A({\bf k},E)
&=\frac{1}{(E-\varepsilon_{\bf k})^2+(\hbar/2\tau)^2}\\
&\simeq\frac{2\pi\tau}{\hbar}\,\delta(E-\varepsilon_{\bf k})
\end{aligned}
\tag{C28}
$$

で、これは式 (C16) のローレンツ関数の幅 $\delta$ を $\hbar/2\tau$ にしたものにあたる。したがって式 (C19) は

$$
\sigma_{xx}=2e^2\tau\,\frac{\overline{v_x^2}\,D(E_F)}{\Omega}
\tag{C29}
$$

となり、式 (C18) と同じ Drude の式（自由電子では $ne^2\tau/m$）が出る。
素朴に $\delta\to 0$ とすると抵抗が 0 になったが、不純物による散乱の幅が $\delta$ の役を担い、不純物は寿命 $\tau$ を通して抵抗を与える。

## バーテックス補正: 輸送の緩和時間 {#vertex-correction}

積の平均 (C20) を、不純物の線で二本の $G$ を結ぶはしご型のダイヤグラムまで取る（Born 近似の自己エネルギー (C26) と組になる近似）。
電流のバーテックス $\Gamma_x({\bf k})$ は

$$
\begin{aligned}
\Gamma_x({\bf k})={}&v_x({\bf k})+n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\\
&\times\bar G^R({\bf k}',E_F)\bar G^A({\bf k}',E_F)\,\Gamma_x({\bf k}')
\end{aligned}
\tag{C30}
$$

を満たす。$\Gamma_x({\bf k})=\gamma\,v_x({\bf k})$ と置き、(C28) を使うと、${\bf k}'$ の向きについての平均で $v_x({\bf k}')$ が $v_x({\bf k})\cos\theta_{\bf kk'}$ になるので、
$\gamma=1+\gamma\,(1-\tau/\tau_{\rm tr})$、すなわち $\gamma=\tau_{\rm tr}/\tau$ を得る。ここで

$$
\frac{1}{\tau_{\rm tr}}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,
\big(1-\cos\theta_{\bf kk'}\big)\,\delta(E_F-\varepsilon_{\bf k'})
\tag{C31}
$$

は輸送の緩和時間である。式 (C29) の片方の $v_x$ を $\Gamma_x$ に替えると

$$
\sigma_{xx}=2e^2\tau_{\rm tr}\,\frac{\overline{v_x^2}\,D(E_F)}{\Omega}
\tag{C32}
$$

となる（自由電子では $ne^2\tau_{\rm tr}/m$）。前方への散乱（$\cos\theta\simeq 1$）は電流をほとんど減らさないので、抵抗に効くのは $1-\cos\theta$ の重みを付けた散乱である。
寿命 $\tau$（一本の $G$ の平均で決まる）と、抵抗を決める $\tau_{\rm tr}$（積の平均で決まる）は一般に違う。
$u$ が点状（$u({\bf q})$ が ${\bf q}$ によらない）なら $\cos\theta$ の平均は 0 で、バーテックス補正は消えて $\tau_{\rm tr}=\tau$ になる。
この結果は Boltzmann 方程式で得られるものと同じである。

はしご型でなく、不純物の線が交差するダイヤグラムを足すと、弱局在の補正が出る（ここでは扱わない）。

## メモ: フェルミ面項とフェルミ海項 {#fermi-sea-surface}

> **予備: ブリルアン域の幾何学（トポロジー）**
>
> 二次元のトーラスで理解しておけば、直感的にそれなりに理解できる。これがトポロジカルな物理量を議論するときのキモである。いきなり計算に入ってしまうと、よく分からなくなる。
>
> ${\bf k}$ と ${\bf k}+{\bf G}$（${\bf G}$ は逆格子ベクトル）は同じ状態を表すので、ブリルアン域は向かい合う面を貼り合わせた三次元のトーラス $T^3$ である。
> $T^2$（二次元のトーラス）とは、ドーナツの表面のような二次元の空間のことである。$T^3$（三次元のトーラス）とは、$T^2$ の三次元版で、境界の無い閉じた三次元の空間の一種である。
> $T^2$ には、一点に縮めることのできない周回路（$T^1$）が 2 種類ある（連続変形で移り合う周回路は同じとみなす）。ドーナツの大きい円を「輪」、断面の小さい円を「管」と呼ぶと、輪に沿って回るものと、管を回るものである。ほかの周回路はこの 2 種類の組み合わせで、例えば、管を 2 周回ったあと輪に沿って 1 周回る経路などである。
> $T^3$ では、そのような周回路（$T^1$）が $k_x$、$k_y$、$k_z$ の向きに巻く 3 種類あり、さらに一点に縮めることのできない閉じた面（$T^2$）が 3 種類ある。これらは境界の無い閉じた面（二次元のトーラス）であり、$k_x=$一定 などで指定される。
> 例えば $k_z=$一定 の断面 $T^2(k_z)$ は、ブリルアン域を $k_x$、$k_y$ の向きに巻く。このような面の上の積分は、境界の上の積分に書き直すことができない。
> 一方、ブリルアン域の中に収まる領域（下のポケット $S$ のようなもの）は境界を持ち、その上の積分は Stokes の定理で境界の上の積分に書き直せる。

なぜそもそもこんなことを考えないといけないかと言うと、線形応答の計算において（トポロジカルな応答として）問題点が発生する場合があるからである。線形応答では $T^3$ の上の積分量を計算する。

単純化して考えると、$T^3$ の各点 $(x,y,z)$ に位相因子 $e^{i\theta(x,y,z)}$ が指定されている場合である。振幅 × 位相の形で与えられる場合もあるが、振幅の部分は本質的でなく、あとから考えればよい。

位相は連続的かつ滑らかに変化するものとする。ただし、特異線（渦糸。近づく向きによって位相の極限が違うので、そこで位相が定まらない。周りを一周すると位相が $2\pi$ の整数倍まわる。例えば $z$ 軸の上の $\theta=\arg(x+iy)$）や、バンドの縮退する特異点も、場合により考える。三次元では、位相の特異は普通は線になる。

久保公式の静的な伝導度 $\sigma_{ij}$ を、重み $-\partial f/\partial E$ を持つ **フェルミ面項** と、重み $f(E)$ で占有状態全体を積分する **フェルミ海項** に分ける書き方がある（Středa の分け方。異常ホール効果の文献でよく使われる）。
フェルミ海項は $i\leftrightarrow j$ について反対称なので、$\sigma_{xx}$ には効かない。式 (C19) はフェルミ面項のうち $G^RG^A$ の部分である。

ただし、この分け方は一つに決まらない。部分積分で一方を他方に書き直せるからである。
式 (C17) では、占有領域の積分（Drude の $n/m$）とフェルミ面での $\overline{v_x^2}D(E_F)$ が同じ量だった。
異常ホール伝導度でも、占有状態のベリー曲率の積分は、フェルミ面の上のベリー位相に書き直せる（F. D. M. Haldane, Phys. Rev. Lett. **93**, 206602 (2004)）。
したがって「フェルミ面の電子が効くか、全部の電子が効くか」という区別は、物理としては意味を持たない。

例外は整数の部分である。次の設定で式で見る。
三次元の結晶（簡単のため直方晶で、格子定数を $a$、$b$、$c$ とする）で、バンド 1 は完全に占有され（フェルミ面が無い）、バンド 2 はブリルアン域の中の領域 $S$（電子のポケット）だけが占有されているとする。フェルミ面は $S$ の境界の閉曲面 $\partial S$ である。
ブロッホ関数の周期部分 $|u_{n{\bf k}}\rangle$ のベリー接続とベリー曲率を

$$
{\bf a}_n({\bf k})=i\langle u_{n{\bf k}}|\nabla_{\bf k}u_{n{\bf k}}\rangle,\qquad
{\bf\Omega}_n({\bf k})=\nabla_{\bf k}\times{\bf a}_n({\bf k})
\tag{C33}
$$

とする。内因性の異常ホール伝導度は（符号は ${\bf\Omega}$ の向きの約束による）、$k_z$ ごとの断面に分けて

$$
\begin{aligned}
\sigma_{xy}&=\frac{e^2}{\hbar}\sum_n\int_{\rm occ}\frac{d^3k}{(2\pi)^3}\,\Omega_n^z\\
&=\frac{e^2}{h}\,\frac{1}{2\pi}\int dk_z\,\Big[\frac{1}{2\pi}\int_{T^2(k_z)}\Omega_1^z\,d^2k+\frac{1}{2\pi}\int_{S(k_z)}\Omega_2^z\,d^2k\Big]
\end{aligned}
\tag{C34}
$$

と書ける。$T^2(k_z)$ は $k_z$ を止めたブリルアン域の断面、$S(k_z)$ はその断面でのポケットである。

**第 1 項（バンド 1）**: 閉じた面 $T^2(k_z)$ の上の積分 $C_1^z=\frac{1}{2\pi}\int_{T^2(k_z)}\Omega_1^z\,d^2k$ はチャーン数で、整数である。バンド 1 にはギャップがあるので $k_z$ によらない。$T^2(k_z)$ は閉じた面なので、これを境界の上の積分には書き直せない。
$k_z$ の積分は $2\pi/c$ を与えるので、第 1 項は $(e^2/h)\,C_1^z/c$ である。ブリルアン域を他の向きに巻く面のチャーン数 $C_1^x$、$C_1^y$ も合わせると

$$
\sigma_{ij}^{(1)}=\frac{e^2}{2\pi h}\,\epsilon_{ijl}\,G_l,\qquad
{\bf G}=2\pi\Big(\frac{C_1^x}{a},\frac{C_1^y}{b},\frac{C_1^z}{c}\Big)
\tag{C35}
$$

で、${\bf G}$ は逆格子ベクトルである（B. I. Halperin, Jpn. J. Appl. Phys. **26**, Suppl. 26-3, 1913 (1987)）。

**第 2 項（バンド 2）**: 断面ごとに、$S(k_z)$ の上で滑らかなゲージを取れば、Stokes の定理で

$$
\int_{S(k_z)}\Omega_2^z\,d^2k=\oint_{\partial S(k_z)}{\bf a}_2\cdot d{\bf k}\equiv\gamma_2(k_z)
\tag{C36}
$$

となる。輪 $\partial S(k_z)$ はフェルミ面を $k_z$ で切った線（$z$ 向きの磁場の中のサイクロトロン軌道）で、$\gamma_2(k_z)$ はその上のベリー位相である。
ところが、フェルミ面の上の状態だけからは $\gamma_2(k_z)$ は $2\pi$ を法としてしか決まらない。ゲージ変換 $|u_{2{\bf k}}\rangle\to e^{i\phi({\bf k})}|u_{2{\bf k}}\rangle$ で ${\bf a}_2\to{\bf a}_2-\nabla_{\bf k}\phi$ となり、

$$
\gamma_2(k_z)\to\gamma_2(k_z)-2\pi w
\tag{C37}
$$

と変わるからである（$w$ は $\phi$ が輪を一周する間に $2\pi$ だけ何回まわるか、という整数）。どの $w$ が正しいか（$S(k_z)$ の内側で滑らかにつながるか）は、フェルミ面の内側の状態を見ないと決まらない。

まとめると、$\sigma_{ij}$ のうちフェルミ面から決まるのは、式 (C35) の形の項（逆格子ベクトル ${\bf G}$ の分）を法とした値だけで、その整数の部分はフェルミ海の情報である（Haldane の結論）。
例として、バンド 2（フェルミ面とその上の状態）はまったく同じで、バンド 1 だけが違う二つの系 A、B を考える。例えば金属の下に、結合させずに、普通の絶縁体（A、${\bf G}=0$）を置くか、$z$ 向きに周期 $c$ で積んだチャーン絶縁体の層（B、$C_1^z=1$）を置くかである。式 (C34) から

$$
\sigma_{xy}^{\rm B}-\sigma_{xy}^{\rm A}=\frac{e^2}{hc}
\tag{C38}
$$

で、この差はフェルミ面からは区別できない。極端な場合が三次元の量子ホール状態の絶縁体で、$S$ が空（フェルミ面が無い）なのに $\sigma_{xy}=C_1^z e^2/hc$ である。

弱い結合を入れてもフェルミ面とベリー位相は連続に変わるだけで、整数はギャップが閉じない限り変わらない。

ただし、表面のある有限の試料では、チャーン数だけの表面状態がギャップを横切り $E_F$ を通る（バルクとエッジの対応）。整数の部分も $E_F$ の表面状態が運ぶと見れば、「$E_F$ の状態で決まる」という見方は境界まで含めて保てる。
周期境界条件のバルクの式では、部分積分の境界項がブリルアン域の周期性で消えるので、整数の部分だけがフェルミ海に残る。
$\sigma_{xx}$ にはこの整数にあたるものが無いので、式 (C17) の部分積分で完全にフェルミ面に移せる。

## メモ: 位相の巻き付きと Thouless のポンプ {#thouless-pump}

**エッセンス**: 三次元の空間 $X$ の各点 $(x,y,z)$ に位相因子 $e^{i\theta(x,y,z)}$ を割り当てることを考える。値として意味を持つのは単位円（$T^1$）の上の点 $e^{i\theta}$ で、$\theta$ 自身は $2\pi$ の整数倍だけの不定さを持つ。
$X$ の中の周回路 $L$ を一周する間に $e^{i\theta}$ が単位円を何回まわるか、という巻き数

$$
w_L=\frac{1}{2\pi}\oint_L\nabla\theta\cdot d{\bf r}
\tag{C39}
$$

は整数である。$\theta$ は一周で元の値に戻る必要はなく、$e^{i\theta}$ が戻ればよいので、$2\pi w_L$ だけずれてよい。
巻き数は、割り当てや周回路を連続に変えても変わらない（整数は連続には変われない）。変わるのは、途中で $e^{i\theta}$ が定義できなくなる所（振幅が 0 になる所、特異点）を通るときだけである。
特に、周回路を（特異点を通らずに）一点に縮められるなら、巻き数は 0 になるしかないことに注意しておく。一点に縮んだ周回路では $\theta$ は一定で巻き数は 0 であり、縮める間に巻き数は変わらないからである。

ここで周回路が縮められるかどうかが効く。$X$ が普通の三次元空間なら、どの周回路も一点に縮められる。したがって、巻き数が 0 でない周回路は、縮める途中で必ず特異点を通る。つまり特異点（渦糸。三次元では線になる）を巻いている。
ところが $X$ が $T^3$（周期的な箱。[予備](#fermi-sea-surface)）なら、$x$、$y$、$z$ の向きに巻く周回路（$T^1$）は一点に縮められない。したがって、どこにも特異点が無くても、三つの向きの巻き数 $(w_x,w_y,w_z)$ を持てる。例えば辺 $a$、$b$、$c$ の周期的な箱で

$$
e^{i\theta}=e^{i{\bf G}\cdot{\bf r}},\qquad {\bf G}=2\pi\Big(\frac{w_x}{a},\frac{w_y}{b},\frac{w_z}{c}\Big)
\tag{C40}
$$

は、どこで見ても滑らかで特異点も無いのに、巻き数 $(w_x,w_y,w_z)$ を持つ。巻き数の組は逆格子ベクトルそのものである（式 (C35) の ${\bf G}$ と同じ形）。トポロジカルな量子化のエッセンスはこれである。

波動関数に当てはめると、ブロッホ関数の周期部分 $|u_{n{\bf k}}\rangle$ は、ブリルアン域という閉じた三次元空間（$T^3$）の各点 ${\bf k}$ に、振幅 × 位相の形で与えられている（$(x,y,z)$ を ${\bf k}$ と読む）。上の議論は、その振幅を固定して、位相因子 $e^{i\theta({\bf k})}$ だけを考えたものにあたる。ゲージ変換 $|u_{n{\bf k}}\rangle\to e^{i\phi({\bf k})}|u_{n{\bf k}}\rangle$ は振幅を変えずに位相だけを変えるもので、その巻き数が式 (C37) の $w$ である。

$T^3$ の上の関数 $e^{i\theta(x,y,z)}$ は、連続変形で移り合うものを同じとみなすと、三つの向きの周回路（$T^1$）の巻き数 $(w_x,w_y,w_z)$ だけで特徴づけられる。では閉じた面（$T^2$）ではどうか。
波動関数の位相は、どこで見ても滑らかなのに、$T^2$ 全体で一つの滑らかな関数には選べないことがある。そのときの食い違いの巻き数がチャーン数（式 (C34) の後の $C_1^z$）である。
$T^3$ には閉じた $T^2$ が 3 種類あるので、その上のチャーン数も三つ $(C_1^x,C_1^y,C_1^z)$ あり、式 (C35) の ${\bf G}$ を与える。
フェルミ面がある場合は、フェルミ面のトポロジーも関係してくる（例えば閉じたポケットのフェルミ面は、中にバンドの縮退点が無ければチャーン数を持てない）が、ここでは立ち入らない。

**Thouless のポンプ**（D. J. Thouless, Phys. Rev. B **27**, 6083 (1983)）: 一次元の絶縁体（格子定数 $a$）のハミルトニアンを、周期 $T$ でゆっくり変える（$H(t+T)=H(t)$）。占有バンドと非占有バンドの間のギャップは、一周期の間ずっと開いているとする。
各時刻 $t$ に、占有バンドのブロッホ関数の周期部分 $|u_{nk}(t)\rangle$ のベリー接続を $k$ の $T^1$（一次元のブリルアン域）の上で積分した位相（Zak の位相）を割り当てる。分極はこれで書ける（R. D. King-Smith and D. Vanderbilt, Phys. Rev. B **47**, 1651 (1993)）:

$$
\gamma(t)=\sum_{n\in{\rm occ}}\oint_{T^1}i\langle u_{nk}(t)|\partial_k u_{nk}(t)\rangle\,dk,\qquad
P(t)=-\frac{e}{2\pi}\,\gamma(t)
\tag{C41}
$$

（電子の電荷を $-e$ とした）。$\gamma$ はゲージ変換で $2\pi$ の整数倍だけ変わる（式 (C37)）ので、$P$ は $e$ を法としてしか決まらない。電子を一つ隣の単位胞に移しても同じ結晶であることに対応する。
断熱的に変えるとき、流れる電流は $dP/dt$ で、一周期で単位胞の境界を通る電荷は

$$
Q=\int_0^T\frac{dP}{dt}\,dt=-\frac{e}{2\pi}\big[\gamma(T)-\gamma(0)\big]=-eC
\tag{C42}
$$

である（$\gamma(t)$ は $t$ について連続にたどる）。$H(T)=H(0)$ なので $e^{i\gamma(T)}=e^{i\gamma(0)}$ で、$C$ は時間の周回路 $T^1$ から単位円への写像 $t\mapsto e^{i\gamma(t)}$ の巻き数 (C39) である。一周期で運ばれる電荷は $e$ の整数倍に量子化される。

ここには特異点は出てこない。ギャップはずっと開いていて、$|u_{nk}(t)\rangle$ は $(k,t)$ について滑らかである。整数は、$t$ が周期的で周回路が縮められないことだけから来る。
$C$ は $(k,t)$ の作る $T^2$ の上のチャーン数でもあり、式 (C35) の $C_1^z$ で $k_y$ を $t$ に置き換えたものと同じ種類の量である。パラメータの空間に広げて見ると、ハミルトニアンの描くループが、ループの外にあるギャップの閉じる点の周りを何回まわるかを数えている（例えば Rice–Mele 模型では、二量化と副格子のポテンシャルの差の描く円が、両方とも 0 の点を囲むとき $C=\pm 1$）。

端のある有限の鎖では、一周期の間に一方の端から他方の端へ電荷 $-eC$ が運ばれる。バルクにはギャップがあるので、これは端の状態の準位がギャップを横切ることで起きる（一周期の間に $C$ 本）。これが一次元版のバルクとエッジの対応である。二次元のホール伝導度でも、試料を筒にして磁束を通すと $k_y$ が $t$ の役をして同じ議論になる（Laughlin の議論）。
