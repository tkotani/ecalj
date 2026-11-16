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
