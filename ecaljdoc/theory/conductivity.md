# 電気伝導度計算の基本 {#conductivity}

::: info 概要
**一体問題** として、線形応答（久保公式）から不純物による電気抵抗を導く。
ミソは、不純物の配置について平均を取るのが **グリーン関数の積** だという点である（一本ずつ平均してから掛けるのでは足りない）。
グリーン関数の符号の約束は [Basic I](./basic) と同じで、遅延グリーン関数は $G^R(E)=1/(E-H+i\delta)$、先進グリーン関数は $G^A(E)=1/(E-H-i\delta)$。
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

> 楠瀬の本（楠瀬博明『スピンと軌道の電子論』講談社 (2019)、§1.3 の式 (1.106)–(1.109)）では、電場を $-\Omega\,{\bf P}\cdot{\bf E}$（${\bf P}$ は分極）で入れ、時間についての部分積分を使って伝導度を書いている（本の式 (1.107)、(1.108)。本では $\Pi_{\mu\nu}=\chi^R_{JJ}$、$V=\Omega$）。式 (C8) と並べると
>
> $$
> \begin{aligned}
> \sigma(\omega)&=\frac{\chi^R_{JJ}(\omega)-ne^2\Omega/m}{i\omega\Omega}&&\text{（式 (C8)）}\\
> \sigma(\omega)&=\frac{\chi^R_{JJ}(\omega)-\chi^R_{JJ}(0)}{i\omega\Omega}&&\text{（本の式 (1.109)）}
> \end{aligned}
> \tag{C9}
> $$
>
> で、違いは第 2 項だけである。これは、式 (C6) の右辺第 1 項 $-(ne^2\Omega/m)A_x(t)$ を、あえて第 2 項の形 $-\chi^R_{JJ}(\omega=0)\,A_x(t)$ に押し込んで書いているためである。
> $\chi^R_{JJ}(\omega=0)\,A_x$ は、$A_x$ を無限にゆっくり（断熱的に）入れたときの密度行列 $\rho$ の変化による電流である。線形応答の式を導いた道を逆にたどると分かる。静的な $A_x$ を $e^{\delta t'}$ で入れ、式 (C3) を導いた $\delta\rho_{\rm I}$ の式で $F(t')=e^{\delta t'}A_x$ として $t=0$ で見ると
>
> $$
> \begin{aligned}
> \delta\rho&=\frac{i}{\hbar}\int_{-\infty}^{0}dt'\,e^{\delta t'}\,[\hat J^0_{x,{\rm H}}(t'),\rho_0]\,A_x,\\
> \langle n|\delta\rho|m\rangle&=(f_m-f_n)\,\frac{\langle n|\hat J^0_x|m\rangle}{E_n-E_m-i\delta}\,A_x,\\
> {\rm Tr}[\delta\rho\,\hat J^0_x]&=\sum_{nm}(f_n-f_m)\,\frac{|\langle n|\hat J^0_x|m\rangle|^2}{E_m-E_n-i\delta}\,A_x=\chi^R_{JJ}(0)\,A_x
> \end{aligned}
> \tag{C10}
> $$
>
> である（$\hbar\delta$ を改めて $\delta$ と書いた）。$i\delta$ はゆっくり入れたことの名残で、三行目は式 (C4) で $\omega=0$ としたものである。
> つまり本の形は、静的な $A_x$ の下では $\rho$ の変化による電流が反磁性の電流をちょうど打ち消し、全電流が 0 になると仮定していることになる。超伝導のように静的な $A_x$ の下で電流が残る系では、本の形は使えず、式 (C8) を使う。

### 直流の極限: Kubo–Greenwood の式 {#kubo-greenwood}

式 (C8) の実部を取ると、反磁性の項は落ちて

$$
\begin{aligned}
{\rm Re}\,\sigma(\omega)&=\frac{1}{\omega\Omega}\,{\rm Im}\,\chi^R_{JJ}(\omega)\\
&=\frac{1}{\omega\Omega}\,{\rm Im}\sum_{nm}(f_n-f_m)\,\frac{|\langle n|\hat J^0_x|m\rangle|^2}{E_m-E_n-\hbar\omega-i\delta}\\
&=\frac{\pi}{\omega\Omega}\sum_{nm}(f_n-f_m)\,|\langle n|\hat J^0_x|m\rangle|^2\,\delta(E_m-E_n-\hbar\omega)
\end{aligned}
\tag{C11}
$$

となる。2 行目では式 (C4) を $\hat A=\hat B=\hat J^0_x$ として入れ、3 行目では ${\rm Im}\,[1/(x-i\delta)]=\pi\delta(x)$ を使った。
$\omega\to 0$ では、デルタ関数により $E_m=E_n+\hbar\omega$ なので $f_n-f_m\simeq -f'(E_n)\,\hbar\omega$ で、$\omega$ が消える。
低温では $-f'(E)=\delta(E-E_F)$ である。$\hat J^0_x$ の一体の行列要素は $-e\langle n|\hat v_x|m\rangle$ で、スピンの因子を考えて 2 倍することにすると（$n$、$m$ は軌道の状態）、

$$
\sigma_{xx}=\frac{2\pi\hbar e^2}{\Omega}\sum_{nm}|\langle n|\hat v_x|m\rangle|^2\,\delta(E_F-E_n)\,\delta(E_F-E_m)
\tag{C12}
$$

となる。
ここでのデルタ関数は、幅 $\delta$ のローレンツ関数 $\delta_\delta(x)=(1/\pi)\,\delta/(x^2+\delta^2)$（添字の $\delta$ は式 (C11) の分母の $\delta$）の $\delta\to+0$ の極限と考える。$\delta$ を有限に残した考察が必要であるが、それは後で式 (C19)、(C20) で行う。
これは、フェルミ準位 $E_F$ での一体の量だけで決まる。
$\delta(E-H)=-{\rm Im}\,G^R(E)/\pi$ で書き直すと、Kubo–Greenwood の式

$$
\sigma_{xx}=\frac{2\hbar e^2}{\pi\Omega}\,
{\rm Tr}\big[\hat v_x\,{\rm Im}\,G^R(E_F)\,\hat v_x\,{\rm Im}\,G^R(E_F)\big]
\tag{C13}
$$

を得る。$\hat v_x$ は速度の演算子、${\rm Tr}$ は一体の（軌道の）状態についての和、因子 2 はスピンである。

> よりよい理解のために、多体の摂動論の見方との関係を見ておく。一体の系では、リングダイアグラム（$G$ 二つの積で、端点に $J$ を置いて閉じたもの）が厳密な結果を与える。
> その線形応答関数は、時間順序のグリーン関数 $G^T(E)=\sum_m|m\rangle\langle m|/\big(E-E_m+i\delta\,{\rm sgn}(E_m-\mu)\big)$ を使って
>
> $$
> \chi^T_{JJ}(\omega)=\frac{i}{2\pi}\int dE\,{\rm Tr}\big[J\,G^T(E)\,J\,G^T(E+\hbar\omega)\big]
> \tag{C14}
> $$
>
> と書ける（絶対零度。$J=-e\hat v_x$ は一体の演算子）。$E$ の積分で留数を取ると、占有の $n$ と非占有の $m$ の組だけが残る。極の位置は式 (C4) と同じで、違うのは $\omega\lt 0$ の側の $i\delta$ の向きだけである。
>
> ただし、式 (C14) の $G^T$ をそのまま $G^R$ に置き換えても遅延型にはならない。$G^T$ は、占有の状態の分は $G^A$、非占有の状態の分は $G^R$ の形をしているので、遅延型への書き換えは応答関数 $\chi$ の側で行う（$\omega$ の極をすべて下半面に移す。$\omega\to\omega+i\delta$）。

> 式 (C4) をグリーン関数で書き直す。$J_{nm}=\langle n|J|m\rangle$ と書き、$f_n=\int dE\,f(E)\,\delta(E-E_n)$ として、$f_n$ の和では $E=E_n$、$f_m$ の和では $E=E_m$ とおくと
>
> $$
> \begin{aligned}
> \chi^R_{JJ}(\omega)&=\sum_{nm}(f_n-f_m)\,\frac{J_{nm}J_{mn}}{E_m-E_n-\hbar\omega-i\delta}\\
> &=-\int dE\,f(E)\,{\rm Tr}\big[J\,\delta(E-H)\,J\,G^R(E+\hbar\omega)\\
> &\qquad+J\,G^A(E-\hbar\omega)\,J\,\delta(E-H)\big]
> \end{aligned}
> \tag{C15}
> $$
>
> となる（$1/(E_m-E-\hbar\omega-i\delta)=-G^R_m(E+\hbar\omega)$、$1/(E-E_n-\hbar\omega-i\delta)=G^A_n(E-\hbar\omega)$ を使った。添字は固有状態の成分）。
> 式 (C15) の一行目で $\hbar\omega=E_m-E_n$ となるエネルギーが励起エネルギーである（$n$ から $m$ への励起）。$\chi^R_{JJ}$ は、この励起に対して遅延型の形をしている（極 $\hbar\omega=E_m-E_n-i\delta$ は $\omega$ の下半面にある）。式 (C15) の最後の式変形は、この励起エネルギー $\hbar\omega=E_m-E_n$ を、遅延型のグリーン関数の準粒子エネルギーにすり替える形になっている。
> 負の振動数の極（$E_m\lt E_n$ の組）は、位相の回転が逆向きの摂動に対する応答を表している（時間は正の向きに進む）。式 (C15) の一行目を実時間で書き直すと分かる:
>
> $$
> \chi^R_{JJ}(t)=\frac{i}{\hbar}\,\theta(t)\sum_{nm}(f_n-f_m)\,J_{nm}J_{mn}\,e^{-i(E_m-E_n)t/\hbar-\delta t/\hbar}
> \tag{C16}
> $$
>
> 各項は $\theta(t)\,e^{\mp i|E_m-E_n|t/\hbar-\delta t/\hbar}$ の形で、$E_m\gt E_n$ なら正の振動数、$E_m\lt E_n$ なら負の振動数で回るが、どちらも $\theta(t)$ を持ち、$t\gt 0$ にだけ応答する。$\int dt\,e^{i\omega t}$ でフーリエ変換すると式 (C15) の一行目に戻る。
> 式 (C15)、(C16) のように、応答関数は必ず、$1/(\varepsilon-\hbar\omega-i\delta)$（振動数で書いた形）もしくはステップ関数 $\times\,e^{-i\varepsilon t/\hbar}$（実時間で書いた形）に、行列要素と占有数の因子を掛けたものの和で書ける（$\varepsilon=E_m-E_n$ は励起エネルギー）。

> 式 (C15) に $\delta(E-H)=-\big(G^R(E)-G^A(E)\big)/2\pi i$ を入れると、任意の温度で次のようになる:
>
> $$
> \begin{aligned}
> \chi^R_{JJ}(\omega)&=\frac{1}{2\pi i}\int dE\,f(E)\\
> &\quad\times{\rm Tr}\big[J\,\big(G^R(E)-G^A(E)\big)\,J\,G^R(E+\hbar\omega)\\
> &\qquad+J\,G^A(E-\hbar\omega)\,J\,\big(G^R(E)-G^A(E)\big)\big]
> \end{aligned}
> \tag{C17}
> $$
>
> となる。

> $G^R(E)=(E-H+i\delta)^{-1}$ と $G^A(E)=(E-H-i\delta)^{-1}$ の違いは、極（$H$ の固有値）を実軸の下に置くか上に置くかの違いである。その差 $G^R-G^A=-2\pi i\,\delta(E-H)$ は、その極の作るデルタ関数である。
> $G^R={\rm Re}\,G+i\,{\rm Im}\,G$、$G^A={\rm Re}\,G-i\,{\rm Im}\,G$（${\rm Re}\,G$、${\rm Im}\,G$ はどちらもエルミート。${\rm Im}\,G=-\pi\delta(E-H)$ は極だけの部分である）と分けると、虚部は ${\rm Re}\,G\cdot{\rm Im}\,G$ と ${\rm Im}\,G\cdot{\rm Im}\,G$ の組み合わせになる。
> ${\rm Re}\,G$ を含む項のトレースは実数なので虚部には効かず、${\rm Im}\,G$ どうしの積（極の寄与）だけが残る。

> 第 2 項で $E\to E+\hbar\omega$ とずらすと
>
> $$
> \begin{aligned}
> {\rm Re}\,\sigma(\omega)&=\frac{1}{\omega\Omega}\,{\rm Im}\,\chi^R_{JJ}(\omega)\\
> &=\frac{1}{\pi\omega\Omega}\int dE\,\big[f(E)-f(E+\hbar\omega)\big]\\
> &\quad\times{\rm Tr}\big[J\,{\rm Im}\,G^R(E)\,J\,{\rm Im}\,G^R(E+\hbar\omega)\big]
> \end{aligned}
> \tag{C18}
> $$
>
> となる。これは任意の $\omega$ と温度で成り立つ（$\omega\neq 0$ では光学伝導度を与える）。$\omega\to 0$、低温の極限では、式 (C12) と一致する。

式 (C12) のデルタ関数の積の和（積分）は分かりにくいと思うので、以下で説明しておく。
$\delta$ を有限に残すと、式 (C12) のデルタ関数は $\delta_\delta$ のままになる。
$\sum_n\delta_\delta(E-E_n)$ は状態密度 $D(E)$（系全体の、エネルギーあたりの状態の数）である。

すこし考察してみる。もし行列要素 $|\langle n|\hat v_x|m\rangle|^2$ が $E_F$ の近くで一様（平均 $\overline{|v_x|^2}$）なら

$$
\begin{aligned}
\sigma_{xx}&\simeq\frac{2\pi\hbar e^2}{\Omega}\,\overline{|v_x|^2}\,\Big[\sum_n\delta_\delta(E_F-E_n)\Big]^2\\
&=\frac{2\pi\hbar e^2}{\Omega}\,\overline{|v_x|^2}\,D(E_F)^2
\end{aligned}
\tag{C19}
$$

で、$\delta$ によらない。一方、対角成分だけの場合は $n=m$ の項だけが残り、$\int dx\,\delta_\delta(x)^2=1/(2\pi\delta)$ から

$$
\begin{aligned}
\sigma_{xx}&\simeq\frac{2\pi\hbar e^2}{\Omega}\sum_n|\langle n|\hat v_x|n\rangle|^2\,\delta_\delta(E_F-E_n)^2\\
&=\frac{\hbar e^2}{\Omega}\,\overline{v_x^2}\,\frac{D(E_F)}{\delta}
\end{aligned}
\tag{C20}
$$

で、$\delta\to 0$ で発散する（散乱が無ければ電流は減衰しない）。大雑把には、$\delta$ はフォノンなどによる準粒子の寿命 $\tau$ から来る幅 $\hbar/2\tau$ と考えてもよい。
実際には、行列要素は対角成分だけではないにせよ、$n$ とつながる $m$ の数は（完全系から $\sum_m|\langle n|\hat v_x|m\rangle|^2=\langle n|\hat v_x^2|n\rangle$ が有限なので）有限であり、式 (C20) が正しいと考えてよい。

式 (C20) の「フェルミ面の電子が電気伝導に効く」は、Drude の「すべての電子が $\tau$ の間加速される」と、次のようにつながる。
式 (C20) の右辺の次の因子は、一つのバンド、スピンあたりで、$v_x\,\delta(E_F-\varepsilon_{\bf k})=-\hbar^{-1}\partial_{k_x}\theta(E_F-\varepsilon_{\bf k})$（$\theta$ は占有の階段関数、$v_x=\hbar^{-1}\partial\varepsilon_{\bf k}/\partial k_x$）を使って $k_x$ について部分積分すると（境界の項はブリルアン域の周期性で消える）

$$
\begin{aligned}
\frac{\overline{v_x^2}\,D(E_F)}{\Omega}&=\int\frac{d^3k}{(2\pi)^3}\,v_x^2\,\delta(E_F-\varepsilon_{\bf k})\\
&=\frac{1}{\hbar}\int\frac{d^3k}{(2\pi)^3}\,\theta(E_F-\varepsilon_{\bf k})\,\frac{\partial v_x}{\partial k_x}\\
&=\int_{\rm occ}\frac{d^3k}{(2\pi)^3}\,\frac{1}{\hbar^2}\frac{\partial^2\varepsilon_{\bf k}}{\partial k_x^2}
\end{aligned}
\tag{C21}
$$

となる。$v_x^2$ の一方の $v_x$ はデルタ関数と組になって $\theta$ の微分になり（$\theta$ を $k_x$ で微分すると $\partial\varepsilon_{\bf k}/\partial k_x=\hbar v_x$ が出る）、部分積分で微分がもう一方の $v_x$ に移って $\partial v_x/\partial k_x=\hbar^{-1}\partial^2\varepsilon_{\bf k}/\partial k_x^2$（逆有効質量）になる。
左辺はフェルミ面の量、右辺は占有された全状態の逆有効質量の和である。
自由電子では $\partial^2\varepsilon_{\bf k}/\partial k_x^2=\hbar^2/m$ が一定なので、右辺は $n/2m$ となる（$n$ は両スピンを含む電子の数密度で、$n/2$ がスピンあたり）。
したがって、式 (C20) で $\delta=\hbar/2\tau$ とし、式 (C21) を使うと

$$
\sigma_{xx}=2e^2\tau\,\frac{\overline{v_x^2}\,D(E_F)}{\Omega}=\frac{ne^2\tau}{m}
\tag{C22}
$$

となり、[Drude の式](#drude) が出る。
このフェルミ面での有効質量を一定とする近似は、Na、Cu、Al などではよく成り立つはずである。

> 占有された領域の積分とフェルミ面の上の積分が部分積分で結ばれるという式 (C21) の構造は、ベリー曲率を使った異常ホール伝導度の計算（占有状態でのベリー曲率の積分を、フェルミ面の上のベリー位相の積分に書き直す）と同じである。


## 不純物の配置についての平均 {#impurity-average}

式 (C13) は、不純物の配置を一つ決めたときの式である。${\rm Im}\,G^R=(G^R-G^A)/2i$ と書くと、$G^A=(G^R)^\dagger$ から ${\rm Tr}[\hat v_xG^R\hat v_xG^A]$ は実数、${\rm Tr}[\hat v_xG^A\hat v_xG^A]$ は ${\rm Tr}[\hat v_xG^R\hat v_xG^R]$ の複素共役なので、式 (C23) の一行目が配置ごとに成り立つ（近似は使っていない）。観測される伝導度は、式 (C13) を不純物の配置について平均したもので、二行目になる（$G$ はすべて $E=E_F$）:

$$
\begin{aligned}
{\rm Tr}\big[\hat v_x\,{\rm Im}\,G^R\,\hat v_x\,{\rm Im}\,G^R\big]&=\frac{1}{2}\Big\{ {\rm Tr}\big[\hat v_x\,G^R\,\hat v_x\,G^A\big]\\
&\qquad-{\rm Re}\,{\rm Tr}\big[\hat v_x\,G^R\,\hat v_x\,G^R\big]\Big\},\\
\sigma_{xx}&=\frac{\hbar e^2}{\pi\Omega}\Big\{\Big\langle{\rm Tr}\big[\hat v_x\,G^R\,\hat v_x\,G^A\big]\Big\rangle_{\rm imp}\\
&\qquad-{\rm Re}\,\Big\langle{\rm Tr}\big[\hat v_x\,G^R\,\hat v_x\,G^R\big]\Big\rangle_{\rm imp}\Big\}
\end{aligned}
\tag{C23}
$$

二つの項を、以下でまとめて展開する（第 2 項が小さいことは [Drude の式](#drude) の節で見る）。
グリーン関数の積の平均を計算するには、そもそもの一体問題に戻る。$G^R(E)=(E-H_0-V+i\delta)^{-1}$ は、不純物の配置を一つ決めた一体のハミルトニアン $H=H_0+V$（$H_0$ は運動エネルギー、$V$ は式 (C1)）のグリーン関数である。これを $V$ について展開し、積の平均を、$V$ の対の平均の積に分解する。

::: tip 準備: 不純物がランダムに配置されているときの平均
**平均を取る前**（配置を一つ決めたとき）: 平面波 $|{\bf k}\rangle$ での行列要素は $\langle{\bf k}|V|{\bf k}'\rangle=\Omega^{-1}u({\bf k}-{\bf k}')\sum_j e^{-i({\bf k}-{\bf k}')\cdot{\bf R}_j}$ である。

**平均を取った後**: 不純物の位置は乱雑なので、位相 $e^{-i({\bf k}-{\bf k}')\cdot{\bf R}_j}$ は平均すると打ち消し合い（ランダム位相）、一次の平均は ${\bf k}\neq{\bf k}'$ で 0 になり、二次では同じ不純物 $j$ について位相が消える組だけが残る:

$$
\begin{aligned}
&\big\langle \langle{\bf k}|V|{\bf k}'\rangle\big\rangle_{\rm imp}=n_i\,u(0)\,\delta_{ {\bf k}{\bf k}'},\\
&\big\langle \langle{\bf k}|V|{\bf k}'\rangle\langle{\bf k}'|V|{\bf k}\rangle\big\rangle_{\rm imp}=\frac{n_i}{\Omega}\,|u({\bf k}-{\bf k}')|^2
\end{aligned}
\tag{C24}
$$

:::

これを使って式 (C23) の二つの項を展開する。$G^{\pm}=(E-H_0-V\pm i\delta)^{-1}$（$G^+=G^R$、$G^-=G^A$）と書き、${\rm Tr}[\hat v_xG^{\pm}\hat v_xG^{\pm'}]$ を考える（二つの複号は独立に取る。$(+,-)$ が第 1 項、$(+,+)$ が第 2 項）。平面波 $|{\bf k}\rangle$ では $\hat v_x$ は対角（$v_x({\bf k})=\hbar k_x/m$）なので、トレースは

$$
\begin{aligned}
{\rm Tr}\big[\hat v_x G^{\pm}\hat v_x G^{\pm'}\big]
&=\sum_{ {\bf k}{\bf k}'}v_x({\bf k})\,G^{\pm}_{ {\bf k}{\bf k}'}\,v_x({\bf k}')\,G^{\pm'}_{ {\bf k}'{\bf k}}
\end{aligned}
\tag{C25}
$$

である（$G_{ {\bf k}{\bf k}'}=\langle{\bf k}|G|{\bf k}'\rangle$。$E=E_F$ は書かない）。$G^{\pm}$ を $V$ について展開すると

$$
\begin{aligned}
G^{\pm}_{ {\bf k}{\bf k}'}={}&G^{\pm}_0({\bf k})\,\delta_{ {\bf k}{\bf k}'}+G^{\pm}_0({\bf k})\,V_{ {\bf k}{\bf k}'}\,G^{\pm}_0({\bf k}')\\
&+\sum_{ {\bf k}_1}G^{\pm}_0({\bf k})\,V_{ {\bf k}{\bf k}_1}\,G^{\pm}_0({\bf k}_1)\,V_{ {\bf k}_1{\bf k}'}\,G^{\pm}_0({\bf k}')+\cdots
\end{aligned}
\tag{C26}
$$

である。ここで $G^{\pm}_0({\bf k})=1/(E-\varepsilon_{\bf k}\pm i\delta)$、$\varepsilon_{\bf k}=\hbar^2k^2/2m$、$V_{ {\bf k}{\bf k}'}=\langle{\bf k}|V|{\bf k}'\rangle$ である。
これを式 (C25) に入れて平均する。$V$ について 1 次の項は ${\bf k}={\bf k}'$ の成分 $n_iu(0)$ だけが残るが、これは一様なエネルギーのずれで、$E_F$ に繰り込む。2 次の項は、二つの $V$ がどこにあるかで三つに分かれ、式 (C24) を使うと

$$
\begin{aligned}
&\Big\langle{\rm Tr}\big[\hat v_x G^{\pm}\hat v_x G^{\pm'}\big]\Big\rangle_{\rm imp}\\
&\quad=\sum_{\bf k}v_x({\bf k})^2\,G^{\pm}_0({\bf k})\,G^{\pm'}_0({\bf k})\\
&\qquad+\color{#1f6fd0}{\sum_{\bf k}v_x({\bf k})^2\,G^{\pm}_0({\bf k})\,\Sigma^{\pm}_2({\bf k})\,G^{\pm}_0({\bf k})\,G^{\pm'}_0({\bf k})}\\
&\qquad+\color{#1f6fd0}{\sum_{\bf k}v_x({\bf k})^2\,G^{\pm}_0({\bf k})\,G^{\pm'}_0({\bf k})\,\Sigma^{\pm'}_2({\bf k})\,G^{\pm'}_0({\bf k})}\\
&\qquad+\color{#d26a00}{\sum_{ {\bf k}{\bf k}'}v_x({\bf k})\,G^{\pm}_0({\bf k})\,G^{\pm'}_0({\bf k})}\\
&\qquad\quad\color{#d26a00}{\times\frac{n_i}{\Omega}\,|u({\bf k}-{\bf k}')|^2\,G^{\pm}_0({\bf k}')\,G^{\pm'}_0({\bf k}')\,v_x({\bf k}')}\\
&\qquad+O(V^4)
\end{aligned}
\tag{C27}
$$

となる。ここで $\Sigma^{\pm}_2({\bf k})=\sum_{ {\bf k}_1}(n_i/\Omega)\,|u({\bf k}-{\bf k}_1)|^2\,G^{\pm}_0({\bf k}_1)$ である。
式 (C27) の第 2 項と第 3 項（<span style="color:#1f6fd0">青</span>）は、一本の $G$ の中の二つの $V$ を組にした項で、$\Sigma_2$ を一つ挟んだ形である（一本の $G$ の中の平均は、ランダム位相のため ${\bf k}={\bf k}'$ のときだけ残るので、${\bf k}$ について対角になる）。$V$ について高次まで見ると、$\Sigma_2$ を二つ、三つと挟んだものが出てくる。対が互いに交差しないものだけを取る（Born 近似）と、それらは等比級数になって

$$
\begin{aligned}
\bar G^{\pm}({\bf k},E)&=G^{\pm}_0+G^{\pm}_0\Sigma^{\pm}G^{\pm}_0+G^{\pm}_0\Sigma^{\pm}G^{\pm}_0\Sigma^{\pm}G^{\pm}_0+\cdots\\
&=\frac{1}{E-\varepsilon_{\bf k}-\Sigma^{\pm}({\bf k},E)}
\end{aligned}
\tag{C28}
$$

にまとまる（$\Sigma^{\pm}=\Sigma^{\pm}_2$）。$\bar G$ は並進対称で、${\bf k}$ について対角である。
第 4 項（<span style="color:#d26a00">橙</span>）は、$G^{\pm}$ の $V$ と $G^{\pm'}$ の $V$ を組にした項で、二つの $G$ の積を平均したことで残る、同じ不純物で散乱される相関である（一本ずつ平均してから掛けると出てこない）。この型の項を落とすと、式 (C23) は

$$
\begin{aligned}
\sigma_{xx}=\frac{\hbar e^2}{\pi\Omega}\sum_{\bf k}v_x({\bf k})^2\Big\{&\bar G^R({\bf k},E_F)\,\bar G^A({\bf k},E_F)\\
&-{\rm Re}\,\bar G^R({\bf k},E_F)^2\Big\}
\end{aligned}
\tag{C29}
$$

となる。以下、$\Sigma$ を求め（次節）、式 (C29) の二つの項を評価する（[Drude の式](#drude) の節）。第 4 項の型の項は [バーテックス補正](#vertex-correction) の節で取り入れる。対が互いに交差する組み方は弱局在の補正を与えるが、ここでは扱わない。

## 第 2 項と第 3 項の足し上げ: 寿命と Drude の式 {#self-energy-sum}

### 式 (C29) の右辺第 1 項: 寿命と Drude の式 {#averaged-g}

式 (C28) の $\Sigma^{\pm}$ を具体的に求める。$\Sigma_2$ の和を積分にすると、自己エネルギーは

$$
\Sigma^R({\bf k},E)=n_i\int\frac{d^3k'}{(2\pi)^3}\,
\frac{|u({\bf k}-{\bf k}')|^2}{E-\varepsilon_{\bf k'}+i\delta}
\tag{C30}
$$

で、その虚部が時間 $\tau$ を与える（実部はエネルギーのずれで、$E_F$ に繰り込む）:

$$
{\rm Im}\,\Sigma^R=-\frac{\hbar}{2\tau},
\tag{C31}
$$

式 (C28) は、不純物散乱を一体問題として（Born 近似の範囲で）すべて取り入れたグリーン関数の表式になっている。その極は実軸から $\hbar/2\tau$ だけ離れ、$G$ は時間とともに $e^{-t/2\tau}$ で減衰する（確率 $|G|^2$ は $e^{-t/\tau}$）。したがって、この $\tau$ は準粒子寿命とみなせる。
ここで $E=E_F$ とおくと

$$
\frac{1}{\tau}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,\delta(E_F-\varepsilon_{\bf k'})
\tag{C32}
$$

$1/\tau$ はフェルミの黄金則による散乱の速さそのものである（[Basic I](./basic#td-vs-ti-perturbation) の、黄金則と自己エネルギーの虚部の対応）。

これを式 (C29) に入れると、$\bar G^R({\bf k},E_F)=1/(E_F-\varepsilon_{\bf k}+i\hbar/2\tau)$、$\bar G^A$ はその複素共役である。第 1 項の $\bar G^R\bar G^A=1/\big[(E_F-\varepsilon_{\bf k})^2+(\hbar/2\tau)^2\big]$ は $\varepsilon_{\bf k}=E_F$ に幅 $\hbar/2\tau$ の鋭い山を持ち、${\bf k}$ で和を取ると $\tau$ に比例する。第 2 項の $(\bar G^R)^2$ は、極が実軸の同じ側にあるのでそうならない。以下で具体的に評価する。

<a id="drude"></a>

式 (C29) では、$\bar G$ に自己エネルギーによる幅 $\hbar/2\tau$ が入っている。
$\tau$ が長い（$\hbar/\tau\ll E_F$）とき

$$
\begin{aligned}
\bar G^R({\bf k},E)\,\bar G^A({\bf k},E)
&=\frac{1}{(E-\varepsilon_{\bf k})^2+(\hbar/2\tau)^2}\\
&\simeq\frac{2\pi\tau}{\hbar}\,\delta(E-\varepsilon_{\bf k})
\end{aligned}
\tag{C33}
$$

で、これは式 (C20) のローレンツ関数の幅 $\delta$ を $\hbar/2\tau$ にしたものにあたる。したがって式 (C29) の第 1 項は

$$
\sigma_{xx}=2e^2\tau\,\frac{\overline{v_x^2}\,D(E_F)}{\Omega}
\tag{C34}
$$

となり、式 (C22) と同じ Drude の式（自由電子では $ne^2\tau/m$）が出る。

### 式 (C29) の右辺第 2 項 {#rr-term}

式 (C29) の第 2 項を正確に見積もる。$g(\varepsilon)=\sum_{\bf k}v_x({\bf k})^2\,\delta(\varepsilon-\varepsilon_{\bf k})$、$\gamma=\hbar/2\tau$ と書くと、

$$
\begin{aligned}
\sum_{\bf k}v_x({\bf k})^2\,\bar G^R({\bf k},E_F)^2
&=\int d\varepsilon\,\frac{g(\varepsilon)}{(E_F-\varepsilon+i\gamma)^2}\\
&=-\int d\varepsilon\,\frac{g'(\varepsilon)}{E_F-\varepsilon+i\gamma}
\end{aligned}
\tag{C35}
$$

である。二つ目の等号は、$(E_F-\varepsilon+i\gamma)^{-2}=\partial_\varepsilon(E_F-\varepsilon+i\gamma)^{-1}$ を使った部分積分である（境界の項は消える）。右辺は $\gamma\to0$ でも有限で、その実部は $-{\rm P}\!\int d\varepsilon\,g'(\varepsilon)/(E_F-\varepsilon)$ となり、$\tau$ によらない。$g$ は $E_F$ ほどのエネルギーの尺度で変わるので、その大きさは $g(E_F)$ の程度である。
一方、第 1 項は、式 (C33) を $g$ で書くと

$$
\begin{aligned}
\sum_{\bf k}v_x({\bf k})^2\,\bar G^R({\bf k},E_F)\,\bar G^A({\bf k},E_F)&=\int d\varepsilon\,\frac{g(\varepsilon)}{(E_F-\varepsilon)^2+\gamma^2}\\
&\simeq\frac{\pi}{\gamma}\,g(E_F)=\frac{2\pi\tau}{\hbar}\,g(E_F)
\end{aligned}
\tag{C36}
$$

で、$\tau$ に比例する。式 (C35) と (C36) の比は $\hbar/(E_F\tau)$ の程度で、散乱が弱い（$\hbar/\tau\ll E_F$）ときは小さいので、第 2 項を落とした。

### まとめ {#lifetime-summary}

まとめると、第 4 項（橙）の型を落とした近似では、式 (C29) の第 1 項から寿命 $\tau$ の Drude の式 (C34) が出て、第 2 項はそれより $\hbar/(E_F\tau)$ の程度だけ小さい。
素朴に $\delta\to 0$ とすると抵抗が 0 になったが、不純物による散乱の幅が $\delta$ の役を担い、不純物は寿命 $\tau$ を通して抵抗を与える。すなわち、電子は電場で加速され、不純物で散乱されるまでの時間 $\tau$ で寿命が決まるという、素朴な Drude の描像の結果を得る。

## バーテックス補正: 輸送緩和時間 {#vertex-correction}

積の平均を、式 (C27) の第 4 項（橙）の型の項を全次数まで足し上げて取る（はしご型。Born 近似の自己エネルギー (C30) と組になる近似）。
電流のバーテックス $\Gamma_x({\bf k})$ は

$$
\begin{aligned}
\Gamma_x({\bf k})={}&v_x({\bf k})+n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\\
&\times\bar G^R({\bf k}',E_F)\bar G^A({\bf k}',E_F)\,\Gamma_x({\bf k}')
\end{aligned}
\tag{C37}
$$

を満たす。$\Gamma_x({\bf k})=\gamma\,v_x({\bf k})$ と置き、(C33) を使うと、${\bf k}'$ の向きについての平均で $v_x({\bf k}')$ が $v_x({\bf k})\cos\theta_{\bf kk'}$ になるので、
$\gamma=1+\gamma\,(1-\tau/\tau_{\rm tr})$、すなわち $\gamma=\tau_{\rm tr}/\tau$ を得る。ここで

$$
\frac{1}{\tau_{\rm tr}}=\frac{2\pi}{\hbar}\,n_i\int\frac{d^3k'}{(2\pi)^3}\,|u({\bf k}-{\bf k}')|^2\,
\big(1-\cos\theta_{\bf kk'}\big)\,\delta(E_F-\varepsilon_{\bf k'})
\tag{C38}
$$

は輸送緩和時間である。式 (C34) の片方の $v_x$ を $\Gamma_x$ に替えると

$$
\sigma_{xx}=2e^2\tau_{\rm tr}\,\frac{\overline{v_x^2}\,D(E_F)}{\Omega}
\tag{C39}
$$

となる（自由電子では $ne^2\tau_{\rm tr}/m$）。

## 二つの緩和時間: 寿命 $\tau$ と輸送緩和時間 $\tau_{\rm tr}$ {#two-taus}

重要なのは、寿命 $\tau$（式 (C32)）と輸送緩和時間 $\tau_{\rm tr}$（式 (C38)）の違いである。二つを並べると、<strong>違いは被積分関数の重み $1-\cos\theta_{\bf kk'}$ だけである。</strong>

| | 準粒子寿命 $\tau$（式 (C32)） | 輸送緩和時間 $\tau_{\rm tr}$（式 (C38)） |
| --- | --- | --- |
| 別の言い方 | 一粒子緩和時間、量子寿命 | 運動量緩和時間、輸送寿命 |
| 散乱の重み | $1$ | $1-\cos\theta_{\bf kk'}$ |
| どこから出るか | 一本の $G$ の平均（自己エネルギー、式 (C27) の第 2、3 項） | $G$ の積の平均（式 (C27) の第 4 項の型の足し上げ） |
| 何の時間か | 電子が状態 ${\bf k}$ に留まる時間。どの向きに散乱されても ${\bf k}$ を失う | 電流（運動量の向き）が失われる時間 |
| どこに現れるか | スペクトル関数の幅 $\hbar/\tau$ | 伝導度 $\sigma_{xx}$（式 (C39)） |
| 散乱断面積の理論では | 全断面積 $\sigma_{\rm tot}$ | 輸送断面積（運動量移行断面積）$\sigma_{\rm tr}$ |

散乱断面積の理論で書くと、次のようになる。一つの不純物による微分断面積を $d\sigma/d\Omega$ とすると（Born 近似では $d\sigma/d\Omega=(m/2\pi\hbar^2)^2|u({\bf q})|^2$、$|{\bf q}|=2k_F\sin(\theta/2)$）、式 (C32) と式 (C38) は

$$
\begin{aligned}
&\frac{1}{\tau}=n_i\,v_F\,\sigma_{\rm tot},\qquad \sigma_{\rm tot}=\int d\Omega\,\frac{d\sigma}{d\Omega},\\
&\frac{1}{\tau_{\rm tr}}=n_i\,v_F\,\sigma_{\rm tr},\qquad \sigma_{\rm tr}=\int d\Omega\,(1-\cos\theta)\,\frac{d\sigma}{d\Omega}
\end{aligned}
\tag{C40}
$$

と書ける。全断面積 $\sigma_{\rm tot}$ はどの向きへの散乱も数え、ビームから粒子が抜ける率にあたる。輸送断面積 $\sigma_{\rm tr}$ は運動量の移行 $\hbar k_F(1-\cos\theta)$ で重みを付けたもので、抵抗や拡散に出るのはこちらである。$\tau$ はスペクトル関数の幅に出る（例えば de Haas–van Alphen 振動の Dingle 温度から求める寿命で、量子寿命とも呼ばれる）。Born 近似を超えるときは、$|u|^2$ を一つの不純物の T 行列に替えればよく、この対応はそのまま成り立つ。

- **前方への散乱**（$\theta\simeq0$）は、電子を別の状態 ${\bf k}'$ に移すので $1/\tau$ には効くが、運動量の向きはほとんど変えないので $1/\tau_{\rm tr}$ にはほとんど効かない。
- **後方への散乱**（$\theta\simeq\pi$）は、運動量の向きを反転させるので、$1/\tau_{\rm tr}$ に 2 倍の重みで効く。
- したがって $1-\cos\theta\le2$ から $\tau_{\rm tr}\ge\tau/2$ で、前方への散乱が強いほど $\tau_{\rm tr}$ は $\tau$ より長くなる。
- $u$ が点状（$u({\bf q})$ が ${\bf q}$ によらない）なら、散乱は向きによらず $\cos\theta$ の平均は 0 で、$\tau_{\rm tr}=\tau$ になる（バーテックス補正は消える）。遮蔽されたクーロン型の不純物のように $u({\bf q})$ が小さい ${\bf q}$ に偏るほど、前方への散乱が強くなり、$\tau_{\rm tr}\gg\tau$ となりうる。

式 (C29) のように一本ずつ平均してから掛けると、$\sigma_{xx}$ に入るのは $\tau$ になる（式 (C34)）。正しい $\tau_{\rm tr}$ を得るには、$G$ の積の平均（バーテックス補正）が要る。これが、平均を取るのは $G$ の積であることの意味である。
$\tau_{\rm tr}$ を使った式 (C39) は、Boltzmann 方程式で得られるものと同じである。

<hr style="margin-top:4em;margin-bottom:2.5em">

不純物の話は以上である。以下は、異常ホール効果などの計算のためのメモである。

## メモ: ブリルアン域での積分とブリルアン域のトポロジー $T^3$ {#bz-topology}

線形応答では、ブリルアン域の上の積分量を計算する。そこで、ブリルアン域の幾何学（トポロジー）を確認しておく。二次元のトーラスで理解しておけば、直感的にそれなりに理解できる。これがトポロジカルな物理量を議論するときのキモである。いきなり計算に入ってしまうと、よく分からなくなる。

${\bf k}$ と ${\bf k}+{\bf G}$（${\bf G}$ は逆格子ベクトル）は同じ点を表すので、ブリルアン域は向かい合う面を貼り合わせた三次元のトーラス $T^3$ である。
$T^2$（二次元のトーラス）とは、ドーナツの表面のような二次元の空間のことである。$T^3$（三次元のトーラス）とは、$T^2$ の三次元版で、境界の無い閉じた三次元の空間の一種である。
$T^2$ には、一点に縮めることのできない周回路（$T^1$）が 2 種類ある（連続変形で移り合う周回路は同じとみなす）。ドーナツの大きい円を「輪」、断面の小さい円を「管」と呼ぶと、輪に沿って回るものと、管を回るものである。ほかの周回路はこの 2 種類の組み合わせで、例えば、管を 2 周回ったあと輪に沿って 1 周回る経路などである。
$T^3$ では、そのような周回路（$T^1$）が $k_x$、$k_y$、$k_z$ の向きに巻く 3 種類あり、さらに一点に縮めることのできない閉じた面（$T^2$）が 3 種類ある。これらは境界の無い閉じた面（二次元のトーラス）であり、$k_x=$一定 などで指定される。
例えば $k_z=$一定 の断面 $T^2(k_z)$ は、ブリルアン域を $k_x$、$k_y$ の向きに巻く。このような面の上の積分は、境界の上の積分に書き直すことができない。
一方、ブリルアン域の中に収まる領域（下のポケット $S$ のようなもの）は境界を持ち、その上の積分は Stokes の定理で境界の上の積分に書き直せる。

なぜそもそもこんなことを考えないといけないかと言うと、線形応答の計算において（トポロジカルな応答として）問題点が発生する場合があるからである。線形応答の積分は、状態についての和の形では位相の選び方によらない。しかし、部分積分を用いて、その積分をベリー曲率やフェルミ面の上の積分に書き直すときには、位相を滑らかに選べることを使う。$T^3$ の中の閉じた面 $T^2$ の上では、それができないことがある。

単純化して考えると、$T^3$ の各点 $(x,y,z)$ に位相因子 $e^{i\theta(x,y,z)}$ が指定されている場合に着目すればよい。振幅 × 位相の形で与えられる場合もあるが、振幅の部分は本質的でなく、あとから考えればよい。

位相は連続的かつ滑らかに変化するものとする。ただし、特異線（渦糸。近づく向きによって位相の極限が違うので、そこで位相が定まらない。周りを一周すると位相が $2\pi$ の整数倍まわる。例えば $z$ 軸の上の $\theta=\arg(x+iy)$）や、バンドの縮退する特異点も、場合により考える。三次元では、位相の特異は普通は線になる。

周回路 $L$ を一周する間に $e^{i\theta}$ が単位円を何回まわるかを、巻き数

$$
w_L=\frac{1}{2\pi}\oint_L\nabla\theta\cdot d{\bf r}
\tag{C41}
$$

という（整数）。巻き数は連続変形では変わらず、変わるのは特異線を横切るときだけである。したがって、一点に縮められる周回路の巻き数は 0 しかない。普通の三次元空間では、0 でない巻き数は、その周回路が特異線を巻いていることを意味する。
ところが $T^3$ では、$x$、$y$、$z$ の向きに巻く周回路（$T^1$）が縮められないので、特異線が無くても巻き数を持てる。例えば辺 $a$、$b$、$c$ の周期的な箱で

$$
e^{i\theta}=e^{i{\bf G}\cdot{\bf r}},\qquad {\bf G}=2\pi\Big(\frac{w_x}{a},\frac{w_y}{b},\frac{w_z}{c}\Big)
\tag{C42}
$$

を考える。$x$ を $a$ だけ進めると位相は $2\pi w_x$ 進み、$e^{i\theta}$ はもとに戻る（$y$、$z$ も同じ）。したがって、これは周期的な箱（$T^3$）の上の関数として、どこで見ても滑らかである。それでも、$x$ の向きに一周すると位相は $2\pi w_x$ まわるので、巻き数 $(w_x,w_y,w_z)$ を持つ。${\bf G}$ は逆格子ベクトルで、巻き数はその成分（$2\pi/a$、$2\pi/b$、$2\pi/c$ を単位とする）である。トポロジカルな量子化のエッセンスはこれである。$T^3$ から単位円への関数 $e^{i\theta(x,y,z)}$ は、この巻き数で分類される。巻き数の同じ関数どうしだけが、連続変形でつながりうる。
波動関数（$(x,y,z)$ を ${\bf k}$ と読む）では、ゲージ変換 $e^{i\phi({\bf k})}$ の巻き数がこれにあたる（後の式 (C47) の $w$）。さらに、位相を閉じた $T^2$ 全体で一つの滑らかな関数に選べないことがあり、その食い違いの巻き数がチャーン数である（後の式 (C44)、(C45)）。フェルミ面がある場合は、フェルミ面のトポロジーも関係してくる（例えば閉じたポケットのフェルミ面は、中にバンドの縮退点が無ければチャーン数を持てない）が、ここでは立ち入らない。

## メモ: フェルミ面項とフェルミ海項 {#fermi-sea-surface}

久保公式の静的な伝導度 $\sigma_{ij}$ を、重み $-\partial f/\partial E$ を持つ **フェルミ面項** と、重み $f(E)$ で占有状態全体を積分する **フェルミ海項** に分ける書き方がある（Středa の分け方。異常ホール効果の文献でよく使われる）。
フェルミ海項は $i\leftrightarrow j$ について反対称なので、$\sigma_{xx}$ には効かない。式 (C23) の第 1 項は、フェルミ面項のうち $G^RG^A$ の部分である。

ただし、この分け方は一つに決まらない。部分積分で一方を他方に書き直せるからである。
式 (C21) では、占有領域の積分（Drude の $n/m$）とフェルミ面での $\overline{v_x^2}D(E_F)$ が同じ量だった。
異常ホール伝導度でも、占有状態のベリー曲率の積分は、フェルミ面の上のベリー位相に書き直せる（F. D. M. Haldane, Phys. Rev. Lett. **93**, 206602 (2004)）。
したがって「フェルミ面の電子が効くか、全部の電子が効くか」という区別は、物理としては意味を持たない。

例外は整数の部分である。次の設定で式で見る。
三次元の結晶（簡単のため直方晶で、格子定数を $a$、$b$、$c$ とする）で、バンド 1 は完全に占有され（フェルミ面が無い）、バンド 2 はブリルアン域の中の領域 $S$（電子のポケット）だけが占有されているとする。フェルミ面は $S$ の境界の閉曲面 $\partial S$ である。
ブロッホ関数の周期部分 $|u_{n{\bf k}}\rangle$ のベリー接続とベリー曲率を

$$
\begin{aligned}
&{\bf a}_n({\bf k})=i\langle u_{n{\bf k}}|\nabla_{\bf k}u_{n{\bf k}}\rangle,\qquad
{\bf\Omega}_n({\bf k})=\nabla_{\bf k}\times{\bf a}_n({\bf k}),\\
&\Omega_n^z=-2\,{\rm Im}\sum_{m\neq n}\frac{\langle u_{n{\bf k}}|\partial_{k_x}H_{\bf k}|u_{m{\bf k}}\rangle\langle u_{m{\bf k}}|\partial_{k_y}H_{\bf k}|u_{n{\bf k}}\rangle}{(E_{n{\bf k}}-E_{m{\bf k}})^2}
\end{aligned}
\tag{C43}
$$

とする（$H_{\bf k}=e^{-i{\bf k}\cdot{\bf r}}He^{i{\bf k}\cdot{\bf r}}$）。二行目は状態についての和の形で、久保公式からそのまま出る形である。内因性の異常ホール伝導度は（符号は ${\bf\Omega}$ の向きの約束による）、$k_z$ ごとの断面に分けて

$$
\begin{aligned}
\sigma_{xy}&=\frac{e^2}{\hbar}\sum_n\int_{\rm occ}\frac{d^3k}{(2\pi)^3}\,\Omega_n^z\\
&=\frac{e^2}{h}\,\frac{1}{2\pi}\int dk_z\,\Big[\frac{1}{2\pi}\int_{T^2(k_z)}\Omega_1^z\,d^2k+\frac{1}{2\pi}\int_{S(k_z)}\Omega_2^z\,d^2k\Big]
\end{aligned}
\tag{C44}
$$

と書ける。$T^2(k_z)$ は $k_z$ を止めたブリルアン域の断面、$S(k_z)$ はその面でのポケットの断面である。

被積分関数 $\Omega_n^z$ は、どちらの項も位相の選び方によらない。位相を $|u_{n{\bf k}}\rangle\to e^{i\phi({\bf k})}|u_{n{\bf k}}\rangle$ と変えるとベリー接続は ${\bf a}_n\to{\bf a}_n-\nabla_{\bf k}\phi$ と変わるが、$\nabla_{\bf k}\times\nabla_{\bf k}\phi=0$ だからである。式 (C43) の二行目でも、各状態の位相が分子で打ち消し合っている。位相の選び方に依るのは、ベリー接続 ${\bf a}_n$ だけである。
バンドが縮退していると、縮退したバンドを混ぜる選び方があるので、$\Omega_n$ はバンドごとには決まらない。占有バンドについての和 $\sum_n\Omega_n$ は決まる。

**第 1 項（バンド 1）**: 閉じた面 $T^2(k_z)$ の上の積分 $C_1^z=\frac{1}{2\pi}\int_{T^2(k_z)}\Omega_1^z\,d^2k$ はチャーン数で、整数である。バンド 1 にはギャップがあるので $k_z$ によらない。$T^2(k_z)$ は閉じた面なので、これを境界の上の積分には書き直せない。
$k_z$ の積分は $2\pi/c$ を与えるので、第 1 項は $(e^2/h)\,C_1^z/c$ である。ブリルアン域を他の向きに巻く面のチャーン数 $C_1^x$、$C_1^y$ も合わせると

$$
\sigma_{ij}^{(1)}=\frac{e^2}{2\pi h}\,\epsilon_{ijl}\,G_l,\qquad
{\bf G}=2\pi\Big(\frac{C_1^x}{a},\frac{C_1^y}{b},\frac{C_1^z}{c}\Big)
\tag{C45}
$$

で、${\bf G}$ は逆格子ベクトルである（B. I. Halperin, Jpn. J. Appl. Phys. **26**, Suppl. 26-3, 1913 (1987)）。

**第 2 項（バンド 2）**: 断面ごとに、$S(k_z)$ の上で滑らかなゲージを取れば、Stokes の定理で

$$
\int_{S(k_z)}\Omega_2^z\,d^2k=\oint_{\partial S(k_z)}{\bf a}_2\cdot d{\bf k}\equiv\gamma_2(k_z)
\tag{C46}
$$

となる。輪 $\partial S(k_z)$ はフェルミ面を $k_z$ で切った線（$z$ 向きの磁場の中のサイクロトロン軌道）で、$\gamma_2(k_z)$ はその上のベリー位相である。
ところが、フェルミ面の上の状態だけからは $\gamma_2(k_z)$ は $2\pi$ を法としてしか決まらない。ゲージ変換 $|u_{2{\bf k}}\rangle\to e^{i\phi({\bf k})}|u_{2{\bf k}}\rangle$ で ${\bf a}_2\to{\bf a}_2-\nabla_{\bf k}\phi$ となり、

$$
\gamma_2(k_z)\to\gamma_2(k_z)-2\pi w
\tag{C47}
$$

と変わるからである（$w$ は $\phi$ が輪を一周する間に $2\pi$ だけ何回まわるか、という整数）。どの $w$ が正しいか（$S(k_z)$ の内側で滑らかにつながるか）は、フェルミ面の内側の状態を見ないと決まらない。

まとめると、フェルミ面の上の情報（フェルミ面とその上のベリー位相）だけからは、$\sigma_{ij}$ は式 (C45) の形の項（逆格子ベクトル ${\bf G}$ の分）の不定さを残してしか決まらない。その整数の部分は、フェルミ面を持たないバンド 1 のような、フェルミ海の情報である（Haldane の結論）。
例として、バンド 2（フェルミ面とその上の状態）はまったく同じで、バンド 1 だけが違う二つの系 A、B を考える。例えば金属の下に、結合させずに、普通の絶縁体（A、${\bf G}=0$）を置くか、$z$ 向きに周期 $c$ で積んだチャーン絶縁体の層（B、$C_1^z=1$）を置くかである。式 (C44) から

$$
\sigma_{xy}^{\rm B}-\sigma_{xy}^{\rm A}=\frac{e^2}{hc}
\tag{C48}
$$

で、この差はフェルミ面からは区別できない。極端な場合が三次元の量子ホール状態の絶縁体で、$S$ が空（フェルミ面が無い）なのに $\sigma_{xy}=C_1^z e^2/hc$ である。

弱い結合を入れてもフェルミ面とベリー位相は連続に変わるだけで、整数はギャップが閉じない限り変わらない。

周期境界条件のバルクの式では、部分積分の境界項がブリルアン域の周期性で消えるので、整数の部分だけがフェルミ海に残る。
$\sigma_{xx}$ にはこの整数にあたるものが無いので、式 (C21) の部分積分で完全にフェルミ面に移せる。

