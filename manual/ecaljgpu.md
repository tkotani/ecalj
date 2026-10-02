# ecalj GPU version の使い方

> ⚠️ GPU 版 `gwsc` も `ctrlg.<sname>.toml` だけを読みます。Legacy `ctrl.<sname>` / `GWinput` directory は `Legacy2toml.py <sname>` で変換してから実行してください。See [TOML migration](./toml_migration).

ecalj ではGPUによるQSGW計算が可能です。
* GPUは自己エネルギーおよび遮蔽クーロン相互作用の計算に使用されます。
* 現状では`lmf` (DFT計算) はCPUで実行されます。
> [!IMPORTANT]
> CPUで実行される部分とGPUで実行される部分で、それぞれ並列数を指定する必要があります。

* GPU の使用には NVIDIA のコンパイラでecaljのGPUバージョンをコンパイルする必要があります。
Kugui でのインストール、初期設定、GPUテスト計算は[installISSP](../install/installISSP.md) を参照してください。

## おすすめ設定 in ISSP Kugui
* ACCノード(GPU搭載ノード)を1ノード使用する。 1ノードに4GPUが搭載されています。
* CPUで実行される部分のMPI並列数(`-np` で指定)を64。
* GPUで実行される部分のMPI並列数(`-np2` で指定)をGPU数と同じ4とする。

以下qsub ファイルの例。

```bash job.sh
#!/bin/sh
#PBS -q i1accs
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1

ulimit -s unlimited
gwsc -np 64 -np2 4 --gpu 2 inas2gasb2 > lgwsc
```

`i1accs` はテストキュー(最大時間30分)ですので、プロダクトランでは`F1accs`等を使用してください。

`ctrlg.<sname>.toml` の `[gw]` セクション抜粋
```toml
[gw]
t_tetrakbt   = 0            # (K) required: 0 = no smearing of chi0; -T = Im chi0 smeared as wide as Fermi-Dirac at T; T > 0 = finite-T tetrahedron (manual/kBT)
KeepEigen = false
zmel_batch_gb = 4
```

* `[gw].KeepEigen = false`
* `[gw].zmel_batch_gb = 4`

それぞれの意味は [gwinput](../manual/gwinput.md) を参照してください。

### 混合精度: `--mp` と `--fp32`

`--mp` と `--mp --fp32` は `--prec=tf32` と `--prec=fp32` の古い書き方で、今も同じ意味で使えます
（対応は下の「精度の選び方と行列演算の自動選択」の節の表）。

`tf32` で入力の仮数を 10 ビットに落とすのは Σc の最後の積だけで、$W$ を作る側（$\chi_0$ の和、誘電行列の逆）は FP32 で計算します。
$W$ を作る側まで仮数 10 ビットにすると、**重元素 + 分子アニオン (NO3 / N3 / ClO / BrO /
PO4 等) で局所場が強く ε(q,ω) の条件数が大きい系**では ~1e-3 の誤差が
$(1-v\chi_0)^{-1}$ で増幅され、$W$ / $\Sigma^c$ が破壊されて QSGW が発散
(`lmf failed`) もしくは NaN (`lgw` / `lqpe`) になります。その例と再現データは
[`Samples/mptf32problem/`](https://github.com/tkotani/ecalj/tree/main/Samples/mptf32problem)
にあります。ストレージは `tf32` でも `fp32` でも単精度なので、メモリ使用量は同じです。

```bash
gwsc 5 -np 64 -np2 4 --gpu --prec=fp32 inas2gasb2 > lgwsc
```

### 精度の選び方と行列演算の自動選択 (2026-09)

`gwsc --prec=tf32|fp32|fp64` で GW 部分の精度を選ぶ。旧来の `--mp`、`--mp --fp32` と同じもので、どちらの書き方でもよい。

| `--prec` | 実行ファイル | 単精度の積 | 旧来の書き方 |
| --- | --- | --- | --- |
| `fp64` | 倍精度版（`hgw_gpu` など） | — | `--gpu` |
| `fp32` | 混合精度版（`hgw_mp_gpu` など） | FP32 | `--gpu --mp --fp32` |
| `tf32` | 混合精度版 | Σc の積だけ入力の仮数 10 ビット（TF32 か FP16、表が選ぶ）、ほかは FP32 | `--gpu --mp` |

`tf32` で精度を下げるのは Σc の最後の積（W と zmel の積、zsec）だけ。W を作る側（χ0 の和、誘電行列の逆）の精度を下げると、
小さな誤差が $(1-v\chi_0)^{-1}$ で拡大される。LiTi2O4 6³ の `hgw` 1 回（kt1、RTX 5090 ×2、2026-09-27 夜の版）で:

| | 時間 | 倍精度との差（E_F ±1 eV の Re Σc、最大） |
| --- | --- | --- |
| fp32 | 367 秒 | 0.02 meV |
| tf32、Σc の積を FP16 で（RTX 5090 の表が選ぶ） | 173 秒 | 1.3 meV |
| tf32、Σc の積を TF32 で（`--linalg=tf32.cgemm.large=cublas,tf32.cgemm.small=cublas`） | — | 0.8 meV |
| すべての積を TF32（2026-09-27 までの `--mp`） | — | 5 meV（NiO のギャップが 0.1 eV ずれる） |

QSGW の 1 反復（`gwsc`、LiTi2O4 6³、tf32、GPU 2 枚・CPU 60 本）は約 212 秒で、その 8 割が `hgw`。

`--mp` だけの指定も `tf32` と同じにした。tf32 の誤差の大半は Σc の一様な縮み（Re Σc の −5×10⁻⁵〜−1×10⁻⁴ 倍）で、
テンソルコアの和の取り方（途中の和の切り捨て）から来る。入力の丸めは最近接で、これ以上は安くは減らせない。
tf32 は途中の反復に使い、最後の反復は `--prec-final=fp32:N` で fp32 にするのが前提（誤差は乱雑ではなく一様な縮みなので、最後の fp32 の反復で消えるはず。
**2026-09-28 の時点ではこの使い方を実際に走らせて確かめていない**。6³ では tf32 だけの 10 反復でもバンドは fp32 と rms 0.3 meV で一致した）。

`--prec-final=fp32:2` のように書くと、その回の最後の 2 反復だけ別の精度で回す（途中は tf32、最後は fp32、など）。

行列の積と誘電行列の逆行列を GPU で**どの方法で計算するか**は、精度ごとに表で決まる。ユーザーが選ぶのは精度だけ。

- 表は `InstallAll.py --gpu` の最後に `linalgtune_gpu` がその GPU で測って `<bindir>/ecalj_linalg_policy.toml` に書く
  （GPU が他のジョブで使われているときは測らずに飛ばす。あとで手で `linalgtune_gpu` を実行すればよい）
- 表が無い、または GPU の名前が違うときは、すべて cuBLAS（以前と同じ動き）
- 使われた表は `hgw` などの標準出力の最初に 1 回出る（`linalg policy (fp32, from …)`）
- 方法の候補: 複素の積は `cublas`、`realsgemm`（複素の積を実数 SGEMM 1 回に組み替える）、`realhgemm`（その FP16 版。入力は FP16、
  和は FP32。仮数は TF32 と同じ 10 ビットで、RTX 5090 では TF32 の 1.8 倍速い。精度から `tf32` の行でだけ選ばれる）、`gemmul8`（Ozaki 法、INT8）。
  逆行列は `lu64`（FP64 の LU）、`mixed1` / `mixed2`（FP32 の LU ＋ Newton 1 回 / 2 回）。計測ツールは、その精度の誤差の約束
  （FP32: 相対 1e-5、FP64: 1e-12 など）を満たすものの中で、cuBLAS より 1 割以上速いものだけを選ぶ
- 手で上書きするときは `--linalg=fp32.cgemm.large=realsgemm,fp32.epsinv.large=mixed1` のように行を並べる（`gwsc` のオプションに書けば全部に渡る）
- `tf32` の Σc の積は表の `tf32` の行で決まる。`--linalg=tf32.cgemm.large=gemmul8:6` とすると Ozaki の 6 分解になり、
  Σc の積の精度が TF32 の 3e-4 から 1e-5 に上がる（そのぶん遅い）。Ozaki は分解数 1 つで精度が約 1 桁変わる
- 例（RTX 5090 ×2、LiTi2O4 の `hgw` 1 回）: 表は fp32 で「単精度の複素積は組み替え、倍精度の積は GEMMul8 14、逆行列は mixed1」、
  tf32 の行（Σc の積）は FP16 の組み替え（`realhgemm`）になる。6³ は以前の 833 秒が fp32 367 秒・tf32 173 秒（2026-09-27 夜）、
  9³ は 5004 秒が fp32 2990 秒・tf32 1975 秒（同日午前の版、四面体の重みの補助あり）。倍精度の計算との差も Re Σc で 2.4e-4 → 4.2e-5 eV に減る
- GEMMul8 は小さな積（どれかの辺が 64 未満、または m·n·k < 1e8）には表に関係なく使わない（分解の固定費が積より大きい）

### 同じ機械で複数のジョブを走らせるとき

`gwsc` は GPU を使う実行ファイルを走らせる前に、GPU ごとのロック（`/tmp/ecalj_res/gpu<N>.lock`、`flock`）を必要な枚数だけ取り、
取れた GPU だけを `CUDA_VISIBLE_DEVICES` で渡す。空きが足りなければ待つ（何を待っているかを表示する）。

| 環境変数 | 意味 |
| --- | --- |
| `CUDA_VISIBLE_DEVICES` | 設定してあれば、その GPU の中からだけ取る |
| `ECALJ_GPU_WAIT=<秒>` | その秒数で待つのをやめて止まる（既定は待ち続ける） |
| `ECALJ_GPU_LOCK=0` | ロックしない |

ロックはプロセスが終われば外れるので、落ちたジョブのロックは残らない。

### 画面出力とログ (2026-09-27)

GW のプログラム（`hvccfp0`、`hsfp0_sc`、`hgw` など）は、ランク 0 の画面出力をそのまま `gwsc` のログ（`lvcc`、`lsxC`、`lgw` など）に書く。
途中経過（q 点ごとの `Time:`、(q, k) ごとの進み具合）とメモリの使用量（`Memused (CPU) … (GPU) …`）はそこで追える（`tail -f lgw`）。
ほかのランクの画面出力は捨てる。ランクごとに見たいときは `--fullstdo`（`gwsc` に書けば全部に渡る）で `stdout.<rank>.<prog>` が出て、
`gwsc` がそれを `STDOUT/` に移す。どのランクのエラーも、ランク番号付きで標準エラー（`gwsc` の出力）に出る。
（2026-09-27 までは、全ランク（ランク 0 も）が `stdout.<rank>.<prog>` を書いていた。）

### 四面体の重みを CPU で先に計算する

`gwsc --gpu` は、GPU が `hvccfp0` と `hgw` を回している間に、空いている CPU コア（`-np` の数）で同じ `hgw` を GPU なしで走らせ（`--tetwt_write`）、
各 q の四面体の重みを `__TETWT.<iq>.<isp>` に書く。`hgw` はそれを読み、重みの計算（LiTi2O4 6³ で 1 q あたり約 6 秒）を省く。
結果は自分で計算した場合とビット一致する。ファイルが無い・合わないときは自分で計算する。使わないときは `--no-tetwt-helper`。

### MPSの設定

以下を`~/.bashrc` に記載する。
```bash ~/.bashrc
if which nvidia-cuda-mps-control > /dev/null 2>&1 ; then
  export CUDA_MPS_PIPE_DIRECTORY=$(pwd)/nvidia-mps-$(hostname)
  export CUDA_MPS_LOG_DIRECTORY=$(pwd)/nvidia-log-$(hostname)
  echo "start nvidia-cuda-mps-control at" $(hostname)
  nvidia-cuda-mps-control -d
fi
```

--------------

# Using the GPU Version

The GPU can be utilized in QSGW calculations:
- The GPU is used for the calculation of self-energy and screened Coulomb interaction.
- `lmf` is executed on the CPU.

## Requirements:
1. NVIDIA GPU
2. Compiler: NVHPC with GPU-accelerated libraries: cuBLAS, cuSolver

## Compilation
```bash
python3 InstallAll.py --fc nvfortran --gpu
```

## Execution
- The GPU version can be used by specifying the `--gpu` option with `gwsc`.
- The number of MPI processes for the GPU code can be specified with `-np2`. This should typically match the number of GPUs available.
Example:
```
gwsc -np 64 -np2 4 --gpu 1 $id > lgwsc
```

## Tips

1. The GPU version uses fewer MPI processes compared to the CPU version, so the available memory per MPI process is larger
   This allows for larger batch sizes in the calculation, potentially improving computation speed.
   The variable controlling the batch size is `zmel_batch_gb` in the `[gw]` section of `ctrlg.<sname>.toml`.
   For GPU calculations, set the value to around 4 (representing 4GB).
   ```toml
   [gw]
   zmel_batch_gb = 4
   ```
   > [!TIP]
   > If not set (or set below 2), 2 GB is used in GPU calculations.

2. When handling large systems, add the following to `[gw]` to prevent memory exhaustion:
   ```toml
   [gw]
   KeepEigen = false
   ```

## Notes
- It is recommended to run GPU calculations on a single node.
- When using multiple nodes, ensure that MPI processes are correctly distributed across nodes.
