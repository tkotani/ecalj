In japanese.2023-12-13

誘電関数ｑ＝０を計算するには、最近のバージョンの（２０２３−１０−２３以後）の
収束したところからスクリプトepsPP0を使ってください。
ecalj/MATERIALS/Cueps
にサンプルがあるので、Cuを収束させた後、
epsPP0 cu -np 4
などとすると、epsinter.glt, epsintra.glt, epsall.glt
ができます。epsintraはフェルミ面の寄与になります。

いくつかのｑ点で計算してｑ＝＞０を取るのです。そのためにctrlには
ーーーーーーーーーーーーーーーーーー
QforEPSunita on
<QforEPS>
 0d0 0d0 0.00001
 0d0 0d0 0.001
 0d0 0d0 0.0014142
 0d0 0d0 0.002
 0d0 0d0 0.0028284
 0d0 0d0 0.004
</QforEPS>
ーーーーーーーーーーーーーーーーーー
などと書いておく必要があります。


最低でも、2行入ります。
QforEPSunita on
<QforEPS>
 0d0 0d0 0.00001
 0d0 0d0 0.001
</QforEPS>
ーーーーーーーーーーーーーーーーーー
が、いります.できれば後何行かあったほうがいいです。

QforEPSunita onは単位を2pi/alat　とする指定です。
最初の行の0d0 0d0 0.00001　は数値誤差を計算するためにいるんです。

 0d0 0d0 0.001での誘電関数計算にも 0d0 0d0 0.00001で計算した行列要素（誤差の大きさを計算）を使うんです。なので、EPS0002以後が意味のある答えになります。
epsinter.gltなどではqを複数計算して重ねてみています。
Cuだとキレイにかさなっているのが見て取れます。

epsPP0では最後にreadeps.pyを読んでて答えをまとめ上げてます。
（readeps.py汚いです。試行錯誤の結果が残っていて関数fd0はいまは使っていません）。
いぜんよりはだいぶとキレイに求まると思います。

Ag、結果がこれでも変なら教えてください。

誘電関数計算にはバンド吸収端で(E-E0)**.5でキレイに立ち上がらない(GaAsなどの場合）、
という数値計算上の問題があります。メッシュを細かく取るときれいにできるんで、
必要なことだけ細かく取らないといかんのです。


kotani
