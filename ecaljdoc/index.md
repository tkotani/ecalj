---
# https://vitepress.dev/reference/default-theme-home-page
layout: home
pageClass: home-dark

hero:
  name: "ecaljdoc (Japanese)"
  tagline: これは第一原理計算の ecalj パッケージのページです。
  image:
    src: /ecalj_logo.png
    alt: ecalj
  actions:
    - theme: brand
      text: GetStarted with AI!
      link: /manual/getstartedAI
    - theme: alt
      text: What's new (2026-11)
      link: /manual/whatsnew
    - theme: alt
      text: メインドキュメント
      link: /manual/README_tutorial
    - theme: alt
      text: 入力ファイル
      link: /manual/toml_migration
    - theme: alt
      text: Theory
      link: /theory/gw
    - theme: alt
      text: PMT 法
      link: /theory/pmt
    - theme: alt
      text: GW1500 database (QSGW80)
      link: https://github.com/tkotani/DOSnpSupplement#readme
    - theme: brand
      text: Download ecalj package
      link: https://github.com/tkotani/ecalj
    - theme: brand
      text: Source of this document (ecalj/ecaljdoc)
      link: https://github.com/tkotani/ecalj/tree/main/ecaljdoc

features:
  - title: Quasiparticle self-consistent GW methods (QSGW)
    details: For electronic excitations. With GPUs.
  - title: Material properties on top of QSGW
    details: Dielectric functions, Spin fluctuations, Impact ionization rate, Spectrum functions of the Green's function, Modelling.
  - title: Easy to use
    details: Based on PMT=LMTO+LAPW. Virtually no headaching settings by hand.
---

> **加筆中**: このサイトにはまだ不完全な部分が多く、書き足している途中です。

## For developers / 開発者へ

The development notes are not on this site. Start from `MD/README.md` in the [ecalj repository](https://github.com/tkotani/ecalj)
(reading order, rules, and jumps by topic: finite temperature, MLO, magnons, GPU).
Claude Code starts from `CLAUDE.md` at the top of the repository. An AI assistant guiding a user starts from [GetStarted with AI](/manual/getstartedAI).

ecalj は小谷岳生と共同研究者が開発している（2026-02 から Claude Code と一緒に）。

開発の記録はこのサイトには無い。[ecalj のリポジトリ](https://github.com/tkotani/ecalj)の `MD/README.md` から読む
（読む順、決まり、話題ごとのジャンプ: 有限温度、MLO、マグノン、GPU）。Claude Code は最上位の `CLAUDE.md` から始める。

## ライセンス / License

ecalj は [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.html) で配布する。条件、ecalj に含まれる他のソフトウェア、論文での引用は [ライセンス](/manual/license) のページを見てほしい（このページが正本。リポジトリの [LICENCE](https://github.com/tkotani/ecalj/blob/main/LICENCE) はその要約）。AI（AI アシスタント、コーディングエージェント）にも例外なく適用される。
