---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "ecaljdoc"
  tagline: This is for ecalj package, a first-principles electronic-structure calculations
  actions:
    - theme: alt
      text: What's new (2026-10)
      link: /manual/whatsnew
    - theme: alt
      text: Tutorial
      link: /manual/README_tutorial
    - theme: alt
      text: TOML migration (2026-05)
      link: /manual/toml_migration
    - theme: alt
      text: Theory
      link: /theory/gw
    - theme: alt
      text: QSGW database (under construction)
      link: https://github.com/tkotani/DOSnpSupplement/blob/main/bandpng.md#band-structure--total-dos
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

## For developers / 開発者へ

The development notes are not on this site. Start from `MD/README.md` in the [ecalj repository](https://github.com/tkotani/ecalj)
(reading order, rules, and jumps by topic: finite temperature, MLO, magnons, GPU).
Claude Code starts from `CLAUDE.md` at the top of the repository.

開発の記録はこのサイトには無い。[ecalj のリポジトリ](https://github.com/tkotani/ecalj)の `MD/README.md` から読む
（読む順、決まり、話題ごとのジャンプ: 有限温度、MLO、マグノン、GPU）。Claude Code は最上位の `CLAUDE.md` から始める。
