---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "ecaljdoc"
  tagline: This is for ecalj package, a first-principles electronic-structure calculations
  actions:
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
      text: Download ecaljdoc (this document)
      link: https://github.com/ecalj/ecaljdoc

features:
  - title: Quasiparticle self-consistent GW methods (QSGW)
    details: For electronic excitations. With GPUs.
  - title: Material properties on top of QSGW
    details: Dielectric functions, Spin fluctuations, Impact ionization rate, Spectrum functions of the Green's function, Modelling.
  - title: Easy to use
    details: Based on PMT=LMTO+LAPW. Virtually no headaching settings by hand.
