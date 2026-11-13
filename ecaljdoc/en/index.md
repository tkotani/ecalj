---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "ecaljdoc (English)"
  tagline: This is for ecalj package, a first-principles electronic-structure calculations. The Japanese pages are the reference; pages not yet translated open in the Japanese site.
  image:
    src: /ecalj_logo.png
    alt: ecalj
  actions:
    - theme: brand
      text: GetStarted with AI (New!)
      link: /manual/getstartedAI
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
      text: PMT method (English)
      link: /en/theory/pmt
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

## For developers

The development notes are not on this site. Start from `MD/README.md` in the [ecalj repository](https://github.com/tkotani/ecalj)
(reading order, rules, and jumps by topic: finite temperature, MLO, magnons, GPU).
Claude Code starts from `CLAUDE.md` at the top of the repository. An AI assistant guiding a user starts from [GetStarted with AI](/manual/getstartedAI).
