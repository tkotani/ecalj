---
# https://vitepress.dev/reference/default-theme-home-page
layout: home
pageClass: home-dark

hero:
  name: "ecaljdoc (English)"
  tagline: This is for ecalj package, a first-principles electronic-structure calculations. Every page has a Japanese and an English version; the Japanese pages are the reference.
  image:
    src: /ecalj_logo.png
    alt: ecalj
  actions:
    - theme: brand
      text: GetStarted with AI!
      link: /en/manual/getstartedAI
    - theme: alt
      text: What's new (2026-11)
      link: /en/manual/whatsnew
    - theme: alt
      text: MainDocument
      link: /en/manual/README_tutorial
    - theme: alt
      text: Input file
      link: /en/manual/toml_migration
    - theme: alt
      text: Theory
      link: /en/theory/gw
    - theme: alt
      text: PMT method
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

> **Work in progress**: many parts of this site are still incomplete and are being written.

## For developers

ecalj is developed by T. Kotani and co-workers, with Claude Code since 2026-02.

The development notes are not on this site. Start from `MD/README.md` in the [ecalj repository](https://github.com/tkotani/ecalj)
(reading order, rules, and jumps by topic: finite temperature, MLO, magnons, GPU).
Claude Code starts from `CLAUDE.md` at the top of the repository. An AI assistant guiding a user starts from [GetStarted with AI](/en/manual/getstartedAI).

## License

ecalj is distributed under the [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.html). For the terms, the third-party software included in ecalj, and how to cite ecalj in publications, see the [License](/en/manual/license) page (the reference; [LICENCE](https://github.com/tkotani/ecalj/blob/main/LICENCE) of the repository is its summary). It applies to AI systems (AI assistants, coding agents) without exception.
