import { defineConfig } from 'vitepress'
import { withMermaid } from "vitepress-plugin-mermaid";
import MarkdownItMathjax3 from 'markdown-it-mathjax3';

export default withMermaid({
  base : '/ecaljdoc/',
  title: "ecaljdoc: ",
  vite: { publicDir: 'public' },
  description: "document of ecalj",
  cleanUrls: true,
  ignoreDeadLinks: true,
  lastUpdated: true,
  // markdown: {
  //   math: true,
  // },
  markdown: {
    config: (md) => {
      md.use(MarkdownItMathjax3, {
        tex: {
          inlineMath: [['$', '$'], ['\\(', '\\)']],
          displayMath: [['$$', '$$'], ['\\[', '\\]']],
        },
        options: {
          useCDN: true, // MathJaxをCDNから読み込む
        },
      });
    },
  },
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: 'Home', link: '/' },
    ],

    search: {
      provider: 'local'
    },

    sidebar: [
      {
        items: [
          { text: 'MainDocument', link: '/manual/README_tutorial' },
        ]
      },
      {
        text: 'Install',
        items: [
          { text: 'Install', link: '/install/install' },
          { text: 'Install@ISSP', link: '/install/installISSP' }
        ]
      },
      {
        text: 'Manual',
        items: [
          { text: 'DFT  calculation : lmf', link: '/manual/lmf' },
          { text: 'QSGW calculation : gwsc', link: '/manual/gwsc' },
          { text: 'GPU version of QSGW',link: '/manual/ecaljgpu.md' },
          // { text: 'Density of states', link: '/manual/dos' },
          // { text: 'Band dispersion plot', link: '/manual/band' },
          // { text: 'Spin-orbit interaction', link: '/manual/soi' },
          // { text: 'Spin response function', link: '/manual/chipm' },
          { text: 'Dielecric function', link: '/manual/optical' },
          // { text: 'Structural optimization', link: '/manual/opts' },
          // { text: 'Quasi-particle\'s life time', link: '/manual/lifetime' },
          // { text: 'Spectrum function', link: '/manual/spec' },
          { text: 'ecalj auto', link: '/manual/auto' }
        ]
      },
      {
        text: 'Input files',
        items: [
          { text: 'ctrl.foobar', link: '/manual/lmf' },
          { text: 'GWinput', link: '/manual/gwinput' },
          { text: 'syml.foobar', link: '/manual/syml' }
        ]
      },
      {
        text: 'Theory',
        items: [
          // { text: 'PMT', link: '/theory/pmt' },
          // { text: 'Spin-orbit interaction', link: '/theory/soi' },
          { text: 'GW approximation', link: '/theory/gw' },
          { text: 'Quasi-particle self-consistent GW', link: '/theory/qsgw' },
          { text: 'Basis', link: '/theory/gw_basis' },
          { text: 'Product basis', link: '/theory/gw_product_basis' },
          { text: 'Projection of product basis', link: '/theory/gw_zmel' },
          { text: 'Optical properties', link: '/theory/optical_properties' },
          { text: 'Basic  I', link: '/theory/basic' },
          { text: 'Basic II', link: '/theory/basic2' },
        ]
      },
      {
        text: 'Implementation',
        items: [
          { text: 'hrcxq/hx0fp', link: '/implementation/hx0fp' },
          { text: 'hsfp0_sc', link: '/implementation/hsfp0' },
        ]
      },
    ],

    socialLinks: [
      { icon: 'github', link: 'https://github.com/tkotani/ecalj' },
      { icon: 'github', link: 'https://github.com/ecalj/ecaljdoc' },
    ]
  },
})
