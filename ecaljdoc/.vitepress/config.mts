import { defineConfig } from 'vitepress'
import { withMermaid } from "vitepress-plugin-mermaid";
import MarkdownItMathjax3 from 'markdown-it-mathjax3';

// The sidebar is written once (2026-11-13); toEn() makes the English one from it.
const sidebarRoot = [
  {
    items: [
      { text: 'GetStarted with AI (New!)', link: '/manual/getstartedAI' },
      { text: 'MainDocument', link: '/manual/README_tutorial' },
    ]
  },
  {
    text: 'Install',
    items: [
      { text: 'Install', link: '/install/install' },
      { text: 'Install@ISSP', link: '/install/installISSP' },
      { text: 'Uninstall', link: '/install/uninstall' }
    ]
  },
  {
    text: 'Manual',
    items: [
      { text: '更新履歴（要約）/ What\'s new', link: '/manual/whatsnew' },
      { text: 'TOML migration (2026-05)', link: '/manual/toml_migration' },
      { text: 'Samples (where what lives)', link: '/manual/samples' },
      { text: 'DFT  calculation : lmf part', link: '/manual/lmf' },
      { text: 'MT 半径の上限と偽の準位', link: '/manual/mtradius' },
      { text: 'Command-line options (--foo) catalogue', link: '/manual/cmdopts' },
      { text: 'QSGW calculation: gwsc part', link: '/manual/gwsc' },
      { text: 'GPU version of QSGW',link: '/manual/ecaljgpu.md' },
      { text: 'kBT — 有限温度と均し方 (t_tetrakbt / t_sigmaw)', link: '/manual/kBT' },
      { text: 'MLO — 最大局在化などに代わる自動モデル化法', link: '/manual/mlo' },
      { text: 'MLO-gwsc — Sigma を MLO 表現で内挿する QSGW（開発中）', link: '/manual/mlo_gwsc' },
      { text: 'UsageDetailed ', link: '/manual/UsageDetailed' },
      // { text: 'Density of states', link: '/manual/dos' },
      // { text: 'Band dispersion plot', link: '/manual/band' },
      // { text: 'Spin-orbit interaction', link: '/manual/soi' },
      // { text: 'Spin response function', link: '/manual/chipm' },
      { text: 'Dielecric function', link: '/manual/optical' },
      // { text: 'Structural optimization', link: '/manual/opts' },
      // { text: 'Quasi-particle\'s life time', link: '/manual/lifetime' },
      // { text: 'Spectrum function', link: '/manual/spec' },
      // developer pages (ForDevelopers, developer, auto, kBT_history, mlo_backup, implementation) moved to ecalj/MD on 2026-10-02 (not on the site)
    ]
  },
  {
    text: 'Input files (TOML, 2026-05+)',
    items: [
      { text: 'ctrlg.<sname>.toml (TOML input)', link: '/manual/toml_migration' },
      { text: 'ctrl.<sname> (legacy)', link: '/manual/lmf_legacy_ctrl' },
      { text: 'GWinput (legacy)', link: '/manual/gwinput' },
      { text: 'syml.<sname>', link: '/manual/syml' }
    ]
  },
  {
    text: 'Theory',
    items: [
      { text: 'PMT 法（APW+MTO, DFT）', link: '/theory/pmt' },
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
]

const enText: Record<string, string> = {
  '更新履歴（要約）/ What\'s new': 'What\'s new',
  'MT 半径の上限と偽の準位': 'Upper limit of MT radii and spurious levels',
  'kBT — 有限温度と均し方 (t_tetrakbt / t_sigmaw)': 'kBT — finite temperature and smearing (t_tetrakbt / t_sigmaw)',
  'MLO — 最大局在化などに代わる自動モデル化法': 'MLO — automatic model construction in place of maximal localization',
  'MLO-gwsc — Sigma を MLO 表現で内挿する QSGW（開発中）': 'MLO-gwsc — QSGW interpolating Sigma in the MLO representation (in development)',
  'PMT 法（APW+MTO, DFT）': 'PMT method (APW+MTO, DFT)',
}
function toEn(items: any[]): any[] {
  return items.map((it) => ({
    ...it,
    ...(it.text !== undefined ? { text: enText[it.text] ?? it.text } : {}),
    ...(typeof it.link === 'string' && it.link.startsWith('/') ? { link: '/en' + it.link } : {}),
    ...(it.items ? { items: toEn(it.items) } : {}),
  }))
}

export default withMermaid({
  base : '/ecaljdoc/',
  title: "ecaljdoc: ",
  vite: { publicDir: 'public' },
  description: "document of ecalj",
  cleanUrls: true,
  ignoreDeadLinks: true,
  srcExclude: ['trash/**'],
  lastUpdated: true,
  appearance: 'dark',   // dark by default; the reader can still switch (2026-11-13)
  sitemap: { hostname: 'https://ecalj.github.io/ecaljdoc/' },   // sitemap.xml for search engines (2026-11-10)
  head: [['link', { rel: 'icon', type: 'image/png', href: '/ecaljdoc/favicon.png' }]],
  // Language switch (2026-11-13): the Japanese pages at the root are the reference; English translations live under /en/
  // with the same file names (every page twice; TOOLS/ecaljdoc_en.py copies and lists stale ones). The English sidebar is
  // made from the same list: /en is put in front of the links and the Japanese labels are replaced by enText.
  locales: {
    root: { label: '日本語', lang: 'ja' },
    en: {
      label: 'English', lang: 'en', link: '/en/',
      themeConfig: {
        nav: [ { text: 'Home', link: '/en/' } ],
        sidebar: toEn(sidebarRoot),
      },
    },
  },
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
    logo: '/ecalj_logo_small.png',
    nav: [
      { text: 'Home', link: '/' },
    ],

    search: {
      provider: 'local'
    },

    sidebar: sidebarRoot,

    socialLinks: [
      { icon: 'github', link: 'https://github.com/tkotani/ecalj' },
      { icon: 'github', link: 'https://github.com/ecalj/ecaljdoc' },
    ]
  },
})
