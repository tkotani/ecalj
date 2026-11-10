import DefaultTheme from 'vitepress/theme'
import { h } from 'vue'
import { withBase } from 'vitepress'
import './custom.css'
import FileList from './components/filelist.vue'
import { setupLightbox } from './lightbox.js'

// The front-page logo with the origin of the name under it (2026-11-10):
// ecalj = electron, calculation, japan; the letters of the name in a slightly heavier weight.
const heroImage = () => h('div', { class: 'ecalj-hero' }, [
  h('img', { src: withBase('/ecalj_logo.png'), alt: 'ecalj' }),
  h('p', { class: 'ecalj-origin' }, [
    h('b', 'e'), 'lectron · ', h('b', 'cal'), 'culation · ', h('b', 'j'), 'apan',
  ]),
])

export default {
 ...DefaultTheme,
 Layout() {
   return h(DefaultTheme.Layout, null, { 'home-hero-image': heroImage })
 },
 enhanceApp({ app }) {
   app.component('filelist', FileList)
   setupLightbox()
 }
}
