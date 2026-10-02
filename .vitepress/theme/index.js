import DefaultTheme from 'vitepress/theme'
import './custom.css'
import FileList from './components/filelist.vue'
import { setupLightbox } from './lightbox.js'

export default {
 ...DefaultTheme,
 enhanceApp({ app }) {
   app.component('filelist', FileList)
   setupLightbox()
 }
}
