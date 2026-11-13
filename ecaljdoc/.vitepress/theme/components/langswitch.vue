<script setup>
// A plain button in the nav bar that switches between the Japanese page (the reference, at the root) and its English
// translation under /en/ (2026-11-13). Every page has an English twin with the same name (TOOLS/ecaljdoc_en.py).
import { computed } from 'vue'
import { useRoute, useData, withBase } from 'vitepress'

const route = useRoute()
const { site } = useData()
const target = computed(() => {
  const base = site.value.base                       // '/ecaljdoc/'
  let p = route.path.startsWith(base) ? route.path.slice(base.length) : route.path.replace(/^\//, '')
  p = p.replace(/\.html$/, '').replace(/\/$/, '')
  if (p === 'en' || p.startsWith('en/')) {             // English -> Japanese
    return { href: withBase('/' + p.replace(/^en\/?/, '')), label: '日本語' }
  }
  return { href: withBase('/en/' + (p === 'index' ? '' : p)), label: 'English' }
})
</script>

<template>
  <a class="ecalj-langswitch" :href="target.href">{{ target.label }}</a>
</template>

<style scoped>
.ecalj-langswitch {
  margin-left: 12px; padding: 4px 12px; border-radius: 6px; font-size: 13px; font-weight: 600; white-space: nowrap;
  border: 1px solid var(--vp-c-brand-1); color: var(--vp-c-brand-1); text-decoration: none;
}
.ecalj-langswitch:hover { background: var(--vp-c-brand-soft); }
</style>
