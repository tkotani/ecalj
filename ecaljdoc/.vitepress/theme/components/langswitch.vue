<script setup>
// A plain button in the nav bar that switches between the Japanese page (the reference, at the root) and its English
// translation under /en/ (2026-11-13). Pages without a translation go to the English top page.
// When a page is translated, add its path (without .md, relative to ecaljdoc/) to `translated`.
import { computed } from 'vue'
import { useRoute, useData, withBase } from 'vitepress'

const translated = ['', 'index', 'theory/pmt']

const route = useRoute()
const { site } = useData()
const target = computed(() => {
  const base = site.value.base                       // '/ecaljdoc/'
  let p = route.path.startsWith(base) ? route.path.slice(base.length) : route.path.replace(/^\//, '')
  p = p.replace(/\.html$/, '').replace(/\/$/, '')
  if (p === 'en' || p.startsWith('en/')) {             // English -> Japanese
    return { href: withBase('/' + p.replace(/^en\/?/, '')), label: '日本語' }
  }
  const has = translated.includes(p)
  return { href: withBase(has ? '/en/' + p : '/en/'), label: 'English' }
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
