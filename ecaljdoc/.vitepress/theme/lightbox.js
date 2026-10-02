/* ---------------------------------------------------------------------------
   Click-to-enlarge for figures in the docs.

   The doc column is ~688px wide, which is fine for a single band plot but
   useless for a multi-panel figure (the FeMgO comparison is 4 columns x 2 rows
   at 3700px natural width). Clicking any image in .vp-doc opens it over the
   page: first click fits it to the viewport, a second click switches to 1:1 so
   you can pan around a detail. Esc or a click on the backdrop closes it.

   No dependency: one delegated listener, added once on the client.
--------------------------------------------------------------------------- */
export function setupLightbox() {
  if (typeof window === 'undefined') return
  if (window.__ecaljLightbox) return
  window.__ecaljLightbox = true

  let overlay = null

  const close = () => {
    if (!overlay) return
    overlay.remove()
    overlay = null
    document.body.style.overflow = ''
    window.removeEventListener('keydown', onKey)
  }

  const onKey = (e) => {
    if (e.key === 'Escape') close()
  }

  const open = (src, alt) => {
    close()
    overlay = document.createElement('div')
    overlay.className = 'ecalj-lightbox'
    overlay.innerHTML = `
      <button class="ecalj-lightbox-close" aria-label="close">&times;</button>
      <div class="ecalj-lightbox-hint">クリックで原寸 / Esc で閉じる</div>
      <div class="ecalj-lightbox-scroll"><img src="${src}" alt="${alt || ''}"></div>`
    const img = overlay.querySelector('img')

    // first state: fitted to the viewport; clicking the image goes to 1:1
    let full = false
    img.addEventListener('click', (e) => {
      e.stopPropagation()
      full = !full
      img.classList.toggle('full', full)
      overlay.querySelector('.ecalj-lightbox-hint').textContent =
        full ? 'クリックで全体表示 / Esc で閉じる' : 'クリックで原寸 / Esc で閉じる'
    })
    overlay.addEventListener('click', close)
    overlay.querySelector('.ecalj-lightbox-close')
           .addEventListener('click', (e) => { e.stopPropagation(); close() })

    document.body.appendChild(overlay)
    document.body.style.overflow = 'hidden'
    window.addEventListener('keydown', onKey)
  }

  document.addEventListener('click', (e) => {
    const img = e.target.closest('.vp-doc img')
    if (!img) return
    // leave real links (e.g. a badge wrapped in <a>) alone
    if (img.closest('a')) return
    e.preventDefault()
    open(img.currentSrc || img.src, img.alt)
  })
}
