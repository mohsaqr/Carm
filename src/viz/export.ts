/**
 * SVG + PNG export for Carm plots.
 * Exports at publication quality (300 DPI for PNG).
 */

/** Download SVG source of a plot container. */
export function exportSVG(container: HTMLElement, filename = 'carm-plot.svg'): void {
  const svgEl = container.querySelector('svg')
  if (!svgEl) throw new Error('exportSVG: no SVG element found in container')

  // Add xmlns if missing
  if (!svgEl.hasAttribute('xmlns')) {
    svgEl.setAttribute('xmlns', 'http://www.w3.org/2000/svg')
  }

  const serializer = new XMLSerializer()
  const svgStr = serializer.serializeToString(svgEl)
  const blob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' })
  triggerDownload(URL.createObjectURL(blob), filename)
}

/** Export plot container as PNG at specified DPI (default 300). */
export function exportPNG(
  container: HTMLElement,
  filename = 'carm-plot.png',
  dpi = 300
): Promise<void> {
  const svgEl = container.querySelector('svg')
  if (!svgEl) throw new Error('exportPNG: no SVG element found in container')

  const svgW = svgEl.viewBox.baseVal.width || svgEl.clientWidth
  const svgH = svgEl.viewBox.baseVal.height || svgEl.clientHeight
  const scale = dpi / 96  // 96 = screen DPI

  const canvas = document.createElement('canvas')
  canvas.width = svgW * scale
  canvas.height = svgH * scale
  const ctx = canvas.getContext('2d')!
  ctx.scale(scale, scale)

  if (!svgEl.hasAttribute('xmlns')) svgEl.setAttribute('xmlns', 'http://www.w3.org/2000/svg')
  const svgStr = new XMLSerializer().serializeToString(svgEl)
  const blob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' })
  const url = URL.createObjectURL(blob)

  return new Promise((resolve, reject) => {
    const img = new Image()
    img.onload = () => {
      ctx.drawImage(img, 0, 0)
      URL.revokeObjectURL(url)
      canvas.toBlob(pngBlob => {
        if (!pngBlob) { reject(new Error('Canvas toBlob failed')); return }
        triggerDownload(URL.createObjectURL(pngBlob), filename)
        resolve()
      }, 'image/png')
    }
    img.onerror = () => { URL.revokeObjectURL(url); reject(new Error('SVG to Image conversion failed')) }
    img.src = url
  })
}

function triggerDownload(url: string, filename: string): void {
  const a = document.createElement('a')
  a.href = url
  a.download = filename
  document.body.appendChild(a)
  a.click()
  document.body.removeChild(a)
  setTimeout(() => URL.revokeObjectURL(url), 1000)
}
