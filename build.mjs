import { build } from 'esbuild'
import { execSync } from 'child_process'
import { rmSync } from 'fs'

const entry = [
  'src/index.ts',
  'src/core/index.ts',
  'src/stats/index.ts',
  'src/viz/index.ts',
]

// Clean
rmSync('dist', { recursive: true, force: true })

// ESM (with splitting for shared chunks)
await build({
  entryPoints: entry,
  bundle: true,
  format: 'esm',
  outdir: 'dist',
  outbase: 'src',
  splitting: true,
  sourcemap: true,
  external: ['d3'],
  treeShaking: true,
})

// CJS
await build({
  entryPoints: entry,
  bundle: true,
  format: 'cjs',
  outdir: 'dist',
  outbase: 'src',
  sourcemap: true,
  external: ['d3'],
  treeShaking: true,
  outExtension: { '.js': '.cjs' },
})

// DTS
execSync('npx tsc --emitDeclarationOnly --declaration --declarationMap', { stdio: 'inherit' })

console.log('Build complete: ESM + CJS + DTS')
