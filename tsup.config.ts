import { defineConfig } from 'tsup'

export default defineConfig({
  entry: {
    index:          'src/index.ts',
    'core/index':   'src/core/index.ts',
    'stats/index':  'src/stats/index.ts',
    'viz/index':    'src/viz/index.ts',
  },
  format: ['esm', 'cjs'],
  dts: true,
  sourcemap: true,
  clean: true,
  splitting: true,
  treeshake: true,
  external: ['d3'],
})
