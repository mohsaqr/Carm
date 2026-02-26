import { defineConfig } from 'vite'
import { resolve } from 'path'

export default defineConfig({
  base: '/Carm/',
  root: resolve(__dirname),
  resolve: {
    alias: {
      'carm': resolve(__dirname, '../src/index.ts'),
    },
  },
})
