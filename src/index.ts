/**
 * Carm — Statistical analysis and visualization for TypeScript/JavaScript.
 *
 * Usage:
 *   import { describe, tTestIndependent, renderViolinBox } from 'carm'
 *   import { runLMM } from 'carm/stats'
 *   import { renderCorrelogram } from 'carm/viz'
 */

// Core: types, matrix, math utilities, APA formatting
export * from './core/index.js'

// Stats: all statistical computation modules
export * from './stats/index.js'

// Viz: all D3 plot renderers
export * from './viz/index.js'
