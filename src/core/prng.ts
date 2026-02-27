/**
 * Deterministic PRNG (splitmix32).
 * Shared across clustering, factor-analysis, bootstrap, and any stochastic module.
 * Default seed: 42. Every stochastic function accepts `seed?` in options.
 */

export class PRNG {
  private state: number
  constructor(seed: number) { this.state = seed >>> 0 }
  next(): number {
    this.state = (this.state + 0x9E3779B9) | 0
    let t = this.state ^ (this.state >>> 16)
    t = Math.imul(t, 0x21F0AAAD)
    t = t ^ (t >>> 15)
    t = Math.imul(t, 0x735A2D97)
    t = t ^ (t >>> 15)
    return (t >>> 0) / 4294967296
  }
}
