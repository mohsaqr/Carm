/**
 * Deterministic PRNG (splitmix32).
 * Shared across clustering, factor-analysis, bootstrap, and any stochastic module.
 * Default seed: 42. Every stochastic function accepts `seed?` in options.
 */
export declare class PRNG {
    private state;
    constructor(seed: number);
    next(): number;
}
//# sourceMappingURL=prng.d.ts.map