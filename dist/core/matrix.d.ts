/**
 * Matrix class for Carm.
 * Implements multiply, transpose, inverse (via Cholesky or LU), log-determinant, and SVD.
 * Pure computation — no DOM, no D3, no side effects.
 */
export declare class Matrix {
    readonly rows: number;
    readonly cols: number;
    private readonly _data;
    constructor(rows: number, cols: number, data?: readonly number[]);
    /** Build from 2-D array (row-major). */
    static fromArray(arr: readonly (readonly number[])[]): Matrix;
    /** Identity matrix of size n. */
    static identity(n: number): Matrix;
    /** Zero matrix. */
    static zeros(rows: number, cols: number): Matrix;
    /** Get element at (i, j) — 0-indexed. */
    get(i: number, j: number): number;
    /** Return 2-D array representation. */
    toArray(): number[][];
    /** Return flat row-major copy. */
    toFlat(): number[];
    /** Matrix transpose. */
    transpose(): Matrix;
    /** Matrix multiplication: this × other. */
    multiply(other: Matrix): Matrix;
    /** Scalar multiplication. */
    scale(s: number): Matrix;
    /** Element-wise add. */
    add(other: Matrix): Matrix;
    /** Element-wise subtract. */
    subtract(other: Matrix): Matrix;
    /**
     * Cholesky decomposition for symmetric positive-definite matrices.
     * Returns lower-triangular L such that this = L * L^T.
     * Throws if matrix is not SPD.
     */
    cholesky(): Matrix;
    /**
     * Log-determinant via Cholesky: log|A| = 2 * Σ log(L_ii).
     * Only valid for symmetric positive-definite matrices.
     */
    logDet(): number;
    /**
     * Inverse via LU decomposition with partial pivoting.
     * Works for any non-singular square matrix.
     * Formula: Doolittle LU, then forward/back substitution for each column of I.
     */
    inverse(): Matrix;
    /**
     * Singular Value Decomposition: A = U · S · V^T
     * Returns { U, S (diagonal values), V }.
     * Algorithm: Golub-Reinsch (one-sided Jacobi for small matrices).
     * Reference: Golub & Van Loan, "Matrix Computations", 4th ed., Algorithm 8.6.2
     */
    svd(): {
        U: Matrix;
        S: number[];
        V: Matrix;
    };
    /**
     * Pseudo-inverse via SVD: A+ = V · S^{-1} · U^T
     */
    pseudoInverse(tol?: number): Matrix;
    /** Trace (sum of diagonal elements). */
    trace(): number;
    /** Extract diagonal as array. */
    diagonal(): number[];
    /** Column vector as Matrix from array. */
    static colVec(data: readonly number[]): Matrix;
    /** Row vector as Matrix from array. */
    static rowVec(data: readonly number[]): Matrix;
    /**
     * Eigendecomposition for symmetric matrices via Jacobi iterations.
     * Returns { values, vectors } where vectors are columns of the eigenvector matrix.
     * Reference: Golub & Van Loan, Algorithm 8.4.1
     */
    eigen(): {
        values: number[];
        vectors: Matrix;
    };
}
/** Solve linear system A·x = b using the (already computed) inverse.  */
export declare function solveLinear(A: Matrix, b: readonly number[]): number[];
//# sourceMappingURL=matrix.d.ts.map