package org.ejml.sparse;

public enum FillReducing2 {
	NONE, //
	RANDOM, //
	IDENTITY, //
	// https://github.com/PetrKryslUCSD/SymRCM.jl
	/**
	 * Symmetric Reverse Cuthill–McKee (RCM) ordering.
	 *
	 * <p>
	 * Graph-based permutation that tends to reduce matrix bandwidth/profile, often
	 * reducing fill-in and speeding up sparse Cholesky/LDLT factorizations. This
	 * variant sorts each node’s neighbors by degree, which can improve ordering
	 * quality but adds noticeable preprocessing cost on high-degree graphs.
	 * </p>
	 * <p>
	 * <b>Requirements:</b> matrix must be square and <b>structurally symmetric</b>
	 * (pattern symmetric; values do not need to be symmetric).
	 * </p>
	 */
	SYMRCM,

	/**
	 * Symmetric Reverse Cuthill–McKee (RCM) ordering without neighbor sorting.
	 *
	 * <p>
	 * Much faster preprocessing (typically near-linear in nnz) and often similar
	 * quality to {@link #SYMRCM}. Good default for very large problems when
	 * permutation time matters.
	 * </p>
	 */
	SYMRCM_NO_SORT
}