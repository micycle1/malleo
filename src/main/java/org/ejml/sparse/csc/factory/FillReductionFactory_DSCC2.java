package org.ejml.sparse.csc.factory;

import java.util.Arrays;
import java.util.Random;

import org.ejml.UtilEjml;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.IGrowArray;
import org.ejml.sparse.ComputePermutation;
import org.ejml.sparse.FillReducing2;

/**
 * @author Peter Abeles
 */
public class FillReductionFactory_DSCC2 {
	public static final Random rand = new Random(234234);

	public static ComputePermutation<DMatrixSparseCSC> create(FillReducing2 type) {
		switch (type) {
			case NONE :
				return null;

			case RANDOM :
				return new ComputePermutation<>(true, true) {
					@Override
					@SuppressWarnings("NullAway") // constructor parameters ensures these are not null
					public void process(DMatrixSparseCSC m) {
						prow.reshape(m.numRows);
						pcol.reshape(m.numCols);
						fillSequence(prow);
						fillSequence(pcol);
						Random _rand;
						synchronized (rand) {
							_rand = new Random(rand.nextInt());
						}
						UtilEjml.shuffle(prow.data, prow.length, 0, prow.length, _rand);
						UtilEjml.shuffle(pcol.data, pcol.length, 0, pcol.length, _rand);
					}
				};

			case IDENTITY :
				return new ComputePermutation<>(true, true) {
					@Override
					@SuppressWarnings("NullAway") // constructor parameters ensures these are not null
					public void process(DMatrixSparseCSC m) {
						prow.reshape(m.numRows);
						pcol.reshape(m.numCols);
						fillSequence(prow);
						fillSequence(pcol);
					}
				};

			case SYMRCM :
				return new ComputePermutation<>(true, true) {
					@Override
					@SuppressWarnings("NullAway") // constructor parameters ensures these are not null
					public void process(DMatrixSparseCSC m) {
						if (m.numRows != m.numCols) {
							throw new IllegalArgumentException("SymRCM requires a square (structurally symmetric) matrix");
						}
						prow.reshape(m.numRows);
						pcol.reshape(m.numCols);
						symrcm(m, true, prow);
						System.arraycopy(prow.data, 0, pcol.data, 0, prow.length);
					}
				};

			case SYMRCM_NO_SORT :
				return new ComputePermutation<>(true, true) {
					@Override
					@SuppressWarnings("NullAway") // constructor parameters ensures these are not null
					public void process(DMatrixSparseCSC m) {
						if (m.numRows != m.numCols) {
							throw new IllegalArgumentException("SymRCM requires a square (structurally symmetric) matrix");
						}
						prow.reshape(m.numRows);
						pcol.reshape(m.numCols);
						symrcm(m, false, prow);
						System.arraycopy(prow.data, 0, pcol.data, 0, prow.length);
					}
				};
			case SYMRCM_RCMPLUSPLUS :
				return new ComputePermutation<>(true, true) {
					@Override
					@SuppressWarnings("NullAway") // constructor parameters ensures these are not null
					public void process(DMatrixSparseCSC m) {
						if (m.numRows != m.numCols) {
							throw new IllegalArgumentException("SymRCM requires a square (structurally symmetric) matrix");
						}
						prow.reshape(m.numRows);
						pcol.reshape(m.numCols);
						symrcmRcmpp(m, false, prow);
						System.arraycopy(prow.data, 0, pcol.data, 0, prow.length);
					}
				};

			default :
				throw new RuntimeException("Unknown " + type);
		}
	}

	private static void fillSequence(IGrowArray perm) {
		for (int i = 0; i < perm.length; i++) {
			perm.data[i] = i;
		}
	}

	/**
	 * Reverse Cuthill-McKee permutation for a (structurally) symmetric sparse
	 * matrix.
	 */
	private static void symrcm(DMatrixSparseCSC A, boolean sortByDegree, IGrowArray perm) {
		final int n = A.numCols;
		perm.reshape(n);
		if (n == 0) {
			return;
		}

		// degree[v] = number of neighbors in column v (matches the Julia code:
		// colptr[v+1]-colptr[v])
		final int[] degree = new int[n];
		for (int v = 0; v < n; v++) {
			degree[v] = A.col_idx[v + 1] - A.col_idx[v];
		}

		// Neighbor access (optionally sort each neighbor list by neighbor degree)
		final int[] colIdx;
		final int[] nzRows;
		if (sortByDegree) {
			colIdx = A.col_idx.clone();
			nzRows = Arrays.copyOf(A.nz_rows, A.nz_length);
			sortAllColumnsByNeighborDegree(n, colIdx, nzRows, degree);
		} else {
			colIdx = A.col_idx;
			nzRows = A.nz_rows;
		}

		final boolean[] label = new boolean[n];
		final int[] queue = new int[n];
		final int[] component = new int[n];
		final int[] order = perm.data;

		int orderLen = 0;

		// Process each connected component:
		// 1) BFS to get component vertices
		// 2) choose root = min-degree vertex in that component
		// 3) reset labels for the component
		// 4) BFS from root to append Cuthill-McKee ordering for that component
		for (int start = 0; start < n; start++) {
			if (label[start]) {
				continue;
			}

			int compLen = bfsCollectComponent(start, label, queue, component, colIdx, nzRows);
			int root = findMinDegreeVertex(component, compLen, degree);

			for (int i = 0; i < compLen; i++) {
				label[component[i]] = false;
			}

			orderLen = bfsAppendOrder(root, label, queue, order, orderLen, colIdx, nzRows);
		}

		// Reverse for RCM
		for (int i = 0, j = n - 1; i < j; i++, j--) {
			int tmp = order[i];
			order[i] = order[j];
			order[j] = tmp;
		}

		// If input graph has isolated vertices, orderLen should still be n
		// (kept as a sanity check; can be removed if undesired)
		if (orderLen != n) {
			throw new IllegalStateException("symrcm: internal error. orderLen=" + orderLen + " n=" + n);
		}
	}

	private static int bfsCollectComponent(int root, boolean[] label, int[] queue, int[] component, int[] colIdx, int[] nzRows) {
		int head = 0, tail = 0, compLen = 0;

		label[root] = true;
		queue[tail++] = root;

		while (head < tail) {
			int u = queue[head++];
			component[compLen++] = u;

			for (int p = colIdx[u]; p < colIdx[u + 1]; p++) {
				int v = nzRows[p];
				if (!label[v]) {
					label[v] = true;
					queue[tail++] = v;
				}
			}
		}
		return compLen;
	}

	private static int bfsAppendOrder(int root, boolean[] label, int[] queue, int[] order, int orderLen, int[] colIdx, int[] nzRows) {
		int head = 0, tail = 0;

		label[root] = true;
		queue[tail++] = root;

		while (head < tail) {
			int u = queue[head++];
			order[orderLen++] = u;

			for (int p = colIdx[u]; p < colIdx[u + 1]; p++) {
				int v = nzRows[p];
				if (!label[v]) {
					label[v] = true;
					queue[tail++] = v;
				}
			}
		}
		return orderLen;
	}

	private static int findMinDegreeVertex(int[] vertices, int length, int[] degree) {
		int best = vertices[0];
		int bestDeg = degree[best];

		for (int i = 1; i < length; i++) {
			int v = vertices[i];
			int d = degree[v];
			if (d < bestDeg || (d == bestDeg && v < best)) {
				best = v;
				bestDeg = d;
			}
		}
		return best;
	}

	private static void sortAllColumnsByNeighborDegree(int n, int[] colIdx, int[] nzRows, int[] degree) {
		for (int col = 0; col < n; col++) {
			int start = colIdx[col];
			int end = colIdx[col + 1] - 1; // inclusive
			if (end <= start) {
				continue;
			}
			quickSortByDegree(nzRows, start, end, degree);
		}
	}

	// Sort nzRows[lo..hi] by (degree[nzRows[i]], nzRows[i])
	private static void quickSortByDegree(int[] a, int lo, int hi, int[] degree) {
		while (lo < hi) {
			int i = lo, j = hi;
			int pivot = a[lo + ((hi - lo) >>> 1)];
			int pivotDeg = degree[pivot];

			while (i <= j) {
				while (lessByDegreeThenIndex(a[i], pivot, degree, pivotDeg)) {
					i++;
				}
				while (lessByDegreeThenIndex(pivot, a[j], degree, pivotDeg)) {
					j--;
				}
				if (i <= j) {
					int tmp = a[i];
					a[i] = a[j];
					a[j] = tmp;
					i++;
					j--;
				}
			}

			// Recurse on smaller partition first (limits recursion depth)
			if (j - lo < hi - i) {
				if (lo < j) {
					quickSortByDegree(a, lo, j, degree);
				}
				lo = i;
			} else {
				if (i < hi) {
					quickSortByDegree(a, i, hi, degree);
				}
				hi = j;
			}
		}
	}

	private static boolean lessByDegreeThenIndex(int x, int y, int[] degree, int degreeY) {
		int dx = degree[x];
		if (dx != degreeY) {
			return dx < degreeY;
		}
		return x < y;
	}

	/**
	 * RCM++ (BNF-based start-node selection) + Reverse Cuthill-McKee ordering for a
	 * (structurally) symmetric sparse matrix.
	 *
	 * BNF heuristic is adapted from the provided Python code: it iteratively
	 * increases eccentricity (pseudo-peripheral search) while tracking the minimum
	 * "width" (max BFS layer size) seen during the run.
	 */
	private static void symrcmRcmpp(DMatrixSparseCSC A, boolean sortByDegree, IGrowArray perm) {
		final int n = A.numCols;
		perm.reshape(n);
		if (n == 0)
			return;

		// degree[v] = number of neighbors in column v
		final int[] degree = new int[n];
		for (int v = 0; v < n; v++) {
			degree[v] = A.col_idx[v + 1] - A.col_idx[v];
		}

		// Neighbor access (optionally sort each neighbor list by neighbor degree)
		final int[] colIdx;
		final int[] nzRows;
		if (sortByDegree) {
			colIdx = A.col_idx.clone();
			nzRows = Arrays.copyOf(A.nz_rows, A.nz_length);
			sortAllColumnsByNeighborDegree(n, colIdx, nzRows, degree);
		} else {
			colIdx = A.col_idx;
			nzRows = A.nz_rows;
		}

		final boolean[] label = new boolean[n];
		final int[] queue = new int[n];
		final int[] order = perm.data;

		// Stamping array for internal BFS runs inside BNF (avoid Arrays.fill every
		// time)
		final int[] seen = new int[n];
		final int[] tokenRef = new int[] { 1 };
		final BfsLevelsResult bfsRes = new BfsLevelsResult();

		int orderLen = 0;

		// Process each connected component
		for (int start = 0; start < n; start++) {
			if (label[start])
				continue;

			// Pick RCM++ root using BNF, seeded by the first unvisited node in the
			// component
			int root = bnfFindRoot(start, label, colIdx, nzRows, degree, queue, seen, tokenRef, bfsRes);

			// Standard CM BFS from the chosen root to append ordering for this component
			orderLen = bfsAppendOrder(root, label, queue, order, orderLen, colIdx, nzRows);
		}

		// Reverse for RCM
		for (int i = 0, j = n - 1; i < j; i++, j--) {
			int tmp = order[i];
			order[i] = order[j];
			order[j] = tmp;
		}

		if (orderLen != n) {
			throw new IllegalStateException("symrcmRcmpp: internal error. orderLen=" + orderLen + " n=" + n);
		}
	}

	private static class BfsLevelsResult {
		int eccentricity;
		int width;
		int farthestStart;
		int farthestEnd;
	}

	/**
	 * Bi-Criteria Node Finder (BNF) from the Python reference: - run BFS from v,
	 * compute eccentricity l and width w (max layer size) - track node a with
	 * minimal width encountered - if eccentricity stops increasing, return a -
	 * otherwise set v to min-degree node among farthest-layer nodes and continue
	 */
	private static int bnfFindRoot(int start, boolean[] label, int[] colIdx, int[] nzRows, int[] degree, int[] queue, int[] seen, int[] tokenRef,
			BfsLevelsResult bfsRes) {
		int lp = 0;
		int v = start;

		int bestWidth = Integer.MAX_VALUE;
		int a = v;

		while (true) {
			tokenRef[0] = bfsLevels(v, label, colIdx, nzRows, seen, tokenRef[0], queue, bfsRes);

			int l = bfsRes.eccentricity;
			int w = bfsRes.width;

			if (w <= bestWidth) {
				bestWidth = w;
				a = v;
			}

			if (l <= lp) {
				break;
			}
			lp = l;

			// next v = min-degree vertex among farthest nodes (tie by index)
			v = findMinDegreeInQueueSegment(queue, bfsRes.farthestStart, bfsRes.farthestEnd, degree);
		}

		return a;
	}

	/**
	 * BFS that computes: - eccentricity = max distance from root in the (unlabeled)
	 * reachable subgraph - width = maximum BFS layer size - farthest layer segment
	 * [farthestStart, farthestEnd) stored inside 'queue'
	 *
	 * Uses 'seen' stamping with token to avoid clearing arrays each call.
	 */
	private static int bfsLevels(int root, boolean[] label, int[] colIdx, int[] nzRows, int[] seen, int token, int[] queue, BfsLevelsResult out) {
		// Advance token, reset if needed
		if (token == Integer.MAX_VALUE) {
			Arrays.fill(seen, 0);
			token = 1;
		} else {
			token++;
		}

		int tail = 0;
		queue[tail++] = root;
		seen[root] = token;

		int levelStart = 0;
		int levelEnd = tail;

		int distance = -1;
		int width = 0;

		int lastLevelStart = 0;
		int lastLevelEnd = 0;

		while (levelStart < tail) {
			distance++;

			int levelSize = levelEnd - levelStart;
			if (levelSize > width)
				width = levelSize;

			lastLevelStart = levelStart;
			lastLevelEnd = levelEnd;

			for (int idx = levelStart; idx < levelEnd; idx++) {
				int u = queue[idx];

				for (int p = colIdx[u]; p < colIdx[u + 1]; p++) {
					int v = nzRows[p];

					// Don't walk into already-ordered vertices (robustness across components)
					if (label[v])
						continue;

					if (seen[v] != token) {
						seen[v] = token;
						queue[tail++] = v;
					}
				}
			}

			levelStart = levelEnd;
			levelEnd = tail;
		}

		out.eccentricity = distance;
		out.width = width;
		out.farthestStart = lastLevelStart;
		out.farthestEnd = lastLevelEnd;
		return token;
	}

	private static int findMinDegreeInQueueSegment(int[] queue, int start, int end, int[] degree) {
		int best = queue[start];
		int bestDeg = degree[best];

		for (int i = start + 1; i < end; i++) {
			int v = queue[i];
			int d = degree[v];
			if (d < bestDeg || (d == bestDeg && v < best)) {
				best = v;
				bestDeg = d;
			}
		}
		return best;
	}
}