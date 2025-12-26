package org.ejml.sparse.csc.factory;

import java.util.Arrays;
import java.util.Random;

import org.ejml.UtilEjml;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.IGrowArray;
import org.ejml.sparse.ComputePermutation;
import org.ejml.sparse.FillReducing;

/**
 * @author Peter Abeles
 */
public class FillReductionFactory_DSCC {
    public static final Random rand = new Random(234234);

    public static ComputePermutation<DMatrixSparseCSC> create( FillReducing type ) {
        switch (type) {
            case NONE:
                return null;

            case RANDOM:
                return new ComputePermutation<>(true, true) {
                    @Override
                    @SuppressWarnings("NullAway") // constructor parameters ensures these are not null
                    public void process( DMatrixSparseCSC m ) {
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

            case IDENTITY:
                return new ComputePermutation<>(true, true) {
                    @Override
                    @SuppressWarnings("NullAway") // constructor parameters ensures these are not null
                    public void process( DMatrixSparseCSC m ) {
                        prow.reshape(m.numRows);
                        pcol.reshape(m.numCols);
                        fillSequence(prow);
                        fillSequence(pcol);
                    }
                };

            case SYMRCM:
                return new ComputePermutation<>(true, true) {
                    @Override
                    @SuppressWarnings("NullAway") // constructor parameters ensures these are not null
                    public void process( DMatrixSparseCSC m ) {
                        if (m.numRows != m.numCols) {
							throw new IllegalArgumentException("SymRCM requires a square (structurally symmetric) matrix");
						}
                        prow.reshape(m.numRows);
                        pcol.reshape(m.numCols);
                        symrcm(m, true, prow);
                        System.arraycopy(prow.data, 0, pcol.data, 0, prow.length);
                    }
                };

            case SYMRCM_NO_SORT:
                return new ComputePermutation<>(true, true) {
                    @Override
                    @SuppressWarnings("NullAway") // constructor parameters ensures these are not null
                    public void process( DMatrixSparseCSC m ) {
                        if (m.numRows != m.numCols) {
							throw new IllegalArgumentException("SymRCM requires a square (structurally symmetric) matrix");
						}
                        prow.reshape(m.numRows);
                        pcol.reshape(m.numCols);
                        symrcm(m, false, prow);
                        System.arraycopy(prow.data, 0, pcol.data, 0, prow.length);
                    }
                };

            default:
                throw new RuntimeException("Unknown " + type);
        }
    }

    private static void fillSequence( IGrowArray perm ) {
        for (int i = 0; i < perm.length; i++) {
            perm.data[i] = i;
        }
    }

    /**
     * Reverse Cuthill-McKee permutation for a (structurally) symmetric sparse matrix.
     */
    private static void symrcm( DMatrixSparseCSC A, boolean sortByDegree, IGrowArray perm ) {
        final int n = A.numCols;
        perm.reshape(n);
        if (n == 0) {
			return;
		}

        // degree[v] = number of neighbors in column v (matches the Julia code: colptr[v+1]-colptr[v])
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

    private static int bfsCollectComponent( int root, boolean[] label, int[] queue, int[] component,
                                           int[] colIdx, int[] nzRows ) {
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

    private static int bfsAppendOrder( int root, boolean[] label, int[] queue,
                                       int[] order, int orderLen,
                                       int[] colIdx, int[] nzRows ) {
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

    private static int findMinDegreeVertex( int[] vertices, int length, int[] degree ) {
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

    private static void sortAllColumnsByNeighborDegree( int n, int[] colIdx, int[] nzRows, int[] degree ) {
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
    private static void quickSortByDegree( int[] a, int lo, int hi, int[] degree ) {
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

    private static boolean lessByDegreeThenIndex( int x, int y, int[] degree, int degreeY ) {
        int dx = degree[x];
        if (dx != degreeY) {
			return dx < degreeY;
		}
        return x < y;
    }
}