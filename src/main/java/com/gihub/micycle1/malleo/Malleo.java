package com.gihub.micycle1.malleo;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import org.ejml.data.DMatrixRMaj;
import org.ejml.data.DMatrixSparseCSC;
import org.ejml.data.DMatrixSparseTriplet;
import org.ejml.interfaces.linsol.LinearSolverSparse;
import org.ejml.ops.DConvertMatrixStruct;
import org.ejml.sparse.FillReducing;
import org.ejml.sparse.csc.CommonOps_DSCC;
import org.ejml.sparse.csc.factory.LinearSolverFactory_DSCC;
import org.locationtech.jts.algorithm.locate.IndexedPointInAreaLocator;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.Location;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.index.SpatialIndex;
import org.locationtech.jts.index.hprtree.HPRtree;
import org.locationtech.jts.operation.overlayng.CoverageUnion;

/**
 * 2D As-Rigid-As-Possible (ARAP) shape deformation using the edge-based 2-step
 * method from Igarashi & Igarashi, JGT 2009.
 *
 * <p>
 * Workflow:
 * <ol>
 * <li><b>Construction</b> (this class): triangulates the input polygon
 * (optionally refined with Steiner points), extracts edges, and precomputes
 * per-edge terms used by ARAP.</li>
 * <li><b>{@link #compile(List)}</b>: "compilation" step for a chosen handle
 * set. Locates each handle in the mesh, computes barycentric coordinates,
 * assembles the two normal-equation matrices (step1/step2), and factorizes them
 * for fast repeated solves.</li>
 * <li><b>{@link #deform(Compiled, List)}</b>: performs one interactive update
 * given handle target positions.</li>
 * </ol>
 *
 * <p>
 * Input/output:
 * <ul>
 * <li>Handles may be placed at arbitrary positions inside the polygon (not just
 * at mesh vertices) using barycentric constraints.</li>
 * <li>The returned polygon is produced by moving only the original boundary
 * ring vertices (shell + holes); interior triangulation vertices are used for
 * computation only.</li>
 * </ul>
 *
 * <p>
 * This implementation solves sparse normal equations with Cholesky
 * factorization (SPD expected).
 * 
 * @author Michael Carleton
 */
public final class Malleo {

	private final Params params;
	private final Mesh mesh;

	private final DMatrixSparseCSC gram1Base; // edges only, 2V×2V
	private final DMatrixSparseCSC gram2Base; // edges only, V×V

	/**
	 * Creates a deformer for a triangulation.
	 *
	 * @param restPolygon rest shape; may contain holes
	 */
	public Malleo(Geometry triangulation) {
		this((Polygon) CoverageUnion.union(triangulation), triangulation, new Params());
	}

	/**
	 * Creates a deformer for a triangulation.
	 *
	 * <p>
	 * Subsequent calls to {@link #compile(List)} define a particular handle set for
	 * interactive deformation.
	 *
	 * @param restPolygon rest shape; may contain holes
	 * @param params      algorithm and triangulation parameters
	 * @throws NullPointerException if any argument is null
	 */
	public Malleo(Geometry triangulation, Params params) {
		this((Polygon) CoverageUnion.union(triangulation), triangulation, params);
	}

	/**
	 * Creates a deformer from an externally-provided triangulation geometry.
	 *
	 * The provided {@code triangulation} Geometry is expected to contain triangle
	 * polygons (e.g. the result of a triangulation builder). Triangles whose
	 * centroids lie outside {@code restPolygon} are ignored (same behavior as the
	 * internal triangulation path).
	 *
	 * @param restPolygon   rest shape; may contain holes
	 * @param triangulation triangulation geometry (collection of triangle polygons)
	 * @param params        algorithm and triangulation parameters
	 * @throws NullPointerException if any argument is null
	 */
	private Malleo(Polygon restPolygon, Geometry triangulation, Params params) {
		this.params = Objects.requireNonNull(params);
		Objects.requireNonNull(restPolygon, "restPolygon must not be null");
		Objects.requireNonNull(triangulation, "triangulation must not be null");
		this.mesh = buildMesh(restPolygon, triangulation, params);
		int V = mesh.vertexCount();
		this.gram1Base = assembleGramStep1EdgesOnly(V);
		this.gram2Base = assembleGramStep2EdgesOnly(V);
	}

	/**
	 * Prepares an interactive deformation session for a fixed handle set.
	 *
	 * <p>
	 * This computes, for each handle rest position, the containing triangle and
	 * barycentric weights, then assembles and factorizes the two normal-equation
	 * matrices used by ARAP:
	 * <ul>
	 * <li><b>Step 1</b>: similarity solve (rotation + uniform scale)</li>
	 * <li><b>Step 2</b>: scale-adjustment solve (uses normalized per-edge rotation
	 * from step 1)</li>
	 * </ul>
	 *
	 * @param handleRestPoints handle locations in rest-space; each must lie inside
	 *                         (or on) {@code restPolygon}
	 * @return compiled, reusable factorization state for fast calls to
	 *         {@link #deform(Compiled, List)}
	 * @throws IllegalArgumentException if fewer than 2 handles are provided, or if
	 *                                  any handle lies outside the polygon
	 * @throws IllegalStateException    if factorization fails (e.g., matrix not SPD
	 *                                  due to degeneracy)
	 * @throws NullPointerException     if {@code handleRestPoints} is null
	 */
	public Compiled compile(List<Coordinate> handleRestPoints) {
		Objects.requireNonNull(handleRestPoints);
		if (handleRestPoints.size() < 2) {
			throw new IllegalArgumentException("ARAP requires at least 2 handles for a stable solve.");
		}

		// barycentric handles
		List<HandleInfo> handles = new ArrayList<>(handleRestPoints.size());
		for (Coordinate p : handleRestPoints) {
			if (mesh.locator.locate(p) == Location.EXTERIOR) {
//			var point = mesh.gf.createPoint(p);
//				p = DistanceOp.nearestPoints(mesh.restPolygon, point)[0];
				throw new IllegalArgumentException("Handle rest point must be inside polygon: " + p);
			}
			handles.add(locateHandle(p));
		}

		int V = mesh.vertexCount();
		double w2 = params.handleWeight * params.handleWeight;

		// build sparse handle term for step1 (2V×2V)
		DMatrixSparseTriplet t1 = new DMatrixSparseTriplet(2 * V, 2 * V, handles.size() * 36);
		for (HandleInfo h : handles) {
			int ax = 2 * h.a, ay = ax + 1;
			int bx = 2 * h.b, by = bx + 1;
			int cx = 2 * h.c, cy = cx + 1;
			addBaryConstraintGram1D(t1, ax, bx, cx, h.wa, h.wb, h.wc, w2);
			addBaryConstraintGram1D(t1, ay, by, cy, h.wa, h.wb, h.wc, w2);
		}
		DMatrixSparseCSC handle1 = toCsc(t1);

		// build sparse handle term for step2 (V×V)
		DMatrixSparseTriplet t2 = new DMatrixSparseTriplet(V, V, handles.size() * 12);
		for (HandleInfo h : handles) {
			addBaryConstraintGram2D(t2, h.a, h.b, h.c, h.wa, h.wb, h.wc, w2);
		}
		DMatrixSparseCSC handle2 = toCsc(t2);

		// add handle terms to the precomputed bases
		DMatrixSparseCSC gram1 = new DMatrixSparseCSC(gram1Base.numRows, gram1Base.numCols, gram1Base.nz_length + handle1.nz_length);
		CommonOps_DSCC.add(1.0, gram1Base, 1.0, handle1, gram1, null, null);
		CommonOps_DSCC.duplicatesAdd(gram1, null);

		DMatrixSparseCSC gram2 = new DMatrixSparseCSC(gram2Base.numRows, gram2Base.numCols, gram2Base.nz_length + handle2.nz_length);
		CommonOps_DSCC.add(1.0, gram2Base, 1.0, handle2, gram2, null, null);
		CommonOps_DSCC.duplicatesAdd(gram2, null);

		LinearSolverSparse<DMatrixSparseCSC, DMatrixRMaj> solver1 = LinearSolverFactory_DSCC.cholesky(FillReducing.SYMRCM_NO_SORT);
		if (!solver1.setA(gram1)) {
			throw new IllegalStateException("Step1 factorization failed (matrix not SPD?).");
		}
		LinearSolverSparse<DMatrixSparseCSC, DMatrixRMaj> solver2 = LinearSolverFactory_DSCC.cholesky(FillReducing.SYMRCM_NO_SORT);
		if (!solver2.setA(gram2)) {
			throw new IllegalStateException("Step2 factorization failed (matrix not SPD?).");
		}

		return new Compiled(handles, solver1, solver2, V);
	}

	/**
	 * Assembles the Step 1 normal matrix (2V x 2V): edge similarity energy + handle
	 * constraints.
	 */
	private DMatrixSparseCSC assembleGramStep1EdgesOnly(int V) {
		int N = 2 * V;
		DMatrixSparseTriplet T = new DMatrixSparseTriplet(N, N, mesh.edges.size() * 80);
		for (Edge e : mesh.edges) {
			int m = e.boundary ? 3 : 4;
			int cols = 2 * m;
			int[] g = new int[cols];
			for (int t = 0; t < m; t++) {
				g[2 * t] = 2 * e.neigh[t];
				g[2 * t + 1] = 2 * e.neigh[t] + 1;
			}
			double we = e.weight;
			double[][] hk = e.hk;
			for (int a = 0; a < cols; a++) {
				for (int b = a; b < cols; b++) {
					double v = (hk[0][a] * hk[0][b] + hk[1][a] * hk[1][b]) * we;
					if (Math.abs(v) < 1e-30)
						continue;
					addSym(T, g[a], g[b], v);
				}
			}
		}
		return toCsc(T);
	}

	/**
	 * Assembles the Step 2 normal matrix (V x V): weighted Laplacian + handle
	 * constraints.
	 */
	private DMatrixSparseCSC assembleGramStep2EdgesOnly(int V) {
		DMatrixSparseTriplet T = new DMatrixSparseTriplet(V, V, mesh.edges.size() * 6);
		for (Edge e : mesh.edges) {
			double we = e.weight;
			int i = e.vi, j = e.vj;
			addSym(T, i, i, we);
			addSym(T, j, j, we);
			addSym(T, i, j, -we);
		}
		return toCsc(T);
	}

	/**
	 * Performs one deformation update for the given handle target positions.
	 *
	 * <p>
	 * The handle list must match the handle rest points used during
	 * {@link #compile(List)} (same count and order). This method is intended to be
	 * called repeatedly (e.g., during dragging), reusing the sparse factorizations
	 * stored in {@link Compiled}.
	 * 
	 *
	 * @param compiled      result of {@link #compile(List)} for the desired handle
	 *                      set
	 * @param handleTargets target positions (world/deformed space) for each handle,
	 *                      same order as compilation
	 * @return deformed polygon (shell and holes), with boundary vertices updated.
	 *         Can self-intersect.
	 * @throws IllegalArgumentException if the handle count/order does not match
	 *                                  {@code compiled}
	 * @throws NullPointerException     if any argument is null
	 */
	public Polygon deform(Compiled compiled, List<Coordinate> handleTargets) {
		Objects.requireNonNull(compiled);
		Objects.requireNonNull(handleTargets);
		if (handleTargets.size() != compiled.handles.size()) {
			throw new IllegalArgumentException("handleTargets.size must equal compiled handle count.");
		}

		int V = compiled.vertexCount;
		double w = params.handleWeight;
		double w2 = w * w;

		// ---------- Step 1: similarity solve ----------
		// rhs1 = C1^T b1, where C1 row coeffs are w*bary, and b1 entries are w*target
		DMatrixRMaj rhs1 = new DMatrixRMaj(2 * V, 1);
		for (int hi = 0; hi < compiled.handles.size(); hi++) {
			HandleInfo h = compiled.handles.get(hi);
			Coordinate tgt = handleTargets.get(hi);

			addHandleRhsStep1(rhs1, h.a, h.wa, tgt.x, tgt.y, w2);
			addHandleRhsStep1(rhs1, h.b, h.wb, tgt.x, tgt.y, w2);
			addHandleRhsStep1(rhs1, h.c, h.wc, tgt.x, tgt.y, w2);
		}

		DMatrixRMaj sol1 = new DMatrixRMaj(2 * V, 1);
		compiled.solverStep1.solve(rhs1, sol1);

		double[] px = new double[V];
		double[] py = new double[V];
		for (int i = 0; i < V; i++) {
			px[i] = sol1.get(2 * i, 0);
			py[i] = sol1.get(2 * i + 1, 0);
		}

		// ---------- Step 2: scale adjustment solve (x and y together) ----------
		// rhs2 = [rhsX rhsY] (V x 2)
		DMatrixRMaj rhs2 = new DMatrixRMaj(V, 2);

		// edge contributions
		for (Edge e : mesh.edges) {
			double[] rot = computeNormalizedRotation(e, px, py, params.eps);
			double ck = rot[0], sk = rot[1];

			double bx = ck * e.ex + sk * e.ey;
			double by = -sk * e.ex + ck * e.ey;

			double we = e.weight;
			int i = e.vi, j = e.vj;

			rhs2.add(i, 0, -we * bx);
			rhs2.add(j, 0, +we * bx);

			rhs2.add(i, 1, -we * by);
			rhs2.add(j, 1, +we * by);
		}

		// handle contributions
		for (int hi = 0; hi < compiled.handles.size(); hi++) {
			HandleInfo h = compiled.handles.get(hi);
			Coordinate tgt = handleTargets.get(hi);

			rhs2.add(h.a, 0, w2 * h.wa * tgt.x);
			rhs2.add(h.b, 0, w2 * h.wb * tgt.x);
			rhs2.add(h.c, 0, w2 * h.wc * tgt.x);

			rhs2.add(h.a, 1, w2 * h.wa * tgt.y);
			rhs2.add(h.b, 1, w2 * h.wb * tgt.y);
			rhs2.add(h.c, 1, w2 * h.wc * tgt.y);
		}

		DMatrixRMaj sol2 = new DMatrixRMaj(V, 2);
		compiled.solverStep2.solve(rhs2, sol2);

		double[] fx = new double[V];
		double[] fy = new double[V];
		for (int i = 0; i < V; i++) {
			fx[i] = sol2.get(i, 0);
			fy[i] = sol2.get(i, 1);
		}

		// Build output polygon by moving only the original ring vertices
		LinearRing shell = mesh.gf.createLinearRing(buildRing(mesh.shellRing, fx, fy));
		LinearRing[] holes = new LinearRing[mesh.holeRings.length];
		for (int h = 0; h < holes.length; h++) {
			holes[h] = mesh.gf.createLinearRing(buildRing(mesh.holeRings[h], fx, fy));
		}
		return mesh.gf.createPolygon(shell, holes);
	}

	/**
	 * Adds a single vertex contribution to the Step 1 normal-equation RHS for one
	 * handle.
	 *
	 * <p>
	 * Step 1 solves for (x,y) stacked into a single vector of length 2V. Handle
	 * constraints are applied via barycentric weights scaled by {@code w^2}.
	 */
	private static void addHandleRhsStep1(DMatrixRMaj rhs1, int vertexId, double bary, double targetX, double targetY, double w2) {
		int ix = 2 * vertexId;
		int iy = 2 * vertexId + 1;
		rhs1.add(ix, 0, w2 * bary * targetX);
		rhs1.add(iy, 0, w2 * bary * targetY);
	}

	/**
	 * Adds the Gram contribution of a single barycentric constraint in 1D (x or y)
	 * to a symmetric triplet matrix.
	 */
	private static void addBaryConstraintGram1D(DMatrixSparseTriplet T, int ia, int ib, int ic, double wa, double wb, double wc, double w2) {
		// outer product of [wa wb wc] scaled by w^2
		double aa = w2 * wa * wa;
		double bb = w2 * wb * wb;
		double cc = w2 * wc * wc;
		double ab = w2 * wa * wb;
		double ac = w2 * wa * wc;
		double bc = w2 * wb * wc;

		addSym(T, ia, ia, aa);
		addSym(T, ib, ib, bb);
		addSym(T, ic, ic, cc);

		addSym(T, ia, ib, ab);
		addSym(T, ia, ic, ac);
		addSym(T, ib, ic, bc);
	}

	/**
	 * Adds the Gram contribution of a barycentric constraint to the Step 2 system
	 * (same matrix used for x and y).
	 */
	private static void addBaryConstraintGram2D(DMatrixSparseTriplet T, int a, int b, int c, double wa, double wb, double wc, double w2) {
		double aa = w2 * wa * wa;
		double bb = w2 * wb * wb;
		double cc = w2 * wc * wc;
		double ab = w2 * wa * wb;
		double ac = w2 * wa * wc;
		double bc = w2 * wb * wc;

		addSym(T, a, a, aa);
		addSym(T, b, b, bb);
		addSym(T, c, c, cc);

		addSym(T, a, b, ab);
		addSym(T, a, c, ac);
		addSym(T, b, c, bc);
	}

	private static void addSym(DMatrixSparseTriplet T, int r, int c, double v) {
		T.addItem(r, c, v);
		if (r != c) {
			T.addItem(c, r, v);
		}
	}

	/**
	 * Computes the normalized 2D rotation parameters (c,s) for an edge neighborhood
	 * from the Step 1 solution.
	 *
	 * <p>
	 * Returns {@code [c, s]} representing:
	 * 
	 * <pre>
	 * [ c  s]
	 * [-s  c]
	 * </pre>
	 * 
	 * If the magnitude is near-zero, returns the identity rotation.
	 */
	private static double[] computeNormalizedRotation(Edge e, double[] px, double[] py, double eps) {
		int m = e.boundary ? 3 : 4;
		int cols = 2 * m;

		double[] vec = new double[cols];
		for (int t = 0; t < m; t++) {
			int vid = e.neigh[t];
			vec[2 * t] = px[vid];
			vec[2 * t + 1] = py[vid];
		}

		// pk is 2 x cols flattened: row0 then row1
		double ck = 0, sk = 0;
		for (int i = 0; i < cols; i++) {
			ck += e.pk[i] * vec[i];
			sk += e.pk[cols + i] * vec[i];
		}

		double s = Math.sqrt(ck * ck + sk * sk);
		if (s < eps) {
			return new double[] { 1.0, 0.0 };
		}
		return new double[] { ck / s, sk / s };
	}

	/**
	 * Locates a rest-space handle point within the triangulation and returns
	 * barycentric data for constraints.
	 *
	 * <p>
	 * Uses a spatial index over triangle envelopes as a fast candidate filter, then
	 * verifies via barycentric coordinates.
	 *
	 * @throws IllegalStateException if no containing triangle can be found (should
	 *                               be rare if inputs are valid)
	 */
	private HandleInfo locateHandle(Coordinate p) {
		Envelope env = new Envelope(p);
		@SuppressWarnings("unchecked")
		List<Integer> candidates = mesh.triIndex.query(env);

		for (int triId : candidates) {
			HandleInfo hi = tryBarycentric(triId, p);
			if (hi != null) {
				return hi;
			}
		}

		// fallback (rare)
		for (int triId = 0; triId < mesh.triangles.length; triId++) {
			HandleInfo hi = tryBarycentric(triId, p);
			if (hi != null) {
				return hi;
			}
		}
		throw new IllegalStateException("Failed to locate containing triangle for handle point: " + p);
	}

	/**
	 * Computes barycentric coordinates of {@code p} with respect to triangle
	 * {@code triId}.
	 *
	 * @return handle info if inside (with small slack), otherwise {@code null}
	 */
	private HandleInfo tryBarycentric(int triId, Coordinate p) {
		int[] t = mesh.triangles[triId];
		int a = t[0], b = t[1], c = t[2];

		double ax = mesh.rx[a], ay = mesh.ry[a];
		double bx = mesh.rx[b], by = mesh.ry[b];
		double cx = mesh.rx[c], cy = mesh.ry[c];

		double v0x = bx - ax, v0y = by - ay;
		double v1x = cx - ax, v1y = cy - ay;
		double v2x = p.x - ax, v2y = p.y - ay;

		double den = cross(v0x, v0y, v1x, v1y);
		if (Math.abs(den) < params.eps) {
			return null;
		}

		// barycentric with oriented areas
		double wb = cross(v2x, v2y, v1x, v1y) / den;
		double wc = cross(v0x, v0y, v2x, v2y) / den;
		double wa = 1.0 - wb - wc;

		double e = 1e-10; // containment slack
		if (wa >= -e && wb >= -e && wc >= -e && wa <= 1 + e && wb <= 1 + e && wc <= 1 + e) {
			return new HandleInfo(a, b, c, wa, wb, wc);
		}
		return null;
	}

	private static double cross(double ax, double ay, double bx, double by) {
		return ax * by - ay * bx;
	}

	/**
	 * Builds a closed JTS coordinate ring from mesh vertex indices and solved
	 * vertex positions.
	 */
	private static Coordinate[] buildRing(int[] ringIdx, double[] x, double[] y) {
		Coordinate[] out = new Coordinate[ringIdx.length + 1];
		for (int i = 0; i < ringIdx.length; i++) {
			int vid = ringIdx[i];
			out[i] = new Coordinate(x[vid], y[vid]);
		}
		out[ringIdx.length] = new Coordinate(out[0]);
		return out;
	}

	/**
	 * Builds the internal triangle mesh from an externally-provided triangulation
	 * geometry instead of computing one from scratch. This mirrors the logic of
	 * buildMesh(Polygon, Params) but skips refinement and uses the provided
	 * triangulation Geometry.
	 */
	private static Mesh buildMesh(Polygon poly, Geometry triGeom, Params p) {
		GeometryFactory gf = poly.getFactory();
		IndexedPointInAreaLocator locator = new IndexedPointInAreaLocator(poly);

		// collect ring coords (excluding closure)
		List<Coordinate> boundary = new ArrayList<>();
		Coordinate[] shell = stripClosure(poly.getExteriorRing().getCoordinates());
		boundary.addAll(Arrays.asList(shell));

		Coordinate[][] holes = new Coordinate[poly.getNumInteriorRing()][];
		for (int i = 0; i < poly.getNumInteriorRing(); i++) {
			holes[i] = stripClosure(poly.getInteriorRingN(i).getCoordinates());
			boundary.addAll(Arrays.asList(holes[i]));
		}

		// Sites start as boundary vertices (unique)
		LinkedHashMap<CoordKey, Coordinate> sites = new LinkedHashMap<>(boundary.size() * 2);
		for (Coordinate c : boundary) {
			sites.put(new CoordKey(c), new Coordinate(c));
		}

		// Use the provided triangulation geometry to build the Triangulation
		Triangulation tri = triangulateFromGeometry(triGeom, locator, sites.values());

		// map original ring vertices to indices in final mesh
		int[] shellRing = new int[shell.length];
		for (int i = 0; i < shell.length; i++) {
			shellRing[i] = tri.findVertex(shell[i], 1e-8);
		}

		int[][] holeRings = new int[holes.length][];
		for (int h = 0; h < holes.length; h++) {
			holeRings[h] = new int[holes[h].length];
			for (int i = 0; i < holes[h].length; i++) {
				holeRings[h][i] = tri.findVertex(holes[h][i], 1e-8);
			}
		}

		// build triangle spatial index
		HPRtree index = new HPRtree();
		for (int i = 0; i < tri.triangles.size(); i++) {
			int[] t = tri.triangles.get(i);
			Envelope env = new Envelope(new Coordinate(tri.rx[t[0]], tri.ry[t[0]]));
			env.expandToInclude(tri.rx[t[1]], tri.ry[t[1]]);
			env.expandToInclude(tri.rx[t[2]], tri.ry[t[2]]);
			index.insert(env, i);
		}

		// build edges (with neighbor context) + precompute hk, pk, weights
		List<Edge> edges = buildEdgesAndPrecompute(tri.rx, tri.ry, tri.triangles, p);

		// triangles array for handle lookup
		int[][] trianglesArr = tri.triangles.toArray(new int[0][]);

		return new Mesh(gf, tri.rx, tri.ry, trianglesArr, edges, shellRing, holeRings, index, locator);
	}

	/**
	 * Removes the duplicated closing coordinate from a linear ring if present.
	 */
	private static Coordinate[] stripClosure(Coordinate[] ring) {
		if (ring.length >= 2 && ring[0].equals2D(ring[ring.length - 1])) {
			return Arrays.copyOf(ring, ring.length - 1);
		}
		return ring;
	}

	/**
	 * User-tunable parameters affecting stability and constraint weighting.
	 *
	 * <p>
	 * Most callers only need to adjust {@link #targetVertexCount},
	 * {@link #handleWeight}, and {@link #useCotangentWeights}.
	 */
	public static final class Params {
		/** Constraint weight (paper uses 1000). */
		public double handleWeight = 1000.0;
		/**
		 * Use cotangent weights for edges (recommended for irregular triangulations).
		 */
		public boolean useCotangentWeights = true;
		/** Robustness eps for barycentric and degeneracy. */
		public double eps = 1e-12;
		/** If true, clamp cotangent weights to be >= eps to preserve SPD. */
		public boolean clampWeightsPositive = true;
	}

	/**
	 * Immutable, reusable compiled state for a specific handle set.
	 *
	 * <p>
	 * Instances are created by {@link #compile(List)} and can be reused across many
	 * calls to {@link #deform(Compiled, List)} as long as the handle count/order is
	 * unchanged.
	 */
	public static final class Compiled {
		private final List<HandleInfo> handles;
		private final LinearSolverSparse<DMatrixSparseCSC, DMatrixRMaj> solverStep1; // (2V)x(2V)
		private final LinearSolverSparse<DMatrixSparseCSC, DMatrixRMaj> solverStep2; // (V)x(V)
		public final int vertexCount;

		private Compiled(List<HandleInfo> handles, LinearSolverSparse<DMatrixSparseCSC, DMatrixRMaj> solverStep1,
				LinearSolverSparse<DMatrixSparseCSC, DMatrixRMaj> solverStep2, int vertexCount) {
			this.handles = handles;
			this.solverStep1 = solverStep1;
			this.solverStep2 = solverStep2;
			this.vertexCount = vertexCount;
		}
	}

	private static final class Mesh {
		final GeometryFactory gf;

		final double[] rx, ry; // rest vertex positions
		final int[][] triangles; // vertex indices
		final List<Edge> edges; // undirected edges with neighbor context

		// ring vertex indices in mesh, in original ring order (excluding closing point)
		final int[] shellRing;
		final int[][] holeRings;

		// spatial index for triangles (for handle barycentrics)
		final SpatialIndex triIndex;

		final IndexedPointInAreaLocator locator;

		Mesh(GeometryFactory gf, double[] rx, double[] ry, int[][] triangles, List<Edge> edges, int[] shellRing, int[][] holeRings, SpatialIndex triIndex,
				IndexedPointInAreaLocator locator) {
			this.gf = gf;
			this.rx = rx;
			this.ry = ry;
			this.triangles = triangles;
			this.edges = edges;
			this.shellRing = shellRing;
			this.holeRings = holeRings;
			this.triIndex = triIndex;
			this.locator = locator;
		}

		int vertexCount() {
			return rx.length;
		}
	}

	private static final class Edge {
		final int vi, vj; // direction chosen consistently (vi < vj)
		final int[] neigh; // vi,vj,vl,(vr)
		final boolean boundary;
		final double ex, ey; // rest edge vector vj-vi
		final double[] pk; // 2 x (2*m) = [row0... row1...] flattened
		final double[][] hk; // 2 x (2*m) (constant for step1)
		final double weight; // edge weight (1 or cotangent)

		Edge(int vi, int vj, int[] neigh, boolean boundary, double ex, double ey, double[] pk, double[][] hk, double weight) {
			this.vi = vi;
			this.vj = vj;
			this.neigh = neigh;
			this.boundary = boundary;
			this.ex = ex;
			this.ey = ey;
			this.pk = pk;
			this.hk = hk;
			this.weight = weight;
		}
	}

	private static final class HandleInfo {
		final int a, b, c;
		final double wa, wb, wc;

		HandleInfo(int a, int b, int c, double wa, double wb, double wc) {
			this.a = a;
			this.b = b;
			this.c = c;
			this.wa = wa;
			this.wb = wb;
			this.wc = wc;
		}
	}

	private static DMatrixSparseCSC toCsc(DMatrixSparseTriplet triplet) {
		DMatrixSparseCSC csc = new DMatrixSparseCSC(triplet.numRows, triplet.numCols, triplet.nz_length);
		DConvertMatrixStruct.convert(triplet, csc);
//		csc.sortIndices(null); // faster without
		// IMPORTANT: merge duplicates by summing them
		CommonOps_DSCC.duplicatesAdd(csc, null);
		return csc;
	}

	private static final class Triangulation {
		final double[] rx, ry;
		final List<int[]> triangles;
		final Map<CoordKey, Integer> indexMap;

		Triangulation(double[] rx, double[] ry, List<int[]> triangles, LinkedHashMap<CoordKey, Coordinate> sitesExact, Map<CoordKey, Integer> indexMap) {
			this.rx = rx;
			this.ry = ry;
			this.triangles = triangles;
			this.indexMap = indexMap;
		}

		int findVertex(Coordinate c, double tol) {
			Integer idx = indexMap.get(new CoordKey(c));
			if (idx != null) {
				return idx;
			}

			// fallback nearest
			int best = -1;
			double bestD2 = Double.POSITIVE_INFINITY;
			for (int i = 0; i < rx.length; i++) {
				double dx = rx[i] - c.x, dy = ry[i] - c.y;
				double d2 = dx * dx + dy * dy;
				if (d2 < bestD2) {
					bestD2 = d2;
					best = i;
				}
			}
			if (best >= 0 && bestD2 <= tol * tol) {
				return best;
			}
			throw new IllegalStateException("Boundary vertex not found in triangulation: " + c);
		}
	}

	/**
	 * Parses a triangle Geometry into a Triangulation instance.
	 *
	 * The method iterates all component geometries (expected Polygons) in
	 * {@code trisGeom}, computes centroids and filters out triangles whose
	 * centroids are outside the polygon (using the provided locator), and
	 * constructs the vertex arrays, triangle index list and index map.
	 *
	 * @param trisGeom triangle geometry (collection of triangle polygons)
	 * @param locator  locator for testing triangle centroid containment
	 * @param sites    original site coordinates (used to populate sitesExact)
	 */
	private static Triangulation triangulateFromGeometry(Geometry trisGeom, IndexedPointInAreaLocator locator, Collection<Coordinate> sites) {
		Map<CoordKey, Integer> index = new HashMap<>();
		int N = trisGeom.getNumGeometries();
		List<Double> xs = new ArrayList<>(N);
		List<Double> ys = new ArrayList<>(N);
		List<int[]> triangles = new ArrayList<>(N);

		Coordinate centroid = new Coordinate(); // reuse object

		for (int gi = 0; gi < N; gi++) {
			Geometry g = trisGeom.getGeometryN(gi);
			if (!(g instanceof Polygon)) {
				continue;
			}
			Polygon tPoly = (Polygon) g;
			Coordinate[] c = tPoly.getExteriorRing().getCoordinates();
			if (c.length < 4) {
				continue;
			}

			Coordinate aC = c[0], bC = c[1], cC = c[2];

			centroid.x = (aC.x + bC.x + cC.x) / 3.0;
			centroid.y = (aC.y + bC.y + cC.y) / 3.0;

			if (locator.locate(centroid) == Location.EXTERIOR) {
				// just in case
				continue;
			}

			int a = getOrAdd(index, xs, ys, aC);
			int bI = getOrAdd(index, xs, ys, bC);
			int cI = getOrAdd(index, xs, ys, cC);

			if (a == bI || bI == cI || cI == a) {
				continue;
			}
			triangles.add(new int[] { a, bI, cI });
		}

		double[] rx = xs.stream().mapToDouble(d -> d).toArray();
		double[] ry = ys.stream().mapToDouble(d -> d).toArray();

		LinkedHashMap<CoordKey, Coordinate> sitesExact = new LinkedHashMap<>();
		for (Coordinate s : sites) {
			sitesExact.put(new CoordKey(s), s);
		}

		return new Triangulation(rx, ry, triangles, sitesExact, index);
	}

	private static int getOrAdd(Map<CoordKey, Integer> index, List<Double> xs, List<Double> ys, Coordinate c) {
		CoordKey key = new CoordKey(c);
		Integer idx = index.get(key);
		if (idx != null) {
			return idx;
		}
		int id = xs.size();
		xs.add(c.x);
		ys.add(c.y);
		index.put(key, id);
		return id;
	}

	private static final class CoordKey {
		// quantize to stabilize hashing across JTS floating noise
		private static final double Q = 1e-12;
		final long xq, yq;

		CoordKey(Coordinate c) {
			this.xq = Math.round(c.x / Q);
			this.yq = Math.round(c.y / Q);
		}

		@Override
		public boolean equals(Object o) {
			if (this == o) {
				return true;
			}
			if (!(o instanceof CoordKey k)) {
				return false;
			}
			return xq == k.xq && yq == k.yq;
		}

		@Override
		public int hashCode() {
			return Objects.hash(xq, yq);
		}
	}

	private static final class Adj {
		final int opp; // opposite vertex index in this triangle

		Adj(int opp) {
			this.opp = opp;
		}
	}

	private static List<Edge> buildEdgesAndPrecompute(double[] rx, double[] ry, List<int[]> tris, Params p) {
		// map undirected edge -> list of opposite vertices (one per adjacent triangle)
		Map<Long, List<Adj>> adj = new HashMap<>();

		for (int[] t : tris) {
			int a = t[0], b = t[1], c = t[2];
			addAdj(adj, a, b, c);
			addAdj(adj, b, c, a);
			addAdj(adj, c, a, b);
		}

		List<Edge> edges = new ArrayList<>(adj.size());

		for (Map.Entry<Long, List<Adj>> en : adj.entrySet()) {
			long key = en.getKey();
			int vi = (int) (key >>> 32);
			int vj = (int) (key & 0xffffffffL);

			List<Adj> list = en.getValue();
			boolean boundary = list.size() == 1;
			if (!(list.size() == 1 || list.size() == 2)) {
				continue;
			}

			int vl = list.get(0).opp;
			int[] neigh;
			if (boundary) {
				neigh = new int[] { vi, vj, vl };
			} else {
				int vr = list.get(1).opp;
				neigh = new int[] { vi, vj, vl, vr };
			}

			double ex = rx[vj] - rx[vi];
			double ey = ry[vj] - ry[vi];

			// compute Pk and hk
			Precomp ph = precomputePkHk(rx, ry, neigh, ex, ey);

			double weight = 1.0;
			if (p.useCotangentWeights) {
				double cot1 = cotangentOpp(rx, ry, vi, vj, vl, p.eps);
				double cot2 = boundary ? 0.0 : cotangentOpp(rx, ry, vi, vj, neigh[3], p.eps);
				weight = 0.5 * (cot1 + cot2);
				if (p.clampWeightsPositive) {
					weight = Math.max(weight, p.eps);
				}
			}

			edges.add(new Edge(vi, vj, neigh, boundary, ex, ey, ph.pk, ph.hk, weight));
		}

		return edges;
	}

	private static void addAdj(Map<Long, List<Adj>> adj, int u, int v, int opp) {
		int a = Math.min(u, v);
		int b = Math.max(u, v);
		long key = (((long) a) << 32) | (b & 0xffffffffL);
		adj.computeIfAbsent(key, k -> new ArrayList<>(2)).add(new Adj(opp));
	}

	private static double cotangentOpp(double[] rx, double[] ry, int vi, int vj, int vk, double eps) {
		// angle at vk opposite edge (vi,vj)
		double ux = rx[vi] - rx[vk];
		double uy = ry[vi] - ry[vk];
		double vx = rx[vj] - rx[vk];
		double vy = ry[vj] - ry[vk];

		double dot = ux * vx + uy * vy;
		double cr = Math.abs(cross(ux, uy, vx, vy));
		if (cr < eps) {
			return 0.0;
		}
		return dot / cr;
	}

	private static final class Precomp {
		final double[] pk; // 2 x cols flattened
		final double[][] hk; // 2 x cols

		Precomp(double[] pk, double[][] hk) {
			this.pk = pk;
			this.hk = hk;
		}
	}

	private static Precomp precomputePkHk(double[] rx, double[] ry, int[] neigh, double ex, double ey) {
		int m = neigh.length;
		int cols = 2 * m;

		// Build Gk: (2m x 2)
		// rows per vertex v: [vx vy; vy -vx]
		double[][] G = new double[2 * m][2];
		for (int t = 0; t < m; t++) {
			int id = neigh[t];
			double vx = rx[id], vy = ry[id];
			G[2 * t][0] = vx;
			G[2 * t][1] = vy;
			G[2 * t + 1][0] = vy;
			G[2 * t + 1][1] = -vx;
		}

		// GTG = G^T G (2x2)
		double a00 = 0, a01 = 0, a11 = 0;
		for (int r = 0; r < 2 * m; r++) {
			double g0 = G[r][0];
			double g1 = G[r][1];
			a00 += g0 * g0;
			a01 += g0 * g1;
			a11 += g1 * g1;
		}
		double det = a00 * a11 - a01 * a01;
		if (Math.abs(det) < 1e-18) {
			// degenerate neighborhood; fall back to identity-ish
			double[] pk = new double[2 * cols];
			pk[0] = 1; // crude
			double[][] hk = constantM(m);
			return new Precomp(pk, hk);
		}
		double inv00 = a11 / det;
		double inv01 = -a01 / det;
		double inv11 = a00 / det;

		// Compute Pk = (GTG)^-1 G^T : 2 x (2m) (cols = 2m)
		double[] pk = new double[2 * cols];
		for (int c = 0; c < cols; c++) {
			// column of G^T is row c of G
			double gt0 = G[c][0];
			double gt1 = G[c][1];
			pk[c] = inv00 * gt0 + inv01 * gt1; // row0
			pk[cols + c] = inv01 * gt0 + inv11 * gt1; // row1
		}

		// E = [[ex ey],[ey -ex]]
		double e00 = ex, e01 = ey;
		double e10 = ey, e11 = -ex;

		// hk = M - E*Pk
		double[][] M = constantM(m);
		double[][] hk = new double[2][cols];
		for (int c = 0; c < cols; c++) {
			double p0 = pk[c];
			double p1 = pk[cols + c];
			double ep0 = e00 * p0 + e01 * p1;
			double ep1 = e10 * p0 + e11 * p1;
			hk[0][c] = M[0][c] - ep0;
			hk[1][c] = M[1][c] - ep1;
		}

		return new Precomp(pk, hk);
	}

	private static double[][] constantM(int m) {
		// m=3 -> cols=6: [-1 0 1 0 0 0; 0 -1 0 1 0 0]
		// m=4 -> cols=8: [-1 0 1 0 0 0 0 0; 0 -1 0 1 0 0 0 0]
		int cols = 2 * m;
		double[][] M = new double[2][cols];
		M[0][0] = -1;
		M[0][2] = 1;
		M[1][1] = -1;
		M[1][3] = 1;
		return M;
	}
}