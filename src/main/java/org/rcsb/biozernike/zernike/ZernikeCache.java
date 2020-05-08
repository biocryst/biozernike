package org.rcsb.biozernike.zernike;

import org.rcsb.biozernike.complex.Complex;

import java.util.ArrayList;
import java.util.List;

/**
 * Precomputed data-independent coefficients needed to calculate
 * 3D Zernike moments up to the order of 20.
 * Code adapted from the ZernikeMoments C++ library (Novotni-Klein),
 * see http://www.cg.cs.uni-bonn.de/project-pages/3dsearch/
 */
public class ZernikeCache {

	// TODO move this to a parameter and create init method to be called by user when calculation is needed
	private static final int maxN = 21;
	private static final double[] clmCache = new double[maxN * maxN];
	private static final double[] qnlmCache = new double[maxN * maxN * maxN];
	private static final List<List<List<List<ComplexCoeff>>>> gCoefCache = new ArrayList<>(maxN);

	static {
		// c[l,m] init
		for (int l = 0; l < maxN; l++) {
			for (int m = 0; m <= l; m++) {
				double nSqr = (2.0 * l + 1.0) * FactorialCache.getValue(l + 1, l + m);
				double dSqr = FactorialCache.getValue(l - m + 1, l);
				clmCache[l * maxN + m] = Math.sqrt(nSqr / dSqr);
			}
		}

		// q[n,l,mu] init
		for (int n = 0; n < maxN; ++n) {
			int l0 = n % 2;
			for (int l = l0; l <= n; l += 2) {
				int k = (n - l) / 2;
				for (int mu = 0; mu <= k; ++mu) {
					double nom = BinomialCache.getValue(2 * k, k) * BinomialCache.getValue(k, mu) *
							BinomialCache.getValue(2 * (k + l + mu) + 1, 2 * k);
					if ((k + mu) % 2 != 0) {
						nom = -nom;
					}
					double den = Math.pow(2.0, 2.0 * k) * BinomialCache.getValue(k + l + mu, k);
					double n_sqrt = 2.0 * l + 4.0 * k + 3.0;
					qnlmCache[n * maxN * maxN + (l / 2) * maxN + mu] = nom / den * Math.sqrt(n_sqrt / 3.0);
				}
			}
		}

		// gCoeffs init
		for (int n = 0; n < maxN; ++n) {
			List<List<List<ComplexCoeff>>> gCoefNLevel = new ArrayList<>(n / 2 + 1);
			int li = 0, l0 = n % 2;

			for (int l = l0; l <= n; ++li, l += 2) {
				List<List<ComplexCoeff>> gCoefLLevel = new ArrayList<>(l + 1);

				for (int m = 0; m <= l; ++m) {
					List<ComplexCoeff> gCoefMLevel = new ArrayList<>(l + 1);
					double w = 1.0 / Math.pow(2.0, m);
					int k = (n - l) / 2;

					for (int nu = 0; nu <= k; ++nu) {
						double w_Nu = w * getQnlmValue(n, li, nu);

						for (int alpha = 0; alpha <= nu; ++alpha) {
							double w_NuA = w_Nu * BinomialCache.getValue(nu, alpha);

							for (int beta = 0; beta <= nu - alpha; ++beta) {
								double w_NuAB = w_NuA * BinomialCache.getValue(nu - alpha, beta);

								for (int p = 0; p <= m; ++p) {
									double w_NuABP = w_NuAB * BinomialCache.getValue(m, p);

									for (int mu = 0; mu <= (l - m) / 2; ++mu) {
										double w_NuABPMu = w_NuABP *
												BinomialCache.getValue(l, mu) *
												BinomialCache.getValue(l - mu, m + mu) /
												Math.pow(2.0, 2.0 * mu);

										for (int q = 0; q <= mu; ++q) {
											double w_NuABPMuQ = w_NuABPMu * BinomialCache.getValue(mu, q);
											if ((m - p + mu) % 2 != 0) {
												w_NuABPMuQ = -w_NuABPMuQ;
											}

											int rest = p % 4;
											Complex c = null;

											switch (rest) {
												case 0:
													c = new Complex(w_NuABPMuQ, 0);
													break;
												case 1:
													c = new Complex(0, w_NuABPMuQ);
													break;
												case 2:
													c = new Complex(-w_NuABPMuQ, 0);
													break;
												case 3:
													c = new Complex(0, -w_NuABPMuQ);
													break;
											}

											int z_i = l - m + 2 * (nu - alpha - beta - mu);
											int y_i = 2 * (mu - q + beta) + m - p;
											int x_i = 2 * q + p + 2 * alpha;
											ComplexCoeff cc = new ComplexCoeff(x_i, y_i, z_i, c.conj());
											gCoefMLevel.add(cc);
										}
									}
								}
							}
						}
					}
					gCoefLLevel.add(gCoefMLevel);
				}
				gCoefNLevel.add(gCoefLLevel);
			}
			gCoefCache.add(gCoefNLevel);
		}
	}

	public static double getClmValue(int l, int m) {
		return clmCache[l * maxN + m];
	}

	static double getQnlmValue(int n, int l, int m) {
		return qnlmCache[n * maxN * maxN + l * maxN + m];
	}

	public static List<ComplexCoeff> getGCoefs(int n, int l, int m) {
		return gCoefCache.get(n).get(l).get(m);
	}
}


