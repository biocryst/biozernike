package org.rcsb.biozernike;

import org.rcsb.biozernike.complex.Complex;
import org.rcsb.biozernike.complex.PolynomialSolver;
import org.rcsb.biozernike.volume.Volume;
import org.rcsb.biozernike.zernike.BinomialCache;
import org.rcsb.biozernike.zernike.ZernikeCache;
import org.rcsb.biozernike.zernike.ZernikeMoments;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;
import java.util.*;

public class InvariantNorm {
	private static final Logger logger = LoggerFactory.getLogger(InvariantNorm.class);

	private ZernikeMoments moments = null;
	private Vector3d center = null;
	private int maxOrder;

	private Map<Map.Entry<Integer, Integer>, List<MomentTransform>> transformsMap;

	// list of equivalent solutions to the normalizing equations. 4 or 8 (depending on the orders used).
	// transformsMap are scaled
	private PolynomialSolver solver = new PolynomialSolver();

	// invariantsMap - mean moments amplitude in the equivalent solutions.
	private Map<Integer, List<Double>> invariantsMap;

	private InvariantNorm() {
	}

	public InvariantNorm(ZernikeMoments moments) {
		this.moments = moments;
		this.center = new Vector3d(moments.getVolume().getCenterReal());
		this.maxOrder = moments.getMaxOrder();
		this.transformsMap = new HashMap<>();
		this.invariantsMap = new HashMap<>();
	}


	public InvariantNorm(Volume volume, int maxOrder) {
		this.moments = new ZernikeMoments(volume, maxOrder);
		this.center = new Vector3d(volume.getCenterReal());
		this.maxOrder = maxOrder;
		this.transformsMap = new HashMap<>();
		this.invariantsMap = new HashMap<>();
	}


	public InvariantNorm(Map<Map.Entry<Integer, Integer>, List<MomentTransform>> transformsMap, double[] center) {
		this.center = new Vector3d(center);
		this.transformsMap = transformsMap;
		this.invariantsMap = new HashMap<>();
	}

	public static List<Matrix4d> alignMultiple(List<InvariantNorm> invariantNorms) {
		List<Matrix3d> rotations = RotationAlignment.alignMultiple(invariantNorms);
		return addTranslation(invariantNorms, rotations);
	}

	public static List<Matrix4d> alignMultiple(List<InvariantNorm> invariantNorms, int maxNormInd, boolean useExistingTransforms) {
		List<Matrix3d> rotations = RotationAlignment.alignMultiple(invariantNorms, maxNormInd, useExistingTransforms);
		return addTranslation(invariantNorms, rotations);
	}

	private static List<Matrix4d> addTranslation(List<InvariantNorm> invariantNorms, List<Matrix3d> rotations) {
		List<Matrix4d> alignments = new ArrayList<>();

		for (int i = 0; i < invariantNorms.size(); i++) {
			Matrix4d m = new Matrix4d(rotations.get(i), invariantNorms.get(i).getCenter(), 1);
			m.invert();
			alignments.add(m);
		}

		return alignments;

	}

	public int getMaxOrder() {
		return maxOrder;
	}

	public void setMaxOrder(int maxOrderNew) {

		int nMomentsNew = moments.setMaxOrder(maxOrderNew);
		if (maxOrderNew < maxOrder) {
			trimSolutions(nMomentsNew);
		}
		if (maxOrderNew > maxOrder) {
			invariantsMap = new HashMap<>();
			transformsMap = new HashMap<>();
		}
		maxOrder = maxOrderNew;
	}

	private void trimSolutions(int nRemain) {
		for (List<Double> invariants : invariantsMap.values()) {
			invariants.subList(nRemain, invariants.size()).clear();
		}
		for (List<MomentTransform> transforms : transformsMap.values()) {
			for (MomentTransform transform : transforms) {
				transform.trimMoments(nRemain);
			}
		}
	}

	// get moments which correspond to a rotation parametrized by complex numbers a and b
	// see Cayley-Klein parametrization
	private List<List<List<Complex>>> Rotate(Complex a, Complex b) {
		int maxOrder = moments.getMaxOrder();
		List<List<List<Complex>>> zmRotated = new ArrayList<>(maxOrder + 1);

		double aac = a.mul(a.conj()).getReal();
		double bbc = b.mul(b.conj()).getReal();
		double bbcaac = -bbc / aac;
		Complex abc = a.divide(b.conj()).negate();
		Complex ab = a.divide(b);

		double[] bbcaac_pow_k_list = new double[maxOrder + 1];
		double[] aac_pow_l_list = new double[maxOrder + 1];
		Complex[] ab_pow_m_list = new Complex[maxOrder + 1];
		Complex[] abc_pow_mu_list = new Complex[2 * maxOrder + 1];

		bbcaac_pow_k_list[0] = 1.0;
		aac_pow_l_list[0] = 1.0;
		ab_pow_m_list[0] = new Complex(1, 0);
		abc_pow_mu_list[maxOrder] = new Complex(1, 0);

		for (int i = 1; i <= maxOrder; i++) {
			bbcaac_pow_k_list[i] = bbcaac_pow_k_list[i - 1] * bbcaac;
			aac_pow_l_list[i] = aac_pow_l_list[i - 1] * aac;
			ab_pow_m_list[i] = ab_pow_m_list[i - 1].mul(ab);

			abc_pow_mu_list[maxOrder + i] = abc_pow_mu_list[maxOrder + i - 1].mul(abc);
			abc_pow_mu_list[maxOrder - i] = abc_pow_mu_list[maxOrder - i + 1].divide(abc);
		}

		for (int n = 0; n <= maxOrder; ++n) {
			List<List<Complex>> zmLLevel = new ArrayList<>(n / 2 + 1);

			int l0 = n % 2, li = 0;
			for (int l = l0; l <= n; ++li, l += 2) {
				List<Complex> zmMLevel = new ArrayList<>(l + 1);

				for (int m = 0; m <= l; ++m) {

					Complex mu_sum = new Complex(0, 0);
					Complex abm = ab_pow_m_list[m];

					for (int mu = -l; mu <= l; mu++) {

						Complex F_exp;
						if (mu >= 0) {
							F_exp = moments.getOriginalUnscaled(n, l / 2, mu);
						} else {
							F_exp = moments.getOriginalUnscaled(n, l / 2, -mu).conj();
							if (mu % 2 != 0) {
								F_exp = F_exp.negate();
							}
						}

						double k_sum = 0;

						int maxM = Math.max(mu, m);
						int minM = Math.min(0, m + mu);

						for (int k = maxM; k <= l + minM; k++) {
							double bin = BinomialCache.getValue(l - mu, k - mu) * BinomialCache.getValue(l + mu, k - m);
							k_sum += bbcaac_pow_k_list[k] * bin;
						}

						mu_sum = mu_sum.add(F_exp.mul(abc_pow_mu_list[maxOrder + mu]).mul(k_sum));

					}
					double zmValCoef = aac_pow_l_list[l] * ZernikeCache.getClmValue(l, m);
					zmMLevel.add(abm.mul(mu_sum).mul(zmValCoef));
				}
				zmLLevel.add(zmMLevel);
			}
			zmRotated.add(zmLLevel);
		}
		return zmRotated;
	}

	// Populate a list of rotated moments which correspond to
	// M_{Z2}^{2}=0 (even Z), M_{Z1}^{1}=0 (odd Z)
	// Im(M_{R2}^{1})=0 (even rotation), Im(M_{R1}^{1})=0 (odd rotation)
	private List<MomentTransform> computeRotations(int indZero, int indReal) {

		List<MomentTransform> normalizationSolutions = new ArrayList<>();

		if (moments == null) {
			logger.error("Moments are not initialized, cannot compute any new normalizations.");
			return normalizationSolutions;
		}

		if (indZero > moments.getMaxOrder() || indReal > moments.getMaxOrder()) {
			logger.error("Cannot normalize using orders ({}, {}). The moments are calculated up to the order {}.",
					indZero, indReal, moments.getMaxOrder());
			return normalizationSolutions;
		}

		Complex[] abconj_sol;
		int n_abconj;

		if (indZero % 2 == 0) {
			Complex[] abconj_coef = {
					moments.getOriginalUnscaled(indZero, 1, 2),
					moments.getOriginalUnscaled(indZero, 1, 1).negate(),
					moments.getOriginalUnscaled(indZero, 1, 0),
					moments.getOriginalUnscaled(indZero, 1, 1).conj(),
					moments.getOriginalUnscaled(indZero, 1, 2).conj()
			};
			abconj_sol = solver.roots(abconj_coef);
			n_abconj = 4;

		} else {

			Complex[] abconj_coef = {
					moments.getOriginalUnscaled(indZero, 0, 1),
					moments.getOriginalUnscaled(indZero, 0, 0).negate(),
					moments.getOriginalUnscaled(indZero, 0, 1).conj().negate()
			};
			abconj_sol = solver.roots(abconj_coef);
			n_abconj = 2;
		}


		for (int i_abconj = 0; i_abconj < n_abconj; i_abconj++) {
			Complex cur_abconj_sol = abconj_sol[i_abconj];

			double k_re = cur_abconj_sol.getReal();
			double k_im = cur_abconj_sol.getImaginary();
			double k_im2 = k_im * k_im;
			double k_re2 = k_re * k_re;
			double[] bimbre_sol_real = new double[4];

			if (indReal % 2 == 0) {
				double f20 = moments.getOriginalUnscaled(indReal, 1, 0).getReal();
				Complex f21 = moments.getOriginalUnscaled(indReal, 1, 1);
				Complex f22 = moments.getOriginalUnscaled(indReal, 1, 2);
				double f21_im = f21.getImaginary();
				double f21_re = f21.getReal();
				double f22_im = f22.getImaginary();
				double f22_re = f22.getReal();

				double k_im3 = k_im * k_im2;
				double k_im4 = k_im2 * k_im2;
				double k_re4 = k_re2 * k_re2;

				// Magic. Thanks Mathematica!
				double coef4 = 4 * f22_re * k_im * (-1 + k_im2 - 3 * k_re2) - 4 * f22_im * k_re * (1 - 3 * k_im2 + k_re2) - 2 * f21_re * k_im * k_re * (-3 + k_im2 + k_re2) +
						2 * f20 * k_im * (-1 + k_im2 + k_re2) + f21_im * (1 - 6 * k_im2 + k_im2 * k_im2 - k_re2 * k_re2);
				double coef3 = 2 * (-4 * f22_im * (k_im + k_im3 - 3 * k_im * k_re2) + f21_re * (-1 + k_im4 + 6 * k_re2 - k_re4) + 2 * k_re * (f22_re * (2 + 6 * k_im2 - 2 * k_re2) +
						f21_im * k_im * (-3 + k_im2 + k_re2) + f20 * (-1 + k_im2 + k_re2)));

				Complex[] bimbre_coef = {
						new Complex(coef4, 0),
						new Complex(coef3, 0),
						new Complex(0, 0),
						new Complex(coef3, 0),
						new Complex(-coef4, 0)
				};

				Complex[] bimbre_sol = solver.roots(bimbre_coef);
				bimbre_sol_real[0] = bimbre_sol[0].getReal();
				bimbre_sol_real[1] = bimbre_sol[1].getReal();
				bimbre_sol_real[2] = bimbre_sol[2].getReal();
				bimbre_sol_real[3] = bimbre_sol[3].getReal();

			} else {

				double f10 = moments.getOriginalUnscaled(indReal, 0, 0).getReal();
				Complex f11 = moments.getOriginalUnscaled(indReal, 0, 1);
				double f11_im = f11.getImaginary();
				double f11_re = f11.getReal();

				// More magic.
				Complex[] bimbre_coef = {
						new Complex(k_im * (f10 - 2 * f11_re * k_re) + f11_im * (-1 + k_im2 - k_re2), 0),
						new Complex(2 * (k_re * (f10 + 2 * f11_im * k_im) + f11_re * (1 + k_im2 - k_re2)), 0),
						new Complex(2 * f11_re * k_im * k_re - f10 * k_im + f11_im * (1 - k_im2 + k_re2), 0)
				};

				Complex[] bimbre_sol = solver.roots(bimbre_coef);

				bimbre_sol_real[0] = bimbre_sol[0].getReal();
				bimbre_sol_real[1] = bimbre_sol[1].getReal();
				bimbre_sol_real[2] = 0;
				bimbre_sol_real[3] = 0;

			}

			for (int i_bimbre = 0; i_bimbre < 4; i_bimbre++) {
				double cur_bimbre_sol = bimbre_sol_real[i_bimbre];

				if (Math.abs(cur_bimbre_sol) < 0.0000001) {
					continue;
				}

				double bre = 1 / Math.sqrt((1 + k_im2 + k_re2) * (1 + Math.pow(cur_bimbre_sol, 2)));
				double bim = cur_bimbre_sol * bre;
				Complex b = new Complex(bre, bim);
				Complex a = cur_abconj_sol.mul(b.conj());

				MomentTransform transform = new MomentTransform();
				transform.setRotation(a, b);
				transform.setMoments(Rotate(a, b));
				normalizationSolutions.add(transform);
			}
		}

		return normalizationSolutions;
	}

	public List<MomentTransform> getNormalizationSolutions(int indZero, int indReal) {
		Map.Entry<Integer, Integer> normKey = new AbstractMap.SimpleImmutableEntry<>(indZero, indReal);
		if (!transformsMap.containsKey(normKey)) {
			transformsMap.put(normKey, computeRotations(indZero, indReal));
		}
		return transformsMap.get(normKey);
	}

	public Set<Map.Entry<Integer, Integer>> populateNormalisations(int maxInd) {

		for (int indZero = 2; indZero <= maxInd; indZero++) {
			for (int indReal = 2; indReal <= maxInd; indReal++) {
				if (indReal % 2 == 1 && indReal == indZero) {
					continue;
				}
				Map.Entry<Integer, Integer> normKey = new AbstractMap.SimpleImmutableEntry<>(indZero, indReal);
				if (!transformsMap.containsKey(normKey)) {
					List<MomentTransform> rotations = computeRotations(indZero, indReal);
					if (rotations.size() > 0) {
						transformsMap.put(normKey, rotations);
					}
				}
			}
		}
		return transformsMap.keySet();
	}

	// average moment amplitudes of the equivalent solutions
	public List<Double> getInvariants(int normalizationOrder) {

		if (invariantsMap.containsKey(normalizationOrder)) {
			return invariantsMap.get(normalizationOrder);
		}

		List<MomentTransform> altSolutions = getNormalizationSolutions(normalizationOrder, 2);

		List<Double> zmInvariants = new ArrayList<>(500);

		int nSolutions = altSolutions.size();
		int nMoments = altSolutions.get(0).getFlatMoments().length;

		for (int iInvariant = 0; iInvariant < nMoments; iInvariant++) {
			double inv_amp = 0;
			for (MomentTransform altSolution : altSolutions) {
				inv_amp += altSolution.getFlatMoments()[iInvariant].abs();
			}
			zmInvariants.add(inv_amp / nSolutions);
		}

		invariantsMap.put(normalizationOrder, zmInvariants);
		return zmInvariants;
	}

	public double compareInvariants(InvariantNorm other, int indZero) {
		List<Double> invariantsThis = this.getInvariants(indZero);
		List<Double> invariantsOther = other.getInvariants(indZero);
		if (invariantsThis.size() != invariantsOther.size()) {
			throw new IllegalArgumentException("Moments must be calculated to the same order.");
		}

		double sumDiffs = 0;
		for (int k = 0; k < invariantsThis.size(); k++) {
			double diffAbs = Math.abs(invariantsThis.get(k) - invariantsOther.get(k));
			double sumAbs = invariantsThis.get(k) + invariantsOther.get(k) + 1;
			sumDiffs += diffAbs / sumAbs;
		}
		return sumDiffs;
	}

	// superpose using one normalisation
	public Matrix4d alignTo(InvariantNorm other) {

		List<InvariantNorm> zc = new ArrayList<>();

		zc.add(this);
		zc.add(other);

		List<Matrix3d> rotations = RotationAlignment.alignMultiple(zc);

		Matrix4d m1 = new Matrix4d(rotations.get(0), this.getCenter(), 1);
		Matrix4d m2 = new Matrix4d(rotations.get(1), other.getCenter(), 1);
		m1.invert();
		m2.mul(m1);
		return m2;
	}

	public Map<Map.Entry<Integer, Integer>, List<MomentTransform>> getTransformsMap() {
		return transformsMap;
	}

	public Vector3d getCenter() {
		return center;
	}

}
