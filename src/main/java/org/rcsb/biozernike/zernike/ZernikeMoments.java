package org.rcsb.biozernike.zernike;

import org.rcsb.biozernike.complex.Complex;
import org.rcsb.biozernike.volume.Volume;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class ZernikeMoments {

	private static final Logger logger = LoggerFactory.getLogger(ZernikeMoments.class);

	private int maxOrder;
	private Volume volume;
	// zernike moments are calculated through geometric moments
	private GeometricMoments gm;
	// "Original" refers to rotation,
	// i.e. the moments correspond to the original orientation of the object
	// "(Un)scaled" refers to orthonormalization with c_{lm} coefficients.
	// "Unscaled" moments can be rotated, "scaled" moments can be compared.
	private List<List<List<Complex>>> originalMoments;
	private List<List<List<Complex>>> originalMomentsUnscaled;
	//	private NormalizationAlignment normalization;
	private ZernikeMoments() {
	}

	public ZernikeMoments(Volume volume, int maxOrder) {

		if (!volume.isNormalized()) {
			volume.normalize();
		}

		this.volume = volume;
		this.maxOrder = maxOrder;

		reset();
	}

	public ZernikeMoments(List<List<List<Complex>>> moments, boolean isMomentsOrthonormal) {
		this.maxOrder = moments.size()-1;

		List<List<List<Complex>>> momentsMult = scaleMoments(moments, isMomentsOrthonormal);
		if (isMomentsOrthonormal) {
			this.originalMoments = moments;
			this.originalMomentsUnscaled = momentsMult;
		} else {
			this.originalMoments = momentsMult;
			this.originalMomentsUnscaled = moments;
		}
	}

	public Volume getVolume() {
		return volume;
	}

	private void reset() {
		this.gm = new GeometricMoments(volume, volume.getRadiusVarVolume(), maxOrder);
		this.originalMomentsUnscaled = new ArrayList<>(maxOrder + 1);
		this.originalMoments = new ArrayList<>(maxOrder + 1);
		computeMoments();
//		this.normalization = new NormalizationAlignment(this, volume.getCenterReal());
	}

	private void computeMoments() {

		for (int n = 0; n <= maxOrder; ++n) {
			List<List<Complex>> zmLLevelUnscaled = new ArrayList<>(n / 2 + 1);
			List<List<Complex>> zmLLevelScaled = new ArrayList<>(n / 2 + 1);
			int l0 = n % 2, li = 0;
			for (int l = l0; l <= n; ++li, l += 2) {
				List<Complex> zmMLevelUnscaled = new ArrayList<>(l + 1);
				List<Complex> zmMLevelScaled = new ArrayList<>(l + 1);
				for (int m = 0; m <= l; ++m) {
					Complex zm = new Complex(0, 0);
					List<ComplexCoeff> gCoeffsNLM = ZernikeCache.getGCoefs(n, li, m);
					for (ComplexCoeff cc : gCoeffsNLM) {
						zm = zm.add(cc.c.mul(gm.getMoment(cc.p, cc.q, cc.r)));
					}
					zm = zm.mul(3.0 / (4.0 * Math.PI));
					zmMLevelUnscaled.add(zm);
					zmMLevelScaled.add(zm.mul(ZernikeCache.getClmValue(l, m)));
				}
				zmLLevelUnscaled.add(zmMLevelUnscaled);
				zmLLevelScaled.add(zmMLevelScaled);
			}
			originalMomentsUnscaled.add(zmLLevelUnscaled);
			originalMoments.add(zmLLevelScaled);
		}
	}

	public Complex getOriginalUnscaled(int n, int l, int m) {
		return originalMomentsUnscaled.get(n).get(l).get(m);
	}

	public Complex getMoment (int n, int l, int m)
	{
		if (m >= 0) {
			return originalMoments.get(n).get(l).get(m);
		} else {
			Complex moment = originalMoments.get(n).get(l).get(-m).conj();
			if (m%2 != 0) {
				moment = moment.negate();
			}
			return moment;
		}
	}


	public List<List<List<Complex>>> getOriginalMoments() {
		return originalMoments;
	}

	public List<List<List<Complex>>> getOriginalMomentsUnscaled() {
		return originalMomentsUnscaled;
	}

	public int getMaxOrder() {
		return maxOrder;
	}

	public int setMaxOrder(int maxOrderNew) {

		if (maxOrderNew < maxOrder) {
			//trim existing layers of moments
			originalMoments.subList(maxOrderNew + 1, originalMoments.size()).clear();
			originalMomentsUnscaled.subList(maxOrderNew + 1, originalMomentsUnscaled.size()).clear();
			maxOrder = maxOrderNew;
		}

		if (maxOrderNew > maxOrder) {
			// TODO: calc only additional layers in gm and z moments
			logger.warn("Increasing max order is inefficient in the current implementation and should be avoided. Go in the other direction.");
			maxOrder = maxOrderNew;
			reset();
		}

		return originalMoments.stream().flatMap(List::stream).mapToInt(List::size).sum();
	}

	//TODO:
	// 1) move to an util class
	// 2) error checks
	// 3) derive max order from input
	// 4) custom iterator?
	//
	public static List<List<List<Complex>>> scaleMoments(List<List<List<Complex>>> moments, boolean isMomentsOrthonormal) {
		List<List<List<Complex>>> scaledMoments = new ArrayList<>(moments.size());
		int maxOrder = moments.size()-1;
		for (int n = 0; n <= maxOrder; ++n) {
			List<List<Complex>> zmLLevelScaled = new ArrayList<>(n / 2 + 1);
			int l0 = n % 2, li = 0;
			for (int l = l0; l <= n; ++li, l += 2) {
				List<Complex> zmMLevelScaled = new ArrayList<>(l + 1);
				for (int m = 0; m <= l; ++m) {
					Complex zm = moments.get(n).get(li).get(m);
					double clmCoef = ZernikeCache.getClmValue(l, m);
					if (isMomentsOrthonormal) {
						clmCoef = 1/clmCoef;
					}
					zmMLevelScaled.add(zm.mul(clmCoef));
				}
				zmLLevelScaled.add(zmMLevelScaled);
			}
			scaledMoments.add(zmLLevelScaled);
		}
		return scaledMoments;
	}

	public static List<Complex> flattenMomentsComplex(List<List<List<Complex>>> hierarchicalMoments) {
		return hierarchicalMoments.stream().
				flatMap(List::stream).
				flatMap(List::stream).
				collect(Collectors.toList());
	}

	public static List<List<List<Complex>>> unFlattenMomentsComplex(List<Complex> flatMoments, int maxOrder) {
		List<List<List<Complex>>> hierarchicalMoments = new ArrayList<>();
		int flatInd = 0;
		for (int n = 0; n <= maxOrder; ++n) {
			List<List<Complex>> zmLLevel = new ArrayList<>(n / 2 + 1);
			int l0 = n % 2, li = 0;
			for (int l = l0; l <= n; ++li, l += 2) {
				List<Complex> zmMLevel = new ArrayList<>(l + 1);
				for (int m = 0; m <= l; ++m) {
					zmMLevel.add(flatMoments.get(flatInd++));
				}
				zmLLevel.add(zmMLevel);
			}
			hierarchicalMoments.add(zmLLevel);
		}
		return hierarchicalMoments;
	}


	public static List<Double> flattenMomentsDouble(List<List<List<Complex>>> hierarchicalMoments) {
		return flattenMomentsComplex(hierarchicalMoments).stream().
				map(r -> Arrays.asList(r.getReal(),r.getImaginary())).
				flatMap(List::stream).
				collect(Collectors.toList());
	}

	public static List<List<List<Complex>>> unFlattenMomentsDouble(List<Double> flatMoments, int maxOrder) {
		List<Complex> flatMomentsComplex =
				IntStream.range(0, flatMoments.size()/2).
				mapToObj( i -> new Complex(flatMoments.get(2*i),flatMoments.get(2*i+1))).
				collect(Collectors.toList());

		return unFlattenMomentsComplex(flatMomentsComplex, maxOrder);
	}

}
