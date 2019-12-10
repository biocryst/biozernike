package org.rcsb.biozernike.zernike;

import org.rcsb.biozernike.complex.Complex;
import org.rcsb.biozernike.volume.Volume;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.ArrayList;
import java.util.List;


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


}
