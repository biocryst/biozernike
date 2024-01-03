package org.rcsb.biozernike.zernike;

import org.rcsb.biozernike.complex.Complex;
import org.rcsb.biozernike.volume.Volume;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static java.lang.Math.abs;
import static java.lang.Math.pow;

/**
 * Class to hold the 3D Zernike moments of a given volume.
 *
 */
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

	/**
	 * Get the 3D Zernike moment with index n, l, m (Ω_nl^m).
	 * Note that negative m indices are calculated from the positive ones, as in
	 * equation 9 of Novotni and Klein 2004.
	 * @param n the n index
	 * @param l the l index
	 * @param m the m index
	 * @return the Ω_nl^m
	 */
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

	/**
	 * Get the original 3D Zernike moments (on indices n, l, m).
	 * "Original" refers to rotation, i.e. the moments correspond to the original orientation of the object.
	 * <p>
	 * Note that the negative m indices are omitted from the list. They are obtainable with {@link #getMoment(int, int, int)}
	 * @return the Zernike moments
	 */
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


	public static Volume reconstructVolume(ZernikeMoments zm, int dim, int _maxN, boolean useCache, boolean saveCache) throws Exception {
		int dimX = dim;
		int dimY = dim;
		int dimZ = dim;
		int volume_size = dimX*dimY*dimZ;
		double center = dim/2.0;
		//scaling
		double scale = 2.0/dim;

		boolean cache_available = false;

		//translation
		double vx = center;
		double vy = center;
		double vz = center;

		double[] point = new double[3];

		double dist_threshold = 1;
		double[] _grid = new double[dimX*dimY*dimZ];

		if (useCache) {
			Map<Integer, Map<Integer, Map<Integer, Complex[]>>> zpN = ReconstructionCache.readZpCache();
			for (int n=0; n<=_maxN; ++n) {
				int maxK = n / 2;
				for (int k = 0; k <= maxK; ++k) {
					int l = n - 2 * k;
					int li = l / 2;
					for (int m = -l; m <= l; ++m) {
						Complex[] zp_xyz = zpN.get(n).get(l).get(m);
						Complex mom = zm.getMoment(n, li, m);
						for (int flat_ind = 0; flat_ind < volume_size; flat_ind++) {
							_grid[flat_ind] += mom.mul(zp_xyz[flat_ind]).getReal();
						}
					}
				}
			}
		} else {

			Map<Integer, Map<Integer, Map<Integer, Complex[]>>> zpN = new HashMap<>();
			for (int n = 0; n <= _maxN; ++n) {

				Map<Integer, Map<Integer, Complex[]>> zpL = new HashMap<>();

				int maxK = n / 2;
				for (int k = 0; k <= maxK; ++k) {

					int l = n - 2 * k;
					int li = l / 2;

					Map<Integer, Complex[]> zpM = new HashMap<>();
					for (int m = -l; m <= l; ++m) {
						System.out.println("N=" + n + ", L=" + l + " (Li=" + li + "), M=" + m);

						int absM = abs(m);

						List<ComplexCoeff> gCoeffsNLM = ZernikeCache.getGCoefs(n, li, absM);

						int nCoeffs = gCoeffsNLM.size();

						Complex[] zp_xyz = new Complex[volume_size];
						Arrays.fill(zp_xyz, new Complex(0,0));

						for (int x = 0; x < dimX; ++x) {
							point[0] = (vx - x) * scale; // TODO: figure out why x is reversed
							double p0sq = point[0] * point[0];
							if (p0sq > dist_threshold) {
								continue;
							}

							for (int y = 0; y < dimY; ++y) {
								point[1] = (y - vy) * scale;
								double p1sq = point[1] * point[1];
								if (p0sq + p1sq > dist_threshold) {
									continue;
								}

								for (int z = 0; z < dimZ; ++z) {
									// the origin is in the middle of the grid, all voxels are
									// projected into the unit ball
									point[2] = (z - vz) * scale;

									if (p0sq + p1sq + point[2] * point[2] > dist_threshold) {
										continue;
									}

									// zernike polynomial evaluated at point
									Complex zp = new Complex(0, 0);

									for (int i = 0; i < nCoeffs; ++i) {

										ComplexCoeff cc = gCoeffsNLM.get(i);
										double clmCoef = ZernikeCache.getClmValue(l, absM);
										Complex cvalue = cc.c.mul(clmCoef);

										// conjugate if m negative
										if (m < 0) {
											cvalue = cvalue.conj();

											// take care of the sign
											if (m % 2 != 0) {
												cvalue = cvalue.negate();
											}
										}

										zp = zp.add(cvalue.mul(pow(point[0], cc.p)).mul(pow(point[1], cc.q)).mul(pow(point[2], cc.r)));

									}
									zp_xyz[(z * dimY + y) * dimX + x] = zp;
									_grid[(z * dimY + y) * dimX + x] += zp.mul(zm.getMoment(n, li, m)).getReal();
								} //z
							} //y
						} // x

						zpM.put(m, zp_xyz);
					} // m level
					zpL.put(l, zpM);
				} // l level
				zpN.put(n, zpL);
			} // n level
			if (saveCache) {
				ReconstructionCache.writeZpCache(zpN);
			}
		}

		Volume volume = new Volume();
		int[] dimensions = new int[]{dimX, dimY, dimZ};
		volume.createFromData(dimensions, _grid,0.8);
		return volume;
	}
}
