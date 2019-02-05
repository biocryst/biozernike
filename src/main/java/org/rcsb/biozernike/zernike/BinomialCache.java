package org.rcsb.biozernike.zernike;

public class BinomialCache {

	private static int maxI = 60;
	private static final double[] binomial = new double[maxI * maxI];

	static {
		for (int i = 0; i < maxI; i++) {
			for (int j = 0; j < maxI + 1 - i; j++) {
				if (i == 0 || j == 0) {
					binomial[i * maxI + j] = 1.0;
				} else {
					binomial[i * maxI + j] = binomial[i * maxI + j - 1] + binomial[(i - 1) * maxI + j];
				}
			}
		}
	}

	public static double getValue(int i, int j) {
		return binomial[j * maxI + i - j];
	}
}
