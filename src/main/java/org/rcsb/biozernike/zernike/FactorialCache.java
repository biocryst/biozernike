package org.rcsb.biozernike.zernike;

import java.util.Arrays;

class FactorialCache {

	private static final int maxI = 40;
	private static final double[] factorial = new double[maxI * maxI];

	static {
		Arrays.fill(factorial, 1.0);
		for (int i = 1; i < maxI; i++) {
			factorial[i * maxI + i] = (double) i;
			for (int j = i + 1; j < maxI; j++) {
				factorial[i * maxI + j] = factorial[i * maxI + j - 1] * j;
			}
		}
	}

	static double getValue(int i, int j) {
		return factorial[i * maxI + j];
	}
}
