package org.rcsb.biozernike.zernike;

import org.rcsb.biozernike.complex.Complex;

public class ComplexCoeff {
	public int p, q, r;
	public Complex c;

	ComplexCoeff(int p, int q, int r, Complex c) {
		this.p = p;
		this.q = q;
		this.r = r;
		this.c = c;
	}
}
