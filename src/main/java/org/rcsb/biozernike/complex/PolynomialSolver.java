package org.rcsb.biozernike.complex;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * An implementation of the Durand-Kerner-Weierstrass method to
 * determine the roots of a complex univariate polynomial, as described
 * <a href="http://en.wikipedia.org/wiki/Durand-Kerner_method">here</a>.
 *
 * @author John B. Matthews; distribution per LGPL.
 */

public class PolynomialSolver {

	private static final Logger logger = LoggerFactory.getLogger(PolynomialSolver.class);
	private final int MAX_COUNT = 999;
	private double epsilon = 1E-15;

	/**
	 * This implementation uses the Durand-Kerner-Weierstrass method
	 * to find the roots of a polynomial with complex coefficients.
	 * The method requires a monic polynomial; some error may occur
	 * when dividing by the leading coefficient.
	 * The voxelArray should have the highest order coefficient first.
	 *
	 * @param ca a variable number of Complex polynomial coefficients.
	 * @return an voxelArray of the Complex roots of the polynomial.
	 */
	public Complex[] roots(Complex[] ca) {

		Complex[] caOrig = new Complex[ca.length];
		System.arraycopy(ca, 0, caOrig, 0, ca.length);

		Complex[] a0 = new Complex[ca.length - 1];
		Complex[] a1 = new Complex[ca.length - 1];

		// Divide by leading coefficient if not monic
		Complex leading = ca[0];
		if (ca[0].getReal() != 1 || ca[0].getImaginary() != 0) {
			for (int i = 0; i < ca.length; i++) {
				ca[i] = ca[i].divide(leading);
			}
		}

		// Initialize a0
		Complex result = new Complex(0.4, 0.9);
		a0[0] = new Complex(1.0, 0.0);
		for (int i = 1; i < a0.length; i++) {
			a0[i] = a0[i - 1].mul(result);
		}

		// Iterate
		int count = 0;
		while (true) {
			for (int i = 0; i < a0.length; i++) {
				result = new Complex(1.0, 0.0);
				for (int j = 0; j < a0.length; j++) {
					if (i != j) {
						result = a0[i].subtract(a0[j]).mul(result);
					}
				}
				a1[i] = a0[i].subtract(eval(ca, a0[i]).divide(result));
			}
			count++;
			if (done(a0, a1)) {
				break;
			}

			if (count > MAX_COUNT) {
				String coefs = "";
				for (int i = 0; i < caOrig.length; i++) {
					coefs += caOrig[i].toString() + ", ";
				}
				String sol = "";
				for (int i = 0; i < a1.length; i++) {
					sol += a1[i].toString() + ", ";
				}
				logger.warn("Maximum number of iterations reached. The solution may not be accurate. \nCoefs: " + coefs + "\nSol: " + sol);
				break;
			}

			System.arraycopy(a1, 0, a0, 0, a1.length); // a0 := a1
		}

		return a1;
	}

	// Determine if the components of two arrays are unchanging
	private boolean done(Complex[] a, Complex[] b) {
		boolean unchanged = true;
		Complex delta;
		for (int i = 0; i < a.length; i++) {
			delta = a[i].subtract(b[i]);
			unchanged &= Math.abs(delta.getReal()) < epsilon
					&& Math.abs(delta.getImaginary()) < epsilon;
		}
		return unchanged;
	}

/**
 * Evaluate the polynomial at x using
 * <a href="http://en.wikipedia.org/wiki/Horner_scheme">Horner's scheme</a>.
 * The voxelArray should have the highest order coefficient first.
 *
 * @param  ca an voxelArray of Complex polynomial coefficients.
 * @param  x the function argument.
 * @return the Complex value of the function at x.
*/
	private Complex eval(Complex[] ca, Complex x) {
		Complex result = ca[0];
		for (int i = 1; i < ca.length; i++) {
			result = result.mul(x).add(ca[i]);
		}
		return result;
	}

	/** Return the nominal tolerance value. */
	public double getEpsilon() {
		return epsilon;
	}

	/** Set the nominal tolerance value. */
	public void setEpsilon(double e) {
		epsilon = e;
	}
}

