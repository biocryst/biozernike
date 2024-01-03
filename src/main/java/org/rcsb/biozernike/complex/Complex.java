package org.rcsb.biozernike.complex;

import java.io.Serializable;

/**
 * Simple class for Complex numbers.
 *
 */

public class Complex implements Serializable {
	private final double real;
	private final double imaginary;

	public Complex(double real, double imaginary) {
		this.real = real;
		this.imaginary = imaginary;
	}

	public double getReal() {
		return real;
	}

	public double getImaginary() {
		return imaginary;
	}

	public double abs() {
		return Math.sqrt(norm());
	}

	public double norm() {
		return real * real + imaginary * imaginary;
	}

	public Complex conj() {
		return new Complex(real, -imaginary);
	}

	public Complex negate() {
		return new Complex(-real, -imaginary);
	}

	public Complex add(Complex other) {
		return new Complex(this.real + other.real,
				this.imaginary + other.imaginary);
	}

	public Complex subtract(Complex other) {
		return new Complex(this.real - other.real,
				this.imaginary - other.imaginary);
	}

	public Complex mul(Complex other) {
		return new Complex(this.real * other.real - this.imaginary * other.imaginary,
				this.real * other.imaginary + this.imaginary * other.real);
	}

	public Complex mul(double factor) {
		return new Complex(real * factor, imaginary * factor);
	}


	public Complex reciprocal() {
		double denom = real * real + imaginary * imaginary;
		return new Complex(real / denom, -imaginary / denom);
	}

	public Complex divide(Complex other) {
		return this.mul(other.reciprocal());
	}

	/**
	 * Power function (integers only)
	 * implemented after https://stackoverflow.com/a/26691276
	 * @param n
	 * @return
	 */
	public Complex pow(int n) {
		if (n == 0) {
			return new Complex(1.0, 0);
		}
		if (n < 0) {
			return new Complex(1, 0).divide(pow(-n));
		} else {
			// Positive power
			Complex powerOfHalfN = pow(n / 2);
			if (n % 2 == 1) {
				// Odd n
				return this.mul(powerOfHalfN).mul(powerOfHalfN);
			} else {
				// Even n
				return powerOfHalfN.mul(powerOfHalfN);
			}
		}
	}

	public Complex sqrt() {
		double r=Math.sqrt(this.abs());
		double theta=Math.atan2(imaginary,real)/2;
		return new Complex(r*Math.cos(theta),r*Math.sin(theta));
	}

	@Override
	public String toString() {
		return real+", "+imaginary;
	}

	@Override
	public boolean equals(Object o) {
		if (!(o instanceof Complex)) {
			return false;
		}
		Complex c = (Complex) o;
		return Double.compare(this.real,c.real) == 0  &&
				Double.compare(this.imaginary,c.imaginary) == 0;
	}
}
