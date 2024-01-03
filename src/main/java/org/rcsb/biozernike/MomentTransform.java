package org.rcsb.biozernike;

import org.rcsb.biozernike.complex.Complex;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;
import java.io.Serializable;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class MomentTransform implements Serializable {
	private Complex a, b;
	private Complex[] flatMoments;
	private Complex[] flatMomentsDouble;
	private Matrix3d R = null;
	private List<List<List<Complex>>> moments;
	public MomentTransform() {
		flatMoments = null;
	}

	public void setMoments(List<List<List<Complex>>> moments) {
		this.moments = moments;
		this.flatMoments = moments.stream().
				flatMap(List::stream).
				flatMap(List::stream).toArray(Complex[]::new);
	}

	public Complex[] getFlatMoments() {
		return flatMoments;
	}

	public void setFlatMoments(Complex[] moments) {
		this.flatMoments = moments;
	}

	public void trimMoments(int nRemain) {
		Complex[] flatMomentsTrimmed = new Complex[nRemain];
		System.arraycopy(flatMoments, 0, flatMomentsTrimmed, 0, nRemain);
		flatMoments = flatMomentsTrimmed;
	}

	public void setRotation(Complex a, Complex b) {
		this.a = a;
		this.b = b;
	}

	public double distanceTo(MomentTransform other) {
		double sumDiffs = 0;
		Complex[] flatMomentsOther = other.flatMoments;
		for (int indMoment = 0; indMoment < flatMoments.length; indMoment++) {
			double diffAbs = flatMoments[indMoment].subtract(flatMomentsOther[indMoment]).abs();
			double sumAbs = flatMoments[indMoment].abs() + flatMomentsOther[indMoment].abs() + 1;
			sumDiffs += diffAbs / sumAbs;
		}
		return sumDiffs/(flatMoments.length-5);
	}

	public Matrix4d getCoordTransform(Vector3d translation) {
		return new Matrix4d(rotation(), translation, 1);
	}

	public Matrix3d rotation() {
		if (R == null) {

			// See Canterakis 1996, last equation page 6. See also http://mathworld.wolfram.com/Cayley-KleinParameters.html
			Complex a2pb2 = a.pow(2).add(b.pow(2));
			Complex a2mb2 = a.pow(2).subtract(b.pow(2));

			R = new Matrix3d(
					a2pb2.getReal(),
					-a2mb2.getImaginary(),
					2 * a.mul(b).getImaginary(),
					a2pb2.getImaginary(),
					a2mb2.getReal(),
					-2 * a.mul(b).getReal(),
					2 * a.mul(b.conj()).getImaginary(),
					2 * a.mul(b.conj()).getReal(),
					a.mul(a.conj()).getReal() - b.mul(b.conj()).getReal()
			);
		}

		return R;
	}

	public Complex getA() {
		return a;
	}

	public Complex getB() {
		return b;
	}

	public List<List<List<Complex>>> getMoments() {
		return moments;
	}
}
