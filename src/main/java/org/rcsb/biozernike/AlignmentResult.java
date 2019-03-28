package org.rcsb.biozernike;

import javax.vecmath.Matrix4d;
import java.util.*;

public class AlignmentResult {
	private List<Matrix4d> transforms = new ArrayList<>();

	private int indZero = 0;
	private int indReal = 0;
	private double score = Double.POSITIVE_INFINITY;

	public AlignmentResult() {}

	public AlignmentResult(List<Matrix4d> transforms, Map.Entry<Integer, Integer> normKey, double score) {
		this.transforms = transforms;
		this.indZero = normKey.getKey();
		this.indReal = normKey.getValue();
		this.score = score;
	}

	public List<Matrix4d> getTransforms() {
		return transforms;
	}

	public void setTransforms(List<Matrix4d> transforms) {
		this.transforms = transforms;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}

	public int getIndZero() {
		return indZero;
	}

	public void setIndZero(int indZero) {
		this.indZero = indZero;
	}

	public int getIndReal() {
		return indReal;
	}

	public void setIndReal(int indReal) {
		this.indReal = indReal;
	}

}
