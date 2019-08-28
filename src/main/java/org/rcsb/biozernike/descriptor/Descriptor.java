package org.rcsb.biozernike.descriptor;

import org.rcsb.biozernike.InvariantNorm;
import org.rcsb.biozernike.MomentTransform;
import org.rcsb.biozernike.volume.Volume;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.*;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import javax.vecmath.Point3d;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class Descriptor {
	private DescriptorConfig config;

	private List<Double> geometryInvariants;
	private List<Double> momentInvariantsSel;
	private List<Double> momentInvariantsAll;

	private Map<Map.Entry<Integer, Integer>, List<MomentTransform>> transformsMap;
	private double[] volumeCenter;

	public Descriptor(Point3d[] reprPoints, String[] resNames, DescriptorConfig config) {
		this.config = config;

		Volume volume = new Volume();
		volume.create(reprPoints, resNames);

		InvariantNorm normalization = new InvariantNorm(volume, config.maxOrderInvariants);
		volumeCenter = volume.getCenterReal();
		calcGeomInvariants(volume, reprPoints);
		calcMomentInvariants(normalization);
		transformsMap = calcAlignments(normalization);

	}

	private void calcGeomInvariants(Volume volume, Point3d[] reprPoints) {

		geometryInvariants = new ArrayList<>();
		geometryInvariants.add(volume.getRadiusVarReal());
		geometryInvariants.add(volume.getResiduesNominalWeight());

		Point3d centerPoint = new Point3d(0, 0, 0);

		for (Point3d selPoint : reprPoints) {
			centerPoint.add(selPoint);
		}
		centerPoint.scale(1 / (double) reprPoints.length);

		double[] distances = new double[reprPoints.length];
		double[][] centeredPoints = new double[reprPoints.length][3];

		for (int iPoint = 0; iPoint < reprPoints.length; iPoint++) {
			distances[iPoint] = centerPoint.distance(reprPoints[iPoint]);
			Point3d centeredPoint = new Point3d(reprPoints[iPoint]);
			centeredPoint.sub(centerPoint);
			centeredPoints[iPoint][0] = centeredPoint.x;
			centeredPoints[iPoint][1] = centeredPoint.y;
			centeredPoints[iPoint][2] = centeredPoint.z;
		}

		Covariance covariance = new Covariance(centeredPoints);
		RealMatrix cov = covariance.getCovarianceMatrix();
		EigenDecomposition eigenDecomposition = new EigenDecomposition(cov);
		double[] eigenvalues = eigenDecomposition.getRealEigenvalues();

		geometryInvariants.add(Math.sqrt(eigenvalues[0]));
		geometryInvariants.add(Math.sqrt(eigenvalues[1]));
		geometryInvariants.add(Math.sqrt(eigenvalues[2]));


		StandardDeviation standardDeviation = new StandardDeviation();
		Skewness skewness = new Skewness();
		Kurtosis kurtosis = new Kurtosis();
		Percentile percentile = new Percentile();
		percentile.setData(distances);

		for (double p = 10; p < 100; p += 10) {
			geometryInvariants.add(percentile.evaluate(p));
		}

		geometryInvariants.add(standardDeviation.evaluate(distances));
		geometryInvariants.add(skewness.evaluate(distances));
		geometryInvariants.add(kurtosis.evaluate(distances));

	}


	private void calcMomentInvariants(InvariantNorm normalization) {
		momentInvariantsSel = new ArrayList<>();

		momentInvariantsAll = new ArrayList<>(normalization.getFingerprint());

		for (int normOrder : config.normOrders) {
			momentInvariantsAll.addAll(normalization.getInvariants(normOrder));
		}

		for (int i = 0; i < config.searchIndicesZernike.length; i++) {
			momentInvariantsSel.add(momentInvariantsAll.get(config.searchIndicesZernike[i]));
		}
	}

	private Map<Map.Entry<Integer, Integer>, List<MomentTransform>> calcAlignments(InvariantNorm normalization) {
		normalization.setMaxOrder(config.maxOrderAlignment);

		List<Map.Entry<Integer, Integer>> normKeys = new ArrayList<>();
		normKeys.add(new AbstractMap.SimpleImmutableEntry<>(2, 2));
		normKeys.add(new AbstractMap.SimpleImmutableEntry<>(3, 5));
		normKeys.add(new AbstractMap.SimpleImmutableEntry<>(4, 4));
		normKeys.add(new AbstractMap.SimpleImmutableEntry<>(6, 6));

		normalization.setNormalisations(normKeys);
		return normalization.getTransformsMap();
	}

	public List<Double> getGeomInvariants() {
		return geometryInvariants;
	}

	public List<Double> getMomentInvariantsSel() {
		return momentInvariantsSel;
	}
	public List<Double> getMomentInvariantsAll() {
		return momentInvariantsAll;
	}

	public Map<Map.Entry<Integer, Integer>, List<MomentTransform>> getTransformsMap() {
		return transformsMap;
	}

	public double[] getVolumeCenter() {
		return volumeCenter;
	}

}


