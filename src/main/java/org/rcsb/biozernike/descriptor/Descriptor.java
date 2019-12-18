package org.rcsb.biozernike.descriptor;

import org.rcsb.biozernike.InvariantNorm;
import org.rcsb.biozernike.volume.Volume;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.*;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.rcsb.biozernike.zernike.ZernikeMoments;

import javax.vecmath.Point3d;
import java.util.*;

public class Descriptor {
	private DescriptorConfig config;

	private List<Double> geometryDescriptor = null;
	private List<Double> momentDescriptor = null;

	private List<Double> momentsAlign = null;
	private List<List<Double>> momentInvariantsRaw = null;

	private List<Double> volumeCenter = null;

	public Descriptor(Point3d[] reprPoints, String[] resNames, DescriptorConfig config) {
		this.config = config;

		Volume volume = new Volume();
		volume.create(reprPoints, resNames);

		InvariantNorm normalization = new InvariantNorm(volume, config.maxOrderZernike);


		if(config.mode.contains(DescriptorMode.CALCULATE_RAW)) {
			calcGeometryDescriptor(volume, reprPoints);
			calcMomentInvariantsRaw(normalization);
		}

		if(config.mode.contains(DescriptorMode.ALIGN)) {
			momentsAlign = ZernikeMoments.flattenMomentsDouble(
					normalization.getMoments().
							getOriginalMomentsUnscaled().
							subList(0, config.maxOrderZernikeAlign)
			);
			double[] center = volume.getCenterReal();
			volumeCenter = new ArrayList<Double>() {{add(center[0]);add(center[1]);add(center[2]);}};
		}

		if(config.mode.contains(DescriptorMode.CALCULATE_PROCESSED)) {
			calcMomentDescriptor();
		}

	}

	private void calcGeometryDescriptor(Volume volume, Point3d[] reprPoints) {

		geometryDescriptor = new ArrayList<>();
		geometryDescriptor.add(volume.getRadiusVarReal());
		geometryDescriptor.add(volume.getResiduesNominalWeight());

		Point3d centerPoint = new Point3d(0, 0, 0);

		for (Point3d selPoint : reprPoints) {
			centerPoint.add(selPoint);
		}
		centerPoint.scale(1 / (double) reprPoints.length);

		double[] distances = new double[reprPoints.length];
//		double[][] centeredPoints = new double[reprPoints.length][3];

		for (int iPoint = 0; iPoint < reprPoints.length; iPoint++) {
			distances[iPoint] = centerPoint.distance(reprPoints[iPoint]);
//			Point3d centeredPoint = new Point3d(reprPoints[iPoint]);
//			centeredPoint.sub(centerPoint);
//			centeredPoints[iPoint][0] = centeredPoint.x;
//			centeredPoints[iPoint][1] = centeredPoint.y;
//			centeredPoints[iPoint][2] = centeredPoint.z;
		}

		//TODO: make this optional
//		Covariance covariance = new Covariance(centeredPoints);
//		RealMatrix cov = covariance.getCovarianceMatrix();
//		EigenDecomposition eigenDecomposition = new EigenDecomposition(cov);
//		double[] eigenvalues = eigenDecomposition.getRealEigenvalues();
//
//		geometryDescriptor.add(Math.sqrt(eigenvalues[0]));
//		geometryDescriptor.add(Math.sqrt(eigenvalues[1]));
//		geometryDescriptor.add(Math.sqrt(eigenvalues[2]));

		StandardDeviation standardDeviation = new StandardDeviation();
		Skewness skewness = new Skewness();
		Kurtosis kurtosis = new Kurtosis();
		Percentile percentile = new Percentile();
		percentile.setData(distances);

		for (double p = 10; p < 100; p += 10) {
			geometryDescriptor.add(percentile.evaluate(p));
		}

		geometryDescriptor.add(standardDeviation.evaluate(distances));
		geometryDescriptor.add(skewness.evaluate(distances));
		geometryDescriptor.add(kurtosis.evaluate(distances));

	}

	private void calcMomentDescriptor() {
		momentDescriptor = new ArrayList<>();

		for (int i = 0; i < config.normOrders.length; i++) {
			List<Double> normedMoments = momentInvariantsRaw.get(i);
			int[] selIndices = config.indicesZernike.get(i);
			for (int selIndex : selIndices) {
				momentDescriptor.add(normedMoments.get(selIndex));
			}
		}

	}

	private void calcMomentInvariantsRaw(InvariantNorm normalization) {
		momentInvariantsRaw = new ArrayList<>();

		for (int normOrder : config.normOrders) {
			momentInvariantsRaw.add(new ArrayList<>(normalization.getInvariants(normOrder)));
		}
	}

	public List<Double> getGeometryDescriptor() {
		return geometryDescriptor;
	}

	public List<Double> getMomentDescriptor() {
		return momentDescriptor;
	}

	public List<List<Double>> getMomentInvariantsRaw() {
		return momentInvariantsRaw;
	}

	public List<Double> getVolumeCenter() {
		return volumeCenter;
	}

	public List<Double> getMomentsAlign() {
		return momentsAlign;
	}

}


