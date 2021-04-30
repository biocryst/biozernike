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

/**
 * A class to hold BioZernike descriptors data as described in Guzenko et al 2020.
 *
 * @author Dmytro Guzenko
 */
public class Descriptor {

	private static final double[] PERCENTILES_FOR_GEOM = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0};

	private final DescriptorConfig config;

	private double[] geometryDescriptor;
	private double[] momentDescriptor;

	private List<Double> momentsAlign;
	private List<List<Double>> momentInvariantsRaw;

	private List<Double> volumeCenter;

	/**
	 * Construct a new Descriptor object given coordinate data and a configuration
	 * @param reprPoints the representative atom coordinates
	 * @param resNames the residue names for each representative point
	 * @param config the configuration
	 */
	public Descriptor(Point3d[] reprPoints, String[] resNames, DescriptorConfig config) {
		this.config = config;

		Volume volume = new Volume();
		volume.create(reprPoints, resNames);

		init(volume, reprPoints);
	}

	/**
	 * Construct a new Descriptor object given volume data and a configuration.
	 * Only moment descriptors will be available. No geometric descriptors will
	 * be available through {@link #getGeometryDescriptor()}
	 * @param volume the volume
	 * @param config the configuration
	 */
	public Descriptor(Volume volume, DescriptorConfig config) {
		this.config = config;

		init(volume, null);
	}

	/**
	 * Initialise the Descriptor object
	 * @param volume the volume
	 * @param reprPoints the points, if null no geometric descriptors are calculated
	 */
	private void init(Volume volume, Point3d[] reprPoints) {
		InvariantNorm normalization = new InvariantNorm(volume, config.maxOrderZernike);

		if (config.mode.contains(DescriptorMode.CALCULATE_RAW)) {
			if (reprPoints!=null) {
				calcGeometryDescriptor(volume, reprPoints, config.withCovEigenValsInGeom);
			}
			calcMomentInvariantsRaw(normalization);
		}

		if (config.mode.contains(DescriptorMode.ALIGN)) {
			momentsAlign = ZernikeMoments.flattenMomentsDouble(
					normalization.getMoments().
							getOriginalMomentsUnscaled().
							subList(0, config.maxOrderZernikeAlign+1)
			);
			double[] center = volume.getCenterReal();
			volumeCenter = new ArrayList<Double>() {{add(center[0]);add(center[1]);add(center[2]);}};
		}

		if (config.mode.contains(DescriptorMode.CALCULATE_PROCESSED)) {
			calcMomentDescriptor();
		}
	}

	/**
	 * Calculate the BioZernike geometry descriptors. Subsequently call {@link #getGeometryDescriptor()} to retrieve them.
 	 * This results in an array of size 2 + length of{@link #PERCENTILES_FOR_GEOM} + 3 [+ 3, if {@link DescriptorConfig#withCovEigenValsInGeom} is true],
	 * with values:
	 * <li>radius</li>
	 * <li>residues nominal weight</li>
	 * <li>percentiles ({@link #PERCENTILES_FOR_GEOM}) of distance to centroid distribution</li>
	 * <li>standard deviation of distance to centroid distribution</li>
	 * <li>skewness of distance to centroid distribution</li>
	 * <li>kurtosis of distance to centroid distribution</li>
	 * <li>optionally (if {@link DescriptorConfig#withCovEigenValsInGeom} is true) 3 eigenvalues of distances to centroid covariance matrix</li>
	 * @param volume the volume
	 * @param reprPoints the points
	 * @param withCovarianceEigenvalues whether to calculate covariance eigen values (that will be the last 3 descriptors in array) or not
	 *                                  Note that covariance eigenvalue calculation is very expensive compare to the other geometric descriptors
	 */
	private void calcGeometryDescriptor(Volume volume, Point3d[] reprPoints, boolean withCovarianceEigenvalues) {

		List<Double> geomDescriptorList = new ArrayList<>();
		geomDescriptorList.add(volume.getRadiusVarReal());
		geomDescriptorList.add(volume.getResiduesNominalWeight());

		Point3d centerPoint = new Point3d(0, 0, 0);

		for (Point3d selPoint : reprPoints) {
			centerPoint.add(selPoint);
		}
		centerPoint.scale(1 / (double) reprPoints.length);

		double[] distances = new double[reprPoints.length];
		for (int iPoint = 0; iPoint < reprPoints.length; iPoint++) {
			distances[iPoint] = centerPoint.distance(reprPoints[iPoint]);
		}

		StandardDeviation standardDeviation = new StandardDeviation();
		Skewness skewness = new Skewness();
		Kurtosis kurtosis = new Kurtosis();
		Percentile percentile = new Percentile();
		percentile.setData(distances);

		for (double p : PERCENTILES_FOR_GEOM) {
			geomDescriptorList.add(percentile.evaluate(p));
		}

		geomDescriptorList.add(standardDeviation.evaluate(distances));
		geomDescriptorList.add(skewness.evaluate(distances));
		geomDescriptorList.add(kurtosis.evaluate(distances));

		if (withCovarianceEigenvalues) {
			double[][] centeredPoints = new double[reprPoints.length][3];
			for (int iPoint = 0; iPoint < reprPoints.length; iPoint++) {
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

			geomDescriptorList.add(Math.sqrt(eigenvalues[0]));
			geomDescriptorList.add(Math.sqrt(eigenvalues[1]));
			geomDescriptorList.add(Math.sqrt(eigenvalues[2]));
		}

		geometryDescriptor = geomDescriptorList.stream().mapToDouble(d -> d).toArray();
	}

	private void calcMomentDescriptor() {
		List<Double> momentDescriptorList = new ArrayList<>();

		for (int i = 0; i < config.normOrders.length; i++) {
			List<Double> normedMoments = momentInvariantsRaw.get(i);
			int[] selIndices = config.indicesZernike.get(i);
			for (int selIndex : selIndices) {
				momentDescriptorList.add(normedMoments.get(selIndex));
			}
		}
		momentDescriptor = momentDescriptorList.stream().mapToDouble(d -> d).toArray();
	}

	private void calcMomentInvariantsRaw(InvariantNorm normalization) {
		momentInvariantsRaw = new ArrayList<>();

		for (int normOrder : config.normOrders) {
			momentInvariantsRaw.add(new ArrayList<>(normalization.getInvariants(normOrder)));
		}
	}

	/**
	 * Get the BioZernike geometry descriptors, based on the passed {@link DescriptorConfig} configurations:
	 * an array of size 2 + length of{@link #PERCENTILES_FOR_GEOM} + 3 [+ 3, if {@link DescriptorConfig#withCovEigenValsInGeom} is true],
	 * with values:
	 * <li>radius</li>
	 * <li>residues nominal weight</li>
	 * <li>percentiles ({@link #PERCENTILES_FOR_GEOM}) of distance to centroid distribution</li>
	 * <li>standard deviation of distance to centroid distribution</li>
	 * <li>skewness of distance to centroid distribution</li>
	 * <li>kurtosis of distance to centroid distribution</li>
	 * <li>optionally (if {@link DescriptorConfig#withCovEigenValsInGeom} is true) 3 eigenvalues of distances to centroid covariance matrix</li>
	 * @return the array with geometry descriptors
	 */
	public double[] getGeometryDescriptor() {
		return geometryDescriptor;
	}

	/**
	 * Get the BioZernike moment invariant descriptors, based on the passed {@link DescriptorConfig} configurations
	 * @return the array of moment descriptors
	 */
	public double[] getMomentDescriptor() {
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


