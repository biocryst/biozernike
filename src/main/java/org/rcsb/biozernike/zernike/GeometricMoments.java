package org.rcsb.biozernike.zernike;

import org.rcsb.biozernike.volume.Volume;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Computation of geometric moments, needed for 3D Zernike moment calculation.
 * See section 4.1 of Novotni and Klein 2003 (https://cg.cs.uni-bonn.de/aigaion2root/attachments/novotni-2003-3d.pdf)
 */
public class GeometricMoments {

	private double[] voxelArray;
	private double[] center;
	private int[] volumeDims;
	private double scale;
	private int maxOrder;
	private List<List<List<Double>>> moments;
	private double[][] samples;

	public GeometricMoments(Volume volume, double radius, int maxOrder) {
		this.voxelArray = volume.getVoxelArray();
		this.center = volume.getCenterVolume();
		this.volumeDims = volume.getDimensions();
		this.scale = 1.0 / radius;
		this.maxOrder = maxOrder;

		moments = new ArrayList<>(maxOrder + 1);

		for (int i = 0; i <= maxOrder; i++) {
			List<List<Double>> momentsi = new ArrayList<>(maxOrder - i + 1);
			for (int j = 0; j <= maxOrder - i; j++) {
				List<Double> momentsj = new ArrayList<>(Collections.nCopies(maxOrder - i - j + 1, 0.0));
				momentsi.add(momentsj);
			}
			moments.add(momentsi);
		}
		computeSamples();
		compute();
	}

	private void computeSamples() {

		double[] min = {-center[0] * scale, -center[1] * scale, -center[2] * scale};
		int maxDim = volumeDims[0];
		if (volumeDims[1] > maxDim)
			maxDim = volumeDims[1];
		if (volumeDims[2] > maxDim)
			maxDim = volumeDims[2];

		samples = new double[3][maxDim + 1];

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j <= volumeDims[i]; j++) {
				samples[i][j] = min[i] + j * scale;
			}
		}
	}


	private void compute() {

		int arrayDim = volumeDims[2];
		int layerDim = volumeDims[1] * volumeDims[2];

		int diffArrayDim = volumeDims[2] + 1;
		int diffLayerDim = (volumeDims[1] + 1) * volumeDims[2];
		int diffGridDim = (volumeDims[0] + 1) * layerDim;

		double[] diffGrid = new double[diffGridDim];
		double[] diffLayer = new double[diffLayerDim];
		double[] diffArray = new double[diffArrayDim];

		double[] layer = new double[layerDim];
		double[] array = new double[arrayDim];
		double moment;

		int iter = 0, diffIter = 0;

		// generate the diff version of the voxel grid in x direction
		for (int x = 0; x < layerDim; x++) {
			computeDiffFunction(iter, voxelArray, diffIter, diffGrid, volumeDims[0]);
			iter += volumeDims[0];
			diffIter += volumeDims[0] + 1;
		}

		for (int i = 0; i <= maxOrder; i++) {
			diffIter = 0;
			for (int p = 0; p < layerDim; p++) {
				// multiply the diff function with the sample values

				int sampleIter = 0;
				layer[p] = multiply(diffIter, diffGrid, sampleIter, samples[0], volumeDims[0] + 1);
				diffIter += volumeDims[0] + 1;

			}

			iter = 0;
			diffIter = 0;

			for (int y = 0; y < arrayDim; y++) {
				computeDiffFunction(iter, layer, diffIter, diffLayer, volumeDims[1]);
				iter += volumeDims[1];
				diffIter += volumeDims[1] + 1;
			}

			for (int j = 0; j < maxOrder + 1 - i; j++) {
				diffIter = 0;
				for (int p = 0; p < arrayDim; p++) {
					int sampleIter = 0;
					array[p] = multiply(diffIter, diffLayer, sampleIter, samples[1], volumeDims[1] + 1);
					diffIter += volumeDims[1] + 1;
				}


				iter = 0;
				diffIter = 0;
				computeDiffFunction(iter, array, diffIter, diffArray, volumeDims[2]);

				for (int k = 0; k < maxOrder + 1 - i - j; ++k) {
					int sampleIter = 0;
					moment = multiply(diffIter, diffArray, sampleIter, samples[2], volumeDims[2] + 1);
					moments.get(i).get(j).set(k, moment / ((i + 1) * (j + 1) * (k + 1)));
				}
			}
		}

	}

	private double multiply(int diffIter, double[] diffGrid, int sampleIter, double[] sampleGrid, int dim) {
		double sum = 0;
		for (int i = 0; i < dim; ++i) {
			diffGrid[diffIter + i] *= sampleGrid[sampleIter + i];
			sum += diffGrid[diffIter + i];
		}
		return sum;
	}

	private void computeDiffFunction(int iter, double[] iterGrid, int diffIter, double[] diffGrid, int dim) {
		diffGrid[diffIter] = -iterGrid[iter];
		for (int i = 1; i < dim; ++i) {
			diffGrid[diffIter + i] = iterGrid[iter + i - 1] - iterGrid[iter + i];
		}
		diffGrid[diffIter + dim] = iterGrid[iter + dim - 1];
	}


	public double getMoment(int i, int j, int k) {
		return moments.get(i).get(j).get(k);
	}

}
