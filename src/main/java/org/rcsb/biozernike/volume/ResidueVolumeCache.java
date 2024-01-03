package org.rcsb.biozernike.volume;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.HashMap;
import java.util.Map;

/**
 * Pre-computed Gaussian blobs for all (supported) residues at several grid widths: 1,2,4,8 A.
 * Essential functionality follows gmconvert.
 */
public class ResidueVolumeCache {

	private static final Logger logger = LoggerFactory.getLogger(ResidueVolumeCache.class);

	private static final  Map<Double, Map<String, ResidueVolume>> residueBox = new HashMap<>();
	private static final double sdCutoff = 3.0;
	private static final double three_over_2pi_32 = 0.3299226101861591;
	private static final double densityMultiplier = 100;

	public static final double MAX_GRID_WIDTH = 16;
	public static final double MIN_GRID_WIDTH = 0.25;
	public static final double[] GRID_WIDTHS = {1, 0.25, 0.5, 2, 4, 8, 16};

	public static  Map<Double, Integer> maxBoxSize = new HashMap<>();

	static {
		String[] resNames = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LEU", "LYS", "MET", "MSE",
				"PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
//				"A", "T", "C", "G", "U", "I", "DA", "DT", "DG", "DC", "DU", "DI",
				"H", "C", "N", "O", "F", "P", "CL", "CU"
		};

		Map<String, Double> refWeight = new HashMap<>(resNames.length);

		for (double gridWidth : GRID_WIDTHS) {
			maxBoxSize.put(gridWidth, 0);

			Map<String, ResidueVolume> residueBoxGrid = new HashMap<>(resNames.length);

			for (String resName : resNames) {
				double resRadius = VolumeConstants.residueRadius.get(resName);
				double resWeight = VolumeConstants.residueWeight.get(resName);
				double var = resRadius * resRadius / 5.0;
				double sigma = Math.sqrt(var);

				int boxSize = (int) (Math.ceil(2 * sdCutoff * sigma / gridWidth));

				if (boxSize % 2 == 0) {
					boxSize += 1;
				}
				if (boxSize > maxBoxSize.get(gridWidth)) {
					maxBoxSize.put(gridWidth, boxSize);
				}

				int flatBoxSize = boxSize * boxSize * boxSize;
				double[] volume = new double[flatBoxSize];

				int center = boxSize / 2;
				double sqrRadius = center * center;

				double boxSum = 0;
				for (int x = 0; x < boxSize; x++) {
					double mx = x - center;
					for (int y = 0; y < boxSize; y++) {
						double my = y - center;
						for (int z = 0; z < boxSize; z++) {
							double mz = z - center;
							if (mx * mx + my * my + mz * mz <= sqrRadius) {
								double[] dist = {mx * gridWidth, my * gridWidth, mz * gridWidth};
								double xSx = (dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]) / var;
								double gausVal = resWeight * three_over_2pi_32 / (var * sigma) * Math.exp(-0.5 * xSx);
								volume[(z * boxSize + y) * boxSize + x] = gausVal * densityMultiplier;
								boxSum += gausVal;
							}
						}
					}
				}

				if (gridWidth == 1) {
					refWeight.put(resName, boxSum);
				} else {
					double sumCoef = refWeight.get(resName) / (boxSum * gridWidth * gridWidth * gridWidth);
					for (int i = 0; i < flatBoxSize; i++) {
						volume[i] *= sumCoef;
					}
				}

				residueBoxGrid.put(resName, new ResidueVolume(volume, boxSize));
			}
			residueBox.put(gridWidth, residueBoxGrid);
		}
	}

	static boolean isValidGridWidth(double gridWidth) {
		return residueBox.containsKey(gridWidth);
	}

	public static ResidueVolume get(double gridWidth, String resName) {
		if (!residueBox.get(gridWidth).containsKey(resName)) {
			resName = "ASP";
		}
		return residueBox.get(gridWidth).get(resName);
	}
}
