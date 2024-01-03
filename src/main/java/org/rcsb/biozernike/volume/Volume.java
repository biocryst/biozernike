package org.rcsb.biozernike.volume;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.rcsb.biozernike.zernike.GeometricMoments;
import javax.media.j3d.BoundingBox;
import javax.media.j3d.Bounds;
import javax.vecmath.Point3d;
import java.util.Arrays;

public class Volume {

	public static final String DEFAULT_RESIDUE_NAME = "ALA";

	/**
	 * The default radius multiplier as documented in Guzenko et al, PLoS CB 2020. See "BioZernike descriptors" section
	 */
	public static final double DEFAULT_RADIUS_MULTIPLIER = 1.8;

	public static final int DEFAULT_MAX_VOLUME_SIZE = 200;
	public static final int DEFAULT_MIN_VOLUME_SIZE = 50;
	public static final int DEFAULT_MIN_INTERFACE_VOXELS = 200;
	public static final double DEFAULT_ATOM_DISTANCE_PADDING = 5.5;

	/**
	 * if the target volume is larger than maxVolumeSize^3, it will be downscaled
 	 */
	private int maxVolumeSize;

	/**
	 * if the target volume is smaller than minVolumeSize^3, it will be upscaled
 	 */
	private int minVolumeSize;

	private int minInterfaceVoxels;
	private double atomDistancePadding;

	/**
	 * A scaling factor for the radius
	 */
	private double radiusVarMult;

	private double volumeMass = 0;
	private double originalVolumeMass = 0;

	private double[] voxelArray = null;
	private int[] dimensions = {0, 0, 0};
	private double[] center = {0, 0, 0};
	private double[] corner = {0, 0, 0};
	private double radiusVar = 0;
	private double radiusMax = 0;

	private double gridWidth = 0;
	private double residuesNominalWeight = 0;

	private BoundingBox bb = new BoundingBox((Bounds) null);

	public Volume() {
		this.radiusVarMult = DEFAULT_RADIUS_MULTIPLIER;
		this.maxVolumeSize = DEFAULT_MAX_VOLUME_SIZE;
		this.minVolumeSize = DEFAULT_MIN_VOLUME_SIZE;
		this.minInterfaceVoxels = DEFAULT_MIN_INTERFACE_VOXELS;
		this.atomDistancePadding = DEFAULT_ATOM_DISTANCE_PADDING;
	}

	public Volume(Volume other) {
		voxelArray = new double[other.voxelArray.length];
		System.arraycopy( other.voxelArray, 0, this.voxelArray, 0, other.voxelArray.length);
		System.arraycopy( other.dimensions, 0, this.dimensions, 0, other.dimensions.length);
		System.arraycopy( other.center, 0, this.center, 0, other.center.length);
		System.arraycopy( other.corner, 0, this.corner, 0, other.corner.length);

		this.originalVolumeMass = other.originalVolumeMass;
		this.radiusVar = other.radiusVar;
		this.radiusMax = other.radiusMax;
		this.volumeMass = other.volumeMass;
		this.gridWidth = other.gridWidth;
		this.residuesNominalWeight = other.residuesNominalWeight;
		this.bb = new BoundingBox(other.bb);

		this.maxVolumeSize = other.maxVolumeSize;
		this.minVolumeSize = other.minVolumeSize;
		this.minInterfaceVoxels = other.minInterfaceVoxels;
		this.atomDistancePadding = other.atomDistancePadding;
		this.radiusVarMult = other.radiusVarMult;
	}

	private void reset() {
		voxelArray = null;
		originalVolumeMass = 0;
		dimensions[0] = dimensions[1] = dimensions[2] = 0;
		center[0] = center[1] = center[2] = 0;
		corner[0] = corner[1] = corner[2] = 0;
		radiusVar = 0;
		radiusMax = 0;
		volumeMass = 0;
		gridWidth = 0;
		residuesNominalWeight = 0;
		bb = new BoundingBox((Bounds) null);
	}

	public void add(Volume other) {
		if (this.dimensions[0] != other.dimensions[0] ||
				this.dimensions[1] != other.dimensions[1] ||
				this.dimensions[2] != other.dimensions[2]) {
			throw new IllegalArgumentException("Volume dimensions must match");
		}

		for (int flatInd = 0; flatInd < this.voxelArray.length; flatInd++) {
			this.voxelArray[flatInd] += other.voxelArray[flatInd];
		}
		updateCenter();
	}

	/**
	 * Create a volume given the dimensions, the array of voxel values and the grid width.
	 * Note that the center will not be calculated automatically. The caller must call {@link #updateCenter()} subsequently
 	 * @param dimensions array of length 3 with the dimensions
	 * @param voxelArray flattened array with the density values
	 * @param gridWidth the grid width in Angstroms
	 */
	public void createFromData(int[] dimensions, double[] voxelArray, double gridWidth) {
		this.dimensions = dimensions;
		this.voxelArray = voxelArray;
		this.gridWidth = gridWidth;
	}

	public void create(Point3d[] reprCoords) {
		String[] resNames = new String[reprCoords.length];
		Arrays.fill(resNames, DEFAULT_RESIDUE_NAME);

		double[] resCoefs = new double[reprCoords.length];
		Arrays.fill(resCoefs,1.0);

		BoundingBox bb = new BoundingBox((Bounds) null);
		bb.combine(reprCoords);

		create(reprCoords,resNames, resCoefs, bb, getAutoGridWidth(bb));
	}

	public void create(Point3d[] reprCoords, double gridWidth) {
		String[] resNames = new String[reprCoords.length];
		Arrays.fill(resNames, DEFAULT_RESIDUE_NAME);

		double[] resCoefs = new double[reprCoords.length];
		Arrays.fill(resCoefs,1.0);

		BoundingBox bb = new BoundingBox((Bounds) null);
		bb.combine(reprCoords);

		create(reprCoords, resNames, resCoefs, bb, gridWidth);
	}

	public void create(Point3d[] reprCoords, String[] resNames) {
		double[] resCoefs = new double[reprCoords.length];
		Arrays.fill(resCoefs,1.0);

		BoundingBox bb = new BoundingBox((Bounds) null);
		bb.combine(reprCoords);

		create(reprCoords,resNames, resCoefs, bb, getAutoGridWidth(bb));
	}

	public void create(Point3d[] reprCoords, String[] resNames, double gridWidth) {
		double[] resCoefs = new double[reprCoords.length];
		Arrays.fill(resCoefs,1.0);

		BoundingBox bb = new BoundingBox((Bounds) null);
		bb.combine(reprCoords);

		create(reprCoords,resNames, resCoefs, bb, gridWidth);
	}

	public void create(Point3d[] reprCoords, String[] resNames, double[] resCoefs) {
		BoundingBox bb = new BoundingBox((Bounds) null);
		bb.combine(reprCoords);
		create(reprCoords,resNames, resCoefs, bb, getAutoGridWidth(bb));
	}

	public void create(Point3d[] reprCoords, String[] resNames, double[] resCoefs, BoundingBox bb) {
		create(reprCoords,resNames,resCoefs, bb, getAutoGridWidth(bb));
	}

	public void create(Point3d[] reprCoords, String[] resNames, double[] resCoefs, BoundingBox bb, double gridWidth) {
		reset();
		this.bb = bb;
		this.gridWidth = gridWidth;
		if (gridWidth == 0 || !ResidueVolumeCache.isValidGridWidth(gridWidth)) {
			this.gridWidth = getAutoGridWidth(bb);
		}
		this.voxelArray = fillVolume(reprCoords, resNames, resCoefs, bb);
		updateCenter();
	}

	public void createFromInterface(Point3d[] reprCoords1, Point3d[] reprCoords2) {
		String[] resNames1 = new String[reprCoords1.length];
		String[] resNames2 = new String[reprCoords2.length];
		Arrays.fill(resNames1, DEFAULT_RESIDUE_NAME);
		Arrays.fill(resNames2, DEFAULT_RESIDUE_NAME);

		createFromInterface(reprCoords1, reprCoords2, resNames1, resNames2, 0);
	}

	public void createFromInterface(Point3d[] reprCoords1, Point3d[] reprCoords2, String[] resNames1, String[] resNames2) {
		createFromInterface(reprCoords1, reprCoords2, resNames1, resNames2, 0);
	}

	public void createFromInterface(Point3d[] reprCoords1, Point3d[] reprCoords2, String[] resNames1, String[] resNames2, double gridWidth) {

		reset();

		BoundingBox bb2 = new BoundingBox((Bounds) null);

		bb.combine(reprCoords1);
		bb2.combine(reprCoords2);

		Point3d pLower = new Point3d();
		Point3d pUpper = new Point3d();
		bb.getUpper(pUpper);
		bb2.getLower(pLower);

		Point3d pPadding = new Point3d(atomDistancePadding, atomDistancePadding, atomDistancePadding);
		pUpper.add(pPadding);
		pLower.sub(pPadding);

		BoundingBox bb1Ext = new BoundingBox(pLower, pUpper);

		if (!bb1Ext.intersect(bb2)) {
			return;
		}

		bb.combine(bb2);

		this.gridWidth = gridWidth;
		if (gridWidth == 0 || !ResidueVolumeCache.isValidGridWidth(gridWidth)) {
			this.gridWidth = getAutoGridWidth(bb);
		}

		double[] resCoefs = new double[Math.max(reprCoords1.length,reprCoords2.length)];
		Arrays.fill(resCoefs,1.0);

		double[] chain1Array = fillVolume(reprCoords1, resNames1, resCoefs, bb);
		double[] chain2Array = fillVolume(reprCoords2, resNames2, resCoefs, bb);

		this.voxelArray = new double[chain1Array.length];

		int nFilled = 0;
		for (int flatInd = 0; flatInd < chain1Array.length; flatInd++) {
			if (chain1Array[flatInd] > 0 && chain2Array[flatInd] > 0) {
				double sumVol = chain1Array[flatInd] + chain2Array[flatInd];
				this.voxelArray[flatInd] = sumVol;
				nFilled++;
			}
		}
		if (nFilled < minInterfaceVoxels) {
			reset();
		} else {
			updateCenter();
		}
	}

	private double getAutoGridWidth(BoundingBox bb) {
		double gridWidth = 1;
		Point3d pLower = new Point3d();
		Point3d pUpper = new Point3d();
		bb.getUpper(pUpper);
		bb.getLower(pLower);
		pUpper.sub(pLower);

		double bbCubeSize = Math.cbrt(pUpper.x * pUpper.y * pUpper.z);
		while ((bbCubeSize > maxVolumeSize) && (gridWidth < ResidueVolumeCache.MAX_GRID_WIDTH)) {
			gridWidth *= 2;
			bbCubeSize /= 2;
		}
		while ((bbCubeSize < minVolumeSize) && (gridWidth > ResidueVolumeCache.MIN_GRID_WIDTH)) {
			gridWidth /= 2;
			bbCubeSize *= 2;
		}

		return gridWidth;
	}

	private double[] fillVolume(Point3d[] reprCoords, String[] resNames, double[] resCoefs, BoundingBox bb) {

		// make some space for the blobs
		Point3d pMin = new Point3d();
		Point3d pDims = new Point3d();
		bb.getLower(pMin);
		bb.getUpper(pDims);
		pDims.sub(pMin);

		int maxBoxSize = ResidueVolumeCache.maxBoxSize.get(gridWidth);

		this.dimensions = new int[]{
				(int) Math.ceil((pDims.x) / gridWidth) + maxBoxSize,
				(int) Math.ceil((pDims.y) / gridWidth) + maxBoxSize,
				(int) Math.ceil((pDims.z) / gridWidth) + maxBoxSize};

		corner = new double[]{
				pMin.x - maxBoxSize * gridWidth / 2.0,
				pMin.y - maxBoxSize * gridWidth / 2.0,
				pMin.z - maxBoxSize * gridWidth / 2.0
		};

		int flatVolumeSize = this.dimensions[0] * this.dimensions[1] * this.dimensions[2];
		double[] structureVolume = new double[flatVolumeSize];
		int nAtoms = reprCoords.length;
		int dims01 = dimensions[0] * dimensions[1];

		for (int indAtom = 0; indAtom < nAtoms; indAtom++) {
			Point3d atom = reprCoords[indAtom];

			String resName = resNames[indAtom];
			double[] coords = {atom.x, atom.y, atom.z};
			double weightMultiplier = (resCoefs==null)?1:resCoefs[indAtom];

			residuesNominalWeight += VolumeConstants.getWeight(resName)*weightMultiplier;
			ResidueVolume resBox = ResidueVolumeCache.get(gridWidth, resName);
			int boxDim = resBox.size;
			int boxDim2 = boxDim * boxDim;
			int[] coordsBoxCorner = new int[3];

			for (int i = 0; i < 3; i++) {
				coordsBoxCorner[i] = (int) Math.round((coords[i] - corner[i]) / gridWidth - resBox.size/2.0);
			}

			int zVolumeInd = (coordsBoxCorner[2] * dimensions[1] + coordsBoxCorner[1]) * dimensions[0] + coordsBoxCorner[0];
			int zBoxInd = 0;

			for (int zBox = 0; zBox < boxDim; zBox++) {
				int yVolumeInd = zVolumeInd;
				int yBoxInd = zBoxInd;

				for (int yBox = 0; yBox < boxDim; yBox++) {
					int flatVolumeInd = yVolumeInd;
					int flatBoxInd = yBoxInd;

					for (int xBox = 0; xBox < boxDim; xBox++) {
						structureVolume[flatVolumeInd] += resBox.volume[flatBoxInd]*weightMultiplier;
						flatVolumeInd++;
						flatBoxInd++;
					}
					yVolumeInd += dimensions[0];
					yBoxInd += boxDim;
				}
				zVolumeInd += dims01;
				zBoxInd += boxDim2;
			}
		}

		return structureVolume;
	}

	public void normalize() {
		if (isNormalized()) {
			return;
		}
		originalVolumeMass = volumeMass;
		computeRadius();
		// do not recompute the radius after trimming!
		// if you do, extended structures may be trimmed to nothing.
		while (trim(radiusVar) > 0) {
			// need to update the center to accurately compute the geometric moments
			updateCenter();
		}
	}

	public void updateCenter() {
		center[0] = center[1] = center[2] = 0;
		GeometricMoments gm = new GeometricMoments(this, 1, 1);

		volumeMass = gm.getMoment(0, 0, 0);
		center[0] = gm.getMoment(1, 0, 0) / volumeMass;
		center[1] = gm.getMoment(0, 1, 0) / volumeMass;
		center[2] = gm.getMoment(0, 0, 1) / volumeMass;
	}

	private void computeRadius() {

		double weightedRad = 0; // Rg
		double maxRad = 0;

		int zInd = 0;
		int dims01 = dimensions[0] * dimensions[1];

		for (int z = 0; z < dimensions[2]; ++z) {
			double mz2 = z - center[2];
			mz2 = mz2 * mz2;
			int yInd = zInd;

			for (int y = 0; y < dimensions[1]; ++y) {
				double mzy2 = y - center[1];
				mzy2 = mzy2 * mzy2 + mz2;
				int flatInd = yInd;

				for (int x = 0; x < dimensions[0]; ++x) {
					double voxelWeight = voxelArray[flatInd];

					if (voxelWeight > 0) {
						double mx = x - center[0];
						double sqRad = mzy2 + mx * mx;

						if (sqRad > maxRad) {
							maxRad = sqRad;
						}
						weightedRad += sqRad * voxelWeight / volumeMass;
					}
					flatInd++;
				}
				yInd += dimensions[0];
			}
			zInd += dims01;
		}

		radiusVar = Math.sqrt(weightedRad) * radiusVarMult;
		radiusMax = Math.sqrt(maxRad);
	}

	private int trim(double radius) {
		double sqrRadius = radius * radius;
		int nZeroed = 0;
		//(z * dimensions[1] + y) * dimensions[0] + x
		int zInd = 0;
		int dims01 = dimensions[0] * dimensions[1];

		for (int z = 0; z < dimensions[2]; ++z) {

			double mz2 = z - center[2];
			mz2 = mz2 * mz2;
			int yzInd = zInd;

			for (int y = 0; y < dimensions[1]; ++y) {
				double mzy2 = y - center[1];
				mzy2 = mzy2 * mzy2 + mz2;
				int flatInd = yzInd;

				for (int x = 0; x < dimensions[0]; ++x) {
					double voxelWeight = voxelArray[flatInd];
//					assert voxelWeight == getValue(x,y,z);
					if (voxelWeight > 0) {
						double mx = x - center[0];
						double sqRadTmp = mzy2 + mx * mx;
						if (sqRadTmp > sqrRadius) {
							volumeMass -= voxelWeight;
							voxelArray[flatInd] = 0;
							nZeroed++;
						}
					}
					flatInd++;
				}
				yzInd += dimensions[0];
			}
			zInd += dims01;
		}
//		System.out.println("Zeroed: "+nZeroed+", Sum of all new: "+volumeMass);

		return nZeroed;
	}

	/**
	 * Apply a lower threshold, setting to 0 voxels below the threshold. And normalize given a multiplier.
	 * <p>
	 * Note that the center will not be calculated automatically. The caller must call {@link #updateCenter()} subsequently
	 * @param contourThreshold a lower threshold of density values to consider. Values below this threshold will be set to 0,
	 *                         if below 0, then parameter is ignore and instead 3x std-deviation used
	 * @param multiplier a multiplier value to normalise the density values
	 */
	public void applyContourAndNormalize(double contourThreshold, double multiplier) {
		double sumval = 0;
		int nVoxels = 0;

		for (int i = 0; i<voxelArray.length; i++) {
			if(voxelArray[i] <= contourThreshold) {
				voxelArray[i] = 0;
				continue;
			}
			// count voxels with a value > threshold
			nVoxels++;
			sumval += voxelArray[i];
		}

		double normCoef = nVoxels*multiplier/sumval;
		for (int i = 0; i<voxelArray.length; i++) {
			voxelArray[i] *= normCoef;
		}
	}

	/**
	 * Apply a lower threshold expressed in units of standard deviations. And normalize given a multiplier.
	 * <p>
	 * Note that the center will not be calculated automatically. The caller must call {@link #updateCenter()} subsequently
	 * @param stdDevMultiplier value to multiply the standard deviation that will be the contour threshold
	 * @param multiplier a multiplier value to normalise the density values
	 */
	public void applyContourStdAndNormalize(double stdDevMultiplier, double multiplier) {
		StandardDeviation standardDeviation = new StandardDeviation();
		double stdDev = standardDeviation.evaluate(voxelArray);
		applyContourAndNormalize(stdDev * stdDevMultiplier, multiplier);
	}

	/**
	 * Convert all density values to positive (taking absolute value).
	 * Note that EMDB maps usually have negative values. This will convert those to positive, which is better suited for
	 * BioZernike descriptors calculation, especially when comparing to volumes that are calculated from
	 * coordinates (where all density is positive).
	 * <p>
	 * Note that the center will not be calculated automatically. The caller must call {@link #updateCenter()} subsequently
	 */
	public void positivize() {
		for (int i = 0; i<voxelArray.length; i++) {
			voxelArray[i] = Math.abs(voxelArray[i]);
		}
	}

	public DescriptiveStatistics getDescriptiveStatistics() {
		return new DescriptiveStatistics(voxelArray);
	}

	public double getValue(int x, int y, int z) {
		return voxelArray[(z * dimensions[1] + y) * dimensions[0] + x];
	}

	public int getInd(int x, int y, int z) {
		return (z * dimensions[1] + y) * dimensions[0] + x;
	}

	public boolean isNormalized() {
		return radiusVar > 0 && radiusMax > 0;
	}

	public boolean isEmpty() {
		return voxelArray == null;
	}

	public double getRadiusVarMult() {
		return radiusVarMult;
	}

	public void setRadiusVarMult(double radiusVarMult) {
		this.radiusVarMult = radiusVarMult;
	}

	/**
	 * Get the radius of gyration in voxel units
	 * @return the radius of gyration in voxel units
	 */
	public double getRadiusVarVolume() {
		return radiusVar;
	}

	/**
	 * Get the maximum radius in voxel units
	 * @return the maximum radius in voxel units
	 */
	public double getRadiusMaxVolume() {
		return radiusMax;
	}

	public double getTrimmedVolume() {
		return volumeMass / originalVolumeMass;
	}

	/**
	 * Get the radius of gyration in Angstroms
	 * @return the radius of gyration in Angstroms
	 */
	public double getRadiusVarReal() {
		return radiusVar * gridWidth;
	}

	/**
	 * Get the maximum radius in Angstroms
	 * @return the maximum radius in Angstroms
	 */
	public double getRadiusMaxReal() {
		return radiusMax * gridWidth;
	}

	public double getResiduesNominalWeight() {
		return residuesNominalWeight;
	}

	public double[] getVoxelArray() {
		return voxelArray;
	}

	public int[] getDimensions() {
		return dimensions;
	}

	public double[] getCenterVolume() {
		return center;
	}

	/**
	 * Get the volume's center coordinate in Angstrom units
	 * @return the center coordinate
	 */
	public double[] getCenterReal() {
		return new double[]{
				center[0] * gridWidth + corner[0],
				center[1] * gridWidth + corner[1],
				center[2] * gridWidth + corner[2]
		};
	}

	public void setCenterReal(double[] centerReal) {
		center = new double[]{
				(centerReal[0] - corner[0]) / gridWidth ,
				(centerReal[1] - corner[1]) / gridWidth ,
				(centerReal[2] - corner[2]) / gridWidth
		};
		computeRadius();
	}

	/**
	 * Get the grid width in Angstroms
	 * @return the grid width in Angstroms
	 */
	public double getGridWidth() {
		return gridWidth;
	}

	public BoundingBox getBoundingBox() {
		return bb;
	}

	public double getVolumeMass() {
		return volumeMass;
	}

	public double getOriginalVolumeMass() {
		return originalVolumeMass;
	}

	public int getMaxVolumeSize() {
		return maxVolumeSize;
	}

	public void setMaxVolumeSize(int maxVolumeSize) {
		this.maxVolumeSize = maxVolumeSize;
	}

	public int getMinVolumeSize() {
		return minVolumeSize;
	}

	public void setMinVolumeSize(int minVolumeSize) {
		this.minVolumeSize = minVolumeSize;
	}

	public int getMinInterfaceVoxels() {
		return minInterfaceVoxels;
	}

	public void setMinInterfaceVoxels(int minInterfaceVoxels) {
		this.minInterfaceVoxels = minInterfaceVoxels;
	}

	public double getAtomDistancePadding() {
		return atomDistancePadding;
	}

	public void setAtomDistancePadding(double atomDistancePadding) {
		this.atomDistancePadding = atomDistancePadding;
	}

}
