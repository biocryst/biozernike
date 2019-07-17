package org.rcsb.biozernike;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.*;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;

import org.biojava.nbio.structure.*;

import org.biojava.nbio.structure.secstruc.SecStrucState;
import org.rcsb.biozernike.volume.Volume;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Point3d;

import java.util.*;
import java.util.stream.Stream;

public class Descriptor {
	private static final Logger logger = LoggerFactory.getLogger(Descriptor.class);
	final private int MAX_ORDER_INVARIANTS = 20;
	final private int MAX_ORDER_ALIGNMENT = 6;
	final private int[] normOrders = new int[]{2, 3, 4, 5};

	private List<Double> geometryInvariants;
	private List<Double> momentInvariants;
	// tmp
	public double[] distances;

	private Map<Map.Entry<Integer, Integer>, List<MomentTransform>> transformsMap;
	private double[] volumeCenter;

	private Map<Map.Entry<Integer, Integer>, List<MomentTransform>> transformsMapAsym = null;
	private double[] volumeCenterAsym = null;

	public Descriptor(Chain chain) throws Exception {
		// assume ss is calculated
		Atom[] reprAtoms = StructureTools.getRepresentativeAtomArray(chain);
		if (reprAtoms.length == 0) {
			throw new StructureException("Chain " + chain.getId() + " does not have any representative atoms.");
		}
		init(reprAtoms);
	}

	public Descriptor(Structure structure) throws Exception {
		// extract atoms
		Atom[] reprAtoms = StructureTools.getRepresentativeAtomArray(structure);
		if (reprAtoms.length == 0) {
			throw new StructureException("Structure " + structure.getName() + " does not have any representative atoms.");
		}
		init(reprAtoms);
	}

	public Descriptor(Structure structure, String[] asymChainIds) throws Exception {

		Atom[] reprAtoms = StructureTools.getRepresentativeAtomArray(structure);
		if (reprAtoms.length == 0) {
			throw new StructureException("Structure " + structure.getName() + " does not have any representative atoms.");
		}

		if (asymChainIds.length == 0) {
			init(reprAtoms);
			return;
		}

		Point3d[] reprPoints = Calc.atomsToPoints(reprAtoms);
		String[] resNames = Arrays.stream(reprAtoms).map(a -> a.getGroup().getPDBName()).toArray(String[]::new);

		Volume symVolume = new Volume();
		symVolume.create(reprPoints, resNames);
		Volume asymVolume = new Volume(symVolume);
		Volume volumeAU = new Volume();

		List<Atom[]> reprAtomsAUlist = new ArrayList<>();
		for (String asymChainId : asymChainIds) {
			Chain auChain = structure.getPolyChain(asymChainId);
			reprAtomsAUlist.add(StructureTools.getRepresentativeAtomArray(auChain));
		}

		Atom[] reprAtomsAU = reprAtomsAUlist.stream().flatMap(Stream::of).toArray(Atom[]::new);
		String[] resNamesAU = Arrays.stream(reprAtomsAU).map(a -> a.getGroup().getPDBName()).toArray(String[]::new);
		Point3d[] reprPointsAU = Calc.atomsToPoints(reprAtomsAU);

		double[] resCoefsAU = new double[reprAtomsAU.length];
		Arrays.fill(resCoefsAU, 1.0);

		volumeAU.create(reprPointsAU, resNamesAU, resCoefsAU, symVolume.getBoundingBox(), symVolume.getGridWidth());
//		volumeAU.create(reprPointsAU, resNamesAU, getSSCoefs(reprAtomsAU), symVolume.getBoundingBox(), symVolume.getGridWidth());
		asymVolume.add(volumeAU);

		// volume may be cut after this
		// original
		InvariantNorm symNormalization = new InvariantNorm(symVolume, MAX_ORDER_INVARIANTS);

		// TODO: align coord centers after rotating around volume centers?
		// TODO: why translation wiggles?
		// original
		calcGeomInvariants(symVolume, reprPoints); // uses volume radius and weight (has to be non-modified)
		// original
		calcMomentInvariants(symNormalization);

		InvariantNorm asymNormalization = new InvariantNorm(asymVolume, MAX_ORDER_ALIGNMENT);
		// modified
		// TODO: calc a few alignments for the sym volume as well

		transformsMap = calcAlignments(symNormalization);
		volumeCenter = symVolume.getCenterReal();

		// modified, stored for alignment
		transformsMapAsym = calcAlignments(asymNormalization);
		volumeCenterAsym = asymVolume.getCenterReal();
	}


	private void init(Atom[] reprAtoms) {
		Point3d[] reprPoints = Calc.atomsToPoints(reprAtoms);
		String[] resNames = Arrays.stream(reprAtoms).map(a -> a.getGroup().getPDBName()).toArray(String[]::new);

		Volume volume = new Volume();
//		volume.create(reprPoints,resNames,getSSCoefs(reprAtoms));
		volume.create(reprPoints, resNames);

		InvariantNorm normalization = new InvariantNorm(volume, MAX_ORDER_INVARIANTS);
		volumeCenter = volume.getCenterReal();
		calcGeomInvariants(volume, reprPoints);
		calcMomentInvariants(normalization);
		transformsMap = calcAlignments(normalization);
	}

	private double[] getSSCoefs(Atom[] reprAtoms) {
		double[] ssCoef = Arrays.stream(reprAtoms).
				mapToDouble(a -> {
					Character type = ' ';
					try {
						type = ((SecStrucState) a.getGroup().getProperty(Group.SEC_STRUC)).getType().type;
					} catch (Exception ignored) {
					}

					return (type == 'H' || type == 'E') ? 1 : 0.5;
				}).
				toArray();
		return ssCoef;
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

		distances = new double[reprPoints.length];
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
		momentInvariants = new ArrayList<>();

		List<Double> mergedInvariants = new ArrayList<>(normalization.getFingerprint());

		for (int normOrder : normOrders) {
			mergedInvariants.addAll(normalization.getInvariants(normOrder));
		}

		for (int i = 0; i < SearchCoefficients.searchIndsZ.length; i++) {
			momentInvariants.add(mergedInvariants.get(SearchCoefficients.searchIndsZ[i]));
		}
	}

	private Map<Map.Entry<Integer, Integer>, List<MomentTransform>> calcAlignments(InvariantNorm normalization) {
		normalization.setMaxOrder(MAX_ORDER_ALIGNMENT);

		List<Map.Entry<Integer, Integer>> normKeys = new ArrayList<>();
		normKeys.add(new AbstractMap.SimpleImmutableEntry<>(2, 2));
		normKeys.add(new AbstractMap.SimpleImmutableEntry<>(3, 5));
		normKeys.add(new AbstractMap.SimpleImmutableEntry<>(4, 4));
		normKeys.add(new AbstractMap.SimpleImmutableEntry<>(6, 6));

		normalization.setNormalisations(normKeys);
		return normalization.getTransformsMap();
	}

	public static double[] getSearchCoefsZ() {
		return SearchCoefficients.searchCoefsZ;
	}

	public static double getSearchThresholdZ() {
		return SearchCoefficients.thresholdZ;
	}

	public static double[] getSearchCoefsGeom() {
		return SearchCoefficients.searchCoefsGeom;
	}

	public static double getSearchThresholdGeom() {
		return SearchCoefficients.thresholdGeom;
	}


	public List<Double> getGeomInvariants() {
		return geometryInvariants;
	}

	public List<Double> getMomentInvariants() {
		return momentInvariants;
	}

	public Map<Map.Entry<Integer, Integer>, List<MomentTransform>> getTransformsMap() {
		return transformsMap;
	}

	public Map<Map.Entry<Integer, Integer>, List<MomentTransform>> getTransformsMapAsym() {
		return transformsMapAsym;
	}

	public double[] getVolumeCenter() {
		return volumeCenter;
	}

	public double[] getVolumeCenterAsym() {
		return volumeCenterAsym;
	}
}


