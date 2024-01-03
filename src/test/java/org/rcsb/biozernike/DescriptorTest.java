package org.rcsb.biozernike;

import org.apache.commons.lang.ArrayUtils;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.quaternary.BioAssemblyTools;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.junit.BeforeClass;
import org.junit.Test;
import org.rcsb.biozernike.complex.Complex;
import org.rcsb.biozernike.descriptor.Descriptor;
import org.rcsb.biozernike.descriptor.DescriptorConfig;
import org.rcsb.biozernike.descriptor.DescriptorMode;
import org.rcsb.biozernike.volume.MapFileType;
import org.rcsb.biozernike.volume.Volume;
import org.rcsb.biozernike.volume.VolumeIO;
import org.rcsb.biozernike.zernike.ZernikeMoments;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;


public class DescriptorTest {

	@BeforeClass
	public static void setupBioJava() {
		FileParsingParameters params = new FileParsingParameters();
		params.setParseBioAssembly(true);
		AtomCache cache = new AtomCache();
		cache.setFileParsingParams(params);
		StructureIO.setAtomCache(cache);
	}

	@Test
	public void testSpeed() throws Exception {
		int n_iterations = 500;
		List<Double> tmp = new ArrayList<>(n_iterations);
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

		Structure structure = StructureIO.getStructure("1q2w");
		List<BiologicalAssemblyTransformation> transformations =
				structure.getPDBHeader().getBioAssemblies().get(1).getTransforms();

		Structure bioUnitStructure = builder
				.rebuildQuaternaryStructure(BioAssemblyTools.getReducedStructure(structure), transformations, true, false);

		Atom[] reprAtoms = StructureTools.getRepresentativeAtomArray(bioUnitStructure);
		Point3d[] reprPoints = Calc.atomsToPoints(reprAtoms);
		String[] resNames = Arrays.stream(reprAtoms).map(a -> a.getGroup().getPDBName()).toArray(String[]::new);

		EnumSet<DescriptorMode> mode = EnumSet.allOf(DescriptorMode.class);
		DescriptorConfig config = new DescriptorConfig(DescriptorTest.class.getResourceAsStream("/descriptor.properties"), mode);

		long startTime = System.currentTimeMillis();
		for (int i=0;i<n_iterations;i++) {
			Descriptor ssd = new Descriptor(reprPoints,resNames,config);

			 double[] invariantsSearch = ssd.getMomentDescriptor();
			tmp.add(invariantsSearch[0]);
		}
		long endTime = System.currentTimeMillis();
		double duration = endTime - startTime;
		System.out.println("Performance: "+n_iterations/(duration/1000)+" 4HHBs per second");
	}

	@Test
	public void testMoments() throws Exception {
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

		Structure structure = StructureIO.getStructure("4HHB");
		List<BiologicalAssemblyTransformation> transformations =
				structure.getPDBHeader().getBioAssemblies().get(1).getTransforms();

		Structure bioUnitStructure = builder
				.rebuildQuaternaryStructure(BioAssemblyTools.getReducedStructure(structure), transformations, true, false);

		Atom[] reprAtoms = StructureTools.getRepresentativeAtomArray(bioUnitStructure);
		Point3d[] reprPoints = Calc.atomsToPoints(reprAtoms);
		Volume volume = new Volume();
		volume.create(reprPoints);
		ZernikeMoments zernikeMoments1 = new ZernikeMoments(volume,6);

		List<List<List<Complex>>> originalMoments = zernikeMoments1.getOriginalMoments();
		List<Complex> originalMomentsFlatComplex = ZernikeMoments.flattenMomentsComplex(originalMoments);
		List<Double> originalMomentsFlatDouble = ZernikeMoments.flattenMomentsDouble(originalMoments);

		List<List<List<Complex>>> originalMomentsRestoredComplex = ZernikeMoments.unFlattenMomentsComplex(originalMomentsFlatComplex,zernikeMoments1.getMaxOrder());
		List<List<List<Complex>>> originalMomentsRestoredDouble = ZernikeMoments.unFlattenMomentsDouble(originalMomentsFlatDouble,zernikeMoments1.getMaxOrder());

		assert originalMoments.equals(originalMomentsRestoredComplex);
		assert originalMoments.equals(originalMomentsRestoredDouble);

		ZernikeMoments zernikeMoments2 = new ZernikeMoments(zernikeMoments1.getOriginalMomentsUnscaled(), false);
		ZernikeMoments zernikeMoments3 = new ZernikeMoments(zernikeMoments1.getOriginalMoments(), true);

		double[] originalMomentsArr1 = ArrayUtils.toPrimitive(
				ZernikeMoments.flattenMomentsDouble(
						zernikeMoments1.getOriginalMoments()).
						toArray(new Double[0]));
		double[] originalMomentsArr2 = ArrayUtils.toPrimitive(
				ZernikeMoments.flattenMomentsDouble(
						zernikeMoments2.getOriginalMoments()).
						toArray(new Double[0]));
		double[] originalMomentsArr3 = ArrayUtils.toPrimitive(
				ZernikeMoments.flattenMomentsDouble(
						zernikeMoments3.getOriginalMoments()).
						toArray(new Double[0]));

		assertArrayEquals(originalMomentsArr1,originalMomentsArr2,1e-15);
		assertArrayEquals(originalMomentsArr1,originalMomentsArr3,1e-15);


		double[] originalMomentsUnscaledArr1 = ArrayUtils.toPrimitive(
				ZernikeMoments.flattenMomentsDouble(
						zernikeMoments1.getOriginalMomentsUnscaled()).
						toArray(new Double[0]));
		double[] originalMomentsUnscaledArr2 = ArrayUtils.toPrimitive(
				ZernikeMoments.flattenMomentsDouble(
						zernikeMoments2.getOriginalMomentsUnscaled()).
						toArray(new Double[0]));
		double[] originalMomentsUnscaledArr3 = ArrayUtils.toPrimitive(
				ZernikeMoments.flattenMomentsDouble(
						zernikeMoments3.getOriginalMomentsUnscaled()).
						toArray(new Double[0]));

		assertArrayEquals(originalMomentsUnscaledArr1,originalMomentsUnscaledArr2,1e-15);
		assertArrayEquals(originalMomentsUnscaledArr1,originalMomentsUnscaledArr3,1e-15);
	}

	@Test
	public void testEMAlignment() throws Exception {

		// some hardcoded scaling coefficients for the EM volume (as we do not control the density values)
		InputStream is = DescriptorTest.class.getResourceAsStream("/emd_3186.map");
		Volume volumeEM = VolumeIO.read(is, MapFileType.MRC, true);
		volumeEM.positivize();
		volumeEM.applyContourAndNormalize(0.0176, 757);
		volumeEM.updateCenter();
		volumeEM.setRadiusVarMult(1.64);

		InvariantNorm normalizationEM = new InvariantNorm(volumeEM,6);

		// let's fit a similar structure
		Structure structure = StructureIO.getStructure("1HHS.A");
		Volume volumeStructure = new Volume();
		Atom[] reprAtoms = StructureTools.getRepresentativeAtomArray(structure);
		Point3d[] reprPoints = Calc.atomsToPoints(reprAtoms);
		String[] resNames = Arrays.stream(reprAtoms).map(a -> a.getGroup().getPDBName()).toArray(String[]::new);
		volumeStructure.create(reprPoints, resNames);

		InvariantNorm normalizationStructure = new InvariantNorm(volumeStructure, 6);

		double distance = normalizationStructure.compareInvariants(normalizationEM,2);

		// CN of the order 2,2 are not too different
		assertTrue(distance<10);

		List<InvariantNorm> zc = new ArrayList<>();

		zc.add(normalizationEM);
		zc.add(normalizationStructure);
		AlignmentResult alignmentResult = RotationAlignment.alignMultiple(zc);

		// alignment quality is reasonable
		assertTrue(alignmentResult.getScore()<0.4);

		// uncomment the next lines output the model fitted to the density

//		Matrix4d rot1 = alignmentResult.getTransforms().get(0);
//		Matrix4d rot2 = alignmentResult.getTransforms().get(1);
//
//		rot1.invert();
//		rot1.mul(rot2);
//		Calc.transform(structure,rot1);
//
//		try (PrintWriter out = new PrintWriter("src/test/resources/1HHS.A.fitted.pdb")) {
//			out.println(structure.toPDB());
//		}

	}
}
