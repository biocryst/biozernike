package org.rcsb.biozernike;

import org.apache.commons.lang.ArrayUtils;
import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.quaternary.BioAssemblyTools;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.junit.Test;
import org.rcsb.biozernike.complex.Complex;
import org.rcsb.biozernike.descriptor.Descriptor;
import org.rcsb.biozernike.descriptor.DescriptorConfig;
import org.rcsb.biozernike.descriptor.DescriptorMode;
import org.rcsb.biozernike.volume.Volume;
import org.rcsb.biozernike.zernike.ZernikeMoments;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;
import java.util.stream.Collectors;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

//import static org.rcsb.biozernike.zernike.ZernikeMoments.flattenMomentsDouble;

public class DescriptorTest {
	@Test
	public void testSpeed() throws Exception {
		int n_iterations = 500;
		List<Double> tmp = new ArrayList<>(n_iterations);
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

		Structure structure = StructureIO.getStructure("4HHB");
		List<BiologicalAssemblyTransformation> transformations =
				structure.getPDBHeader().getBioAssemblies().get(1).getTransforms();

		Structure bioUnitStructure = builder
				.rebuildQuaternaryStructure(BioAssemblyTools.getReducedStructure(structure), transformations, true, false);

		Atom[] reprAtoms = StructureTools.getRepresentativeAtomArray(bioUnitStructure);
		Point3d[] reprPoints = Calc.atomsToPoints(reprAtoms);
		String[] resNames = Arrays.stream(reprAtoms).map(a -> a.getGroup().getPDBName()).toArray(String[]::new);

		EnumSet<DescriptorMode> mode = EnumSet.allOf(DescriptorMode.class);
		DescriptorConfig config = new DescriptorConfig("src/test/resources/descriptor.properties", mode);

		long startTime = System.currentTimeMillis();
		for (int i=0;i<n_iterations;i++) {
			Descriptor ssd = new Descriptor(reprPoints,resNames,config);
			List<Double> invariantsSearch = ssd.getMomentDescriptor();
			tmp.add(invariantsSearch.get(0));
		}
		long endTime = System.currentTimeMillis();
		double duration = endTime - startTime;
		System.out.println("Performance: "+n_iterations/(duration/1000)+" 4HHBs per second");
	}

	@Test
	public void testMoments() throws Exception {
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

		AtomCache atomCache = new AtomCache();
		atomCache.setUseMmCif(true);
		atomCache.setUseMmtf(false);
		FileParsingParameters params = new FileParsingParameters();
		params.setParseBioAssembly(true);
		atomCache.setFileParsingParams(params);
		StructureIO.setAtomCache(atomCache);

		Structure structure = StructureIO.getStructure("1PRP");
		List<BiologicalAssemblyTransformation> transformations =
				structure.getPDBHeader().getBioAssemblies().get(1).getTransforms();

		Structure bioUnitStructure = builder
				.rebuildQuaternaryStructure(BioAssemblyTools.getReducedStructure(structure), transformations, true, false);

		Atom[] reprAtoms = StructureTools.getRepresentativeAtomArray(bioUnitStructure);
		Point3d[] reprPoints = Calc.atomsToPoints(reprAtoms);
		Volume volume = new Volume();
		String[] resNames = Arrays.stream(reprAtoms).map(a -> a.getGroup().getPDBName()).toArray(String[]::new);
		volume.create(reprPoints, resNames);
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
}
