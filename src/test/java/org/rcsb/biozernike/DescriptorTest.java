package org.rcsb.biozernike;

import org.biojava.nbio.structure.*;
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
import java.util.List;

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

		DescriptorConfig config = new DescriptorConfig("src/test/resources/descriptor.properties", DescriptorMode.COMPARE_ALIGN);

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

		Structure structure = StructureIO.getStructure("4HHB");
		List<BiologicalAssemblyTransformation> transformations =
				structure.getPDBHeader().getBioAssemblies().get(1).getTransforms();

		Structure bioUnitStructure = builder
				.rebuildQuaternaryStructure(BioAssemblyTools.getReducedStructure(structure), transformations, true, false);

		Atom[] reprAtoms = StructureTools.getRepresentativeAtomArray(bioUnitStructure);
		Point3d[] reprPoints = Calc.atomsToPoints(reprAtoms);
		Volume volume = new Volume();
		volume.create(reprPoints);
		ZernikeMoments zernikeMoments1 = new ZernikeMoments(volume,10);

		List<List<List<Complex>>> originalMoments = zernikeMoments1.getOriginalMoments();
		List<Complex> originalMomentsFlatComplex = ZernikeMoments.flattenMomentsComplex(originalMoments);
		List<Double> originalMomentsFlatDouble = ZernikeMoments.flattenMomentsDouble(originalMoments);

		List<List<List<Complex>>> originalMomentsRestoredComplex = ZernikeMoments.unFlattenMomentsComplex(originalMomentsFlatComplex,zernikeMoments1.getMaxOrder());
		List<List<List<Complex>>> originalMomentsRestoredDouble = ZernikeMoments.unFlattenMomentsDouble(originalMomentsFlatDouble,zernikeMoments1.getMaxOrder());

		assert originalMoments.equals(originalMomentsRestoredComplex);
		assert originalMoments.equals(originalMomentsRestoredDouble);

		ZernikeMoments zernikeMoments2 = new ZernikeMoments(zernikeMoments1.getOriginalMomentsUnscaled(), false);
		ZernikeMoments zernikeMoments3 = new ZernikeMoments(zernikeMoments1.getOriginalMoments(), true);

		List<Double> flatMoments1 = ZernikeMoments.flattenMomentsDouble(zernikeMoments1.getOriginalMoments());
		List<Double> flatMoments2 = ZernikeMoments.flattenMomentsDouble(zernikeMoments2.getOriginalMoments());
		List<Double> flatMoments3 = ZernikeMoments.flattenMomentsDouble(zernikeMoments3.getOriginalMoments());

		assert flatMoments1.equals(flatMoments2);
		assert flatMoments1.equals(flatMoments3);
	}
}
