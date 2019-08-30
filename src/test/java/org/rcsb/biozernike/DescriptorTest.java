package org.rcsb.biozernike;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.quaternary.BioAssemblyTools;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.junit.Test;
import org.rcsb.biozernike.descriptor.Descriptor;
import org.rcsb.biozernike.descriptor.DescriptorConfig;
import org.rcsb.biozernike.descriptor.DescriptorMode;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

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
}
