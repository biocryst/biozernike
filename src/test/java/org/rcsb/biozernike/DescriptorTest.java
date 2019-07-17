package org.rcsb.biozernike;

import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureIO;
import org.biojava.nbio.structure.quaternary.BioAssemblyTools;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.junit.Test;

import java.util.ArrayList;
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

		long startTime = System.currentTimeMillis();
		for (int i=0;i<n_iterations;i++) {
			Descriptor ssd = new Descriptor(bioUnitStructure);
			List<Double> invariantsSearch = ssd.getMomentInvariants();
			tmp.add(invariantsSearch.get(0));
		}
		long endTime = System.currentTimeMillis();
		double duration = endTime - startTime;
		System.out.println("Performance: "+n_iterations/(duration/1000)+" 4HHBs per second");
	}

}
