package org.rcsb.biozernike;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.quaternary.BioAssemblyTools;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.junit.Test;
import org.rcsb.biozernike.descriptor.Descriptor;
import org.rcsb.biozernike.descriptor.DescriptorConfig;

import javax.vecmath.Point3d;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;

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

		DescriptorConfig config = readConfig("src/test/resources/descriptor.properties");

		long startTime = System.currentTimeMillis();
		for (int i=0;i<n_iterations;i++) {
			Descriptor ssd = new Descriptor(reprPoints,resNames,config);
			List<Double> invariantsSearch = ssd.getMomentInvariants();
			tmp.add(invariantsSearch.get(0));
		}
		long endTime = System.currentTimeMillis();
		double duration = endTime - startTime;
		System.out.println("Performance: "+n_iterations/(duration/1000)+" 4HHBs per second");
	}

	private DescriptorConfig readConfig(String filename) throws IOException {
		Properties props = new Properties();
		InputStream input = new FileInputStream(filename);
		props.load(input);

		double thresholdZernike =  Double.parseDouble(props.getProperty("threshold.zernike"));
		double thresholdGeometry = Double.parseDouble(props.getProperty("threshold.geometry"));
		int [] searchIndicesZernike = loadIntArrayField(props, "search.indices.zernike");
		double[] searchCoefficientsZernike = loadDoubleArrayField(props, "search.coefficients.zernike");
		double[] searchCoefficientsGeometry = loadDoubleArrayField(props, "search.coefficients.geometry");

		return new DescriptorConfig(
				20,
				6,
				new int[]{2,3,4,5},
				searchIndicesZernike,
				searchCoefficientsZernike,
				searchCoefficientsGeometry,
				thresholdZernike,
				thresholdGeometry
		);
	}

	private static int[] loadIntArrayField(Properties props, String field) {
		String[] tokens = props.getProperty(field).split(",\\s*");
		int[] intArrValue = new int[tokens.length];
		for (int i=0; i<tokens.length; i++) {
			intArrValue[i] = Integer.parseInt(tokens[i]);
		}
		return intArrValue;
	}

	private static double[] loadDoubleArrayField(Properties props, String field) {
		String[] tokens = props.getProperty(field).split(",\\s*");
		double[] intArrValue = new double[tokens.length];
		for (int i=0; i<tokens.length; i++) {
			intArrValue[i] = Double.parseDouble(tokens[i]);
		}
		return intArrValue;
	}

}
