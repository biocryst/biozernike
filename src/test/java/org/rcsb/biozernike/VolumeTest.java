package org.rcsb.biozernike;

import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.StructureFiletype;
import org.biojava.nbio.structure.quaternary.BioAssemblyTools;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyBuilder;
import org.biojava.nbio.structure.quaternary.BiologicalAssemblyTransformation;
import org.junit.Test;
import org.rcsb.biozernike.volume.*;
import org.rcsb.biozernike.zernike.ZernikeMoments;

import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.*;

public class VolumeTest {

	@Test
	public void mass(){
		Point3d[] points1 = {new Point3d(0,0,0)};
		String[] resNames1 = {"ALA"};

		Volume volumeSmall = new Volume();
		volumeSmall.create(points1,resNames1);

		// one residue in the volume, the weight should match
		assertEquals(VolumeConstants.getWeight("ALA"),volumeSmall.getResiduesNominalWeight(),0);
		// very small volume -- minimal grid width
		assertEquals(ResidueVolumeCache.MIN_GRID_WIDTH,volumeSmall.getGridWidth(),0);

		double farAway = volumeSmall.getMaxVolumeSize()*ResidueVolumeCache.MAX_GRID_WIDTH;
		Point3d[] points2 = {new Point3d(0,0,0),
				new Point3d(farAway,farAway,farAway)};
		String[] resNames2 = {"ALA","ALA"};

		Volume volumeLarge = new Volume();
		volumeLarge.create(points2,resNames2);


		assertEquals(VolumeConstants.getWeight("ALA")*points2.length,
				volumeLarge.getResiduesNominalWeight(),0);

		// very large volume -- maximal grid width
		assertEquals(ResidueVolumeCache.MAX_GRID_WIDTH,volumeLarge.getGridWidth(),0);

		// voxel values are scaled by the grid width cubed to preserve relative moment values
		assertEquals(Math.pow(volumeLarge.getGridWidth()/volumeSmall.getGridWidth(),3),
				points2.length*volumeSmall.getVolumeMass()/volumeLarge.getVolumeMass(),0.00001);

	}

	@Test
	public void radius(){
		Point3d[] points = {new Point3d(0,0,0)};
		String[] resNames = {"ALA"};

		Volume volume = new Volume();
		volume.create(points,resNames);
		volume.normalize();

		// volume radius corresponds to the residue radius (+/- half angstroem)
		assertEquals(ResidueVolumeCache.get(volume.getGridWidth(),"ALA").size/2.0,volume.getRadiusMaxVolume(),0.5);

		// max radius is larger than scaled Rg
		assertTrue(volume.getRadiusMaxVolume()>volume.getRadiusVarVolume());
		// scaled Rg is within 10% of the max radius
		assertEquals(volume.getRadiusMaxVolume()/volume.getRadiusVarVolume(),1,0.1);
		// trimmed mass is less then 1%
		assertEquals(volume.getVolumeMass()/volume.getOriginalVolumeMass(),1,0.01);
	}

	@Test
	public void gridding(){
		List<Volume> volumes = new ArrayList<>();
		double shift = 10+Math.random();
		for(double coordStart=0;coordStart<1;coordStart+=0.1) {
			Point3d[] points = {new Point3d(coordStart,coordStart,coordStart),
					new Point3d(coordStart+shift,coordStart+shift,coordStart+shift)};
			String[] resNames = {"G","G"};

			Volume volume = new Volume();
			volume.create(points,resNames);
			int[] dimensions = volume.getDimensions();
			double[] center = volume.getCenterVolume();
			double[] dimCenter = {dimensions[0]/2.0,dimensions[1]/2.0,dimensions[2]/2.0};
			assertArrayEquals(dimCenter,center,1);
			volumes.add(volume);
		}
		// volume depends on the relative distances of atoms
		Volume reprVolume = volumes.get(0);
		for(int i=1;i<volumes.size();i++) {
			assertArrayEquals(reprVolume.getVoxelArray(),volumes.get(i).getVoxelArray(),0.00001);
		}
	}

	@Test
	public void testVolumeBounds() throws Exception {
		Point3d[] points = {new Point3d(-10.123,-10.123,-10.123), new Point3d(10.123,10.123,10.123)};
		String[] resNames = {"DG","DG"}; // largest residues at the corners

		Volume volume = new Volume();
		volume.create(points, resNames);
	}

	@Test
	public void testDnaNegIndexIssue() throws Exception {
		BiologicalAssemblyBuilder builder = new BiologicalAssemblyBuilder();

		// this issue happens only when reading from mmCIF. Comment the following out to read from MMTF and see it works
		// --------------------
		AtomCache atomCache = new AtomCache();
		atomCache.setFiletype(StructureFiletype.CIF);
		FileParsingParameters params = new FileParsingParameters();
		params.setParseBioAssembly(true);
		atomCache.setFileParsingParams(params);
		StructureIO.setAtomCache(atomCache);
		// --------------------

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
	}

	@Test
	public void testLigand() throws Exception {
		Structure structure1 = StructureIO.getStructure("6U9R");
		Structure structure2 = StructureIO.getStructure("6U9N");
		List<Chain> chains1 = structure1.getNonPolyChains();
		List<Chain> chains2 = structure2.getNonPolyChains();
		Atom[] ligand1 = StructureTools.getAllAtomArray(chains1.get(1));
		Atom[] ligand2 = StructureTools.getAllAtomArray(chains2.get(0));

		String[] atomNames1 = Arrays.stream(ligand1).map(a -> a.getElement().name()).toArray(String[]::new);
		String[] atomNames2 = Arrays.stream(ligand2).map(a -> a.getElement().name()).toArray(String[]::new);

		Point3d[] points1 = Calc.atomsToPoints(ligand1);
		Point3d[] points2 = Calc.atomsToPoints(ligand2);

		Volume volume1 = new Volume();
		Volume volume2 = new Volume();

		volume1.create(points1,atomNames1);
		volume2.create(points2,atomNames2);
		volume1.setRadiusVarMult(2);
		volume2.setRadiusVarMult(2);

		VolumeIO.write(volume1,"D:/PT/Ligands/volume1_original.map",MapFileType.CCP4);
		VolumeIO.write(volume2,"D:/PT/Ligands/volume2_original.map",MapFileType.CCP4);

		List<InvariantNorm> zc = new ArrayList<>();

		InvariantNorm normalization1 = new InvariantNorm(volume1, 15);
		InvariantNorm normalization2 = new InvariantNorm(volume2, 15);

		zc.add(normalization1);
		zc.add(normalization2);
		AlignmentResult alignmentResult = RotationAlignment.alignMultiple(zc, 15, false);

		Matrix4d rot1 = alignmentResult.getTransforms().get(0);
		Matrix4d rot2 = alignmentResult.getTransforms().get(1);

		Calc.transform(ligand1,rot1);
		Calc.transform(ligand2,rot2);

		Point3d[] points1_transformed = Calc.atomsToPoints(ligand1);
		Volume volume1_transformed = new Volume();
		volume1_transformed.create(points1_transformed,atomNames1);
		VolumeIO.write(volume1_transformed,"D:/PT/Ligands/volume1_transformed.map",MapFileType.CCP4);

		Point3d[] points2_transformed = Calc.atomsToPoints(ligand2);
		Volume volume2_transformed = new Volume();
		volume2_transformed.create(points2_transformed,atomNames2);
		VolumeIO.write(volume2_transformed,"D:/PT/Ligands/volume2_transformed.map",MapFileType.CCP4);
	}

	@Test
	public void testReconstruction() throws Exception {
		Structure structure1 = StructureIO.getStructure("6U9R");
		List<Chain> chains1 = structure1.getNonPolyChains();
		Atom[] ligand1 = StructureTools.getAllAtomArray(chains1.get(1));

		String[] atomNames1 = Arrays.stream(ligand1).map(a -> a.getElement().name()).toArray(String[]::new);

		Point3d[] points1 = Calc.atomsToPoints(ligand1);

		Volume volume1 = new Volume();

		volume1.create(points1, atomNames1);

		VolumeIO.write(volume1, "D:/PT/Ligands/orders/original.map", MapFileType.MRC);

		ZernikeMoments zm = new ZernikeMoments(volume1, 15);
//		Volume volume_rec = ZernikeMoments.reconstructVolume(zm, 32, 15, false, true);
//		VolumeIO.write(volume_rec, "D:/PT/Ligands/volume1_rec_orig.map", MapFileType.MRC);

		for (int maxN=0;maxN<=15;maxN++) {
			Volume volume_rec2 = ZernikeMoments.reconstructVolume(zm, 32, maxN, true, false);
			VolumeIO.write(volume_rec2, "D:/PT/Ligands/orders/"+maxN+".map", MapFileType.MRC);
		}
	}


//	@Test
//	public void testSaving() throws Exception {
//		Structure structure1 = StructureIO.getStructure("6U9R");
//		List<Chain> chains1 = structure1.getNonPolyChains();
//		Atom[] ligand1 = StructureTools.getAllAtomArray(chains1.get(1));
//
//		String[] atomNames1 = Arrays.stream(ligand1).map(a -> a.getElement().name()).toArray(String[]::new);
//
//		Point3d[] points1 = Calc.atomsToPoints(ligand1);
//
//		Volume volume1 = new Volume();
//
//		volume1.create(points1, atomNames1);
//
//		VolumeIO.write(volume1, "D:/PT/Ligands/volume1_original.map", MapFileType.MRC);
//
//		ZernikeMoments zm = new ZernikeMoments(volume1, 20);
//		List<Double> mom = ZernikeMoments.flattenMomentsDouble(zm.getOriginalMoments());
//
//
//	}
	@Test
	public void testInterface() throws Exception {
		Structure structure1 = StructureIO.getStructure("6U9R");
		List<Chain> chains = structure1.getPolyChains();
		List<Chain> chains1 = structure1.getNonPolyChains();

		Atom[] ligand1 = StructureTools.getAllAtomArray(chains1.get(1));
		String[] atomNames1 = Arrays.stream(ligand1).map(a -> a.getElement().name()).toArray(String[]::new);
		Point3d[] points1 = Calc.atomsToPoints(ligand1);

		Atom[] protein = StructureTools.getAllAtomArray(chains.get(0));
		String[] atomNames = Arrays.stream(protein).map(a -> a.getElement().name()).toArray(String[]::new);
		Point3d[] points = Calc.atomsToPoints(protein);


		Volume ligand_interface = new Volume();

		ligand_interface.createFromInterface(points, points1, atomNames, atomNames1,0.25);

		VolumeIO.write(ligand_interface, "D:/PT/Ligands/6u9r_interface.map", MapFileType.CCP4);


	}
	@Test
	public void testInterfaceProtein() throws Exception {
		Structure structure1 = StructureIO.getStructure("7O29");
		List<Chain> chains = structure1.getPolyChains();

		Atom[] chain1 = StructureTools.getRepresentativeAtomArray(chains.get(0));
		String[] resNames1 = Arrays.stream(chain1).map(a->a.getGroup().getPDBName()).toArray(String[]::new);
		Point3d[] points1 = Calc.atomsToPoints(chain1);

		Volume chain1_volume = new Volume();
		chain1_volume.create(points1, resNames1,0.25);
		VolumeIO.write(chain1_volume, "D:/PT/Ligands/7o29_chain1.map", MapFileType.CCP4);

		Atom[] protein = StructureTools.getRepresentativeAtomArray(chains.get(1));
		String[] resNames = Arrays.stream(protein).map(a->a.getGroup().getPDBName()).toArray(String[]::new);
		Point3d[] points = Calc.atomsToPoints(protein);

		Volume chain2_volume = new Volume();
		chain2_volume.create(points, resNames,0.25);
		VolumeIO.write(chain2_volume, "D:/PT/Ligands/7o29_chain2.map", MapFileType.CCP4);


		Volume chain_interface = new Volume();
		chain_interface.createFromInterface(points, points1, resNames, resNames1,0.25);
		VolumeIO.write(chain_interface, "D:/PT/Ligands/7o29_interface.map", MapFileType.CCP4);



	}

	@Test
	public void testReadWrite() throws Exception {
		Structure structure = StructureIO.getStructure("6U9R");
		List<Chain> chains_nonpoly = structure.getNonPolyChains();

		Atom[] ligandAtoms = StructureTools.getAllAtomArray(chains_nonpoly.get(1));
		String[] ligandAtomNames = Arrays.stream(ligandAtoms).map(a -> a.getElement().name()).toArray(String[]::new);
		Point3d[] ligandPoints = Calc.atomsToPoints(ligandAtoms);

		// volumes for ligand and interface, save as is
		Volume ligandVolume = new Volume();
		ligandVolume.create(ligandPoints, ligandAtomNames,0.25);

		VolumeIO.write(ligandVolume, "D:/PT/Ligands/loading/ligand_original.map", MapFileType.MRC);
		InputStream fs = Files.newInputStream(Paths.get("D:/PT/Ligands/loading/ligand_original.map"));

		Volume readVolume = VolumeIO.read(fs, MapFileType.MRC, false);
		VolumeIO.write(readVolume, "D:/PT/Ligands/loading/ligand_rewrite.map", MapFileType.MRC);


		ZernikeMoments ligandZm = new ZernikeMoments(ligandVolume, 15);
		Volume ligandVolumeRec = ZernikeMoments.reconstructVolume(ligandZm,32, 15, true, false);

		VolumeIO.write(ligandVolumeRec, "D:/PT/Ligands/loading/ligand_rec.map", MapFileType.MRC);

	}

	@Test
	public void testVectors() throws Exception {
		Structure structure = StructureIO.getStructure("6U9R");
		List<Chain> chains_poly = structure.getPolyChains();
		List<Chain> chains_nonpoly = structure.getNonPolyChains();

		Atom[] ligandAtoms = StructureTools.getAllAtomArray(chains_nonpoly.get(1));
		String[] ligandAtomNames = Arrays.stream(ligandAtoms).map(a -> a.getElement().name()).toArray(String[]::new);
		Point3d[] ligandPoints = Calc.atomsToPoints(ligandAtoms);

		Atom[] proteinAtoms = StructureTools.getAllAtomArray(chains_poly.get(0));
		String[] proteinAtomNames = Arrays.stream(proteinAtoms).map(a -> a.getElement().name()).toArray(String[]::new);
		Point3d[] proteinPoints = Calc.atomsToPoints(proteinAtoms);

		// volumes for ligand and interface, save as is
		Volume ligandVolume = new Volume();
		ligandVolume.create(ligandPoints, ligandAtomNames,0.25);

		VolumeIO.write(ligandVolume, "D:/PT/Ligands/orientation/ligand_original.map", MapFileType.MRC);
		ligandVolume.normalize();

		ZernikeMoments ligandZm = new ZernikeMoments(ligandVolume, 15);
		Volume ligandVolumeRec = ZernikeMoments.reconstructVolume(ligandZm,32, 15, true, false);
		VolumeIO.write(ligandVolumeRec, "D:/PT/Ligands/orientation/ligand_original_rec.map", MapFileType.MRC);

		InvariantNorm ligandNormalization = new InvariantNorm(ligandVolume, 15);
		int[] indsPositive = {10,14,17};
		List<MomentTransform> ligandTransforms = ligandNormalization.getConstrainedNormalizationSolution(2,2, indsPositive);
		MomentTransform ligandTransform = ligandTransforms.get(0);
		Matrix3d ligandRotation = ligandTransform.rotation();
		Vector3d ligandCenter = ligandNormalization.getCenter();
		System.out.print(ligandRotation);
		System.out.print(ligandCenter);

//		VolumeIO.write(ligandVolume, "D:/PT/Ligands/orientation/ligand_normalized.map", MapFileType.MRC);

		ligandZm = new ZernikeMoments(ligandTransform.getMoments(),true);
		ligandVolumeRec = ZernikeMoments.reconstructVolume(ligandZm,32, 15, true, false);

		VolumeIO.write(ligandVolumeRec, "D:/PT/Ligands/orientation/ligand_transformed_rec.map", MapFileType.MRC);


		Volume interfaceVolume = new Volume();
		interfaceVolume.createFromInterface(proteinPoints, ligandPoints, proteinAtomNames, ligandAtomNames,0.25);

		VolumeIO.write(interfaceVolume, "D:/PT/Ligands/orientation/interface_original.map", MapFileType.MRC);

		interfaceVolume.normalize();
		InvariantNorm interfaceNormalization = new InvariantNorm(interfaceVolume, 15);

		ZernikeMoments interfaceZm = new ZernikeMoments(interfaceVolume, 15);
		Volume interfaceVolumeRec = ZernikeMoments.reconstructVolume(interfaceZm,32, 15, true, false);
		VolumeIO.write(interfaceVolumeRec, "D:/PT/Ligands/orientation/interface_original_rec.map", MapFileType.MRC);

		List<MomentTransform> interfaceTransforms = interfaceNormalization.getConstrainedNormalizationSolution(2,2, indsPositive);
		MomentTransform interfaceTransform = interfaceTransforms.get(0);
		Matrix3d interfaceRotation = interfaceTransform.rotation();
		Vector3d interfaceCenter = interfaceNormalization.getCenter();

		System.out.print(interfaceRotation);
		System.out.print(interfaceCenter);


		interfaceZm = new ZernikeMoments(interfaceTransform.getMoments(),true);
		interfaceVolumeRec = ZernikeMoments.reconstructVolume(interfaceZm,32, 15, true, false);

		VolumeIO.write(interfaceVolumeRec, "D:/PT/Ligands/orientation/interface_transformed_rec.map", MapFileType.MRC);

	}
	@Test
	public void testCombined() throws Exception {
		// find center and positive rotation of the pocket
		String pocketFileName ="D:\\PT\\Ligands\\Combined\\1a4k_pocket.pdb";
		String ligandFileName ="D:\\PT\\Ligands\\Combined\\1a4k_ligand.pdb";

		Structure structurePocket = StructureIO.getStructure(pocketFileName);
		Atom[] pocket = StructureTools.getAllAtomArray(structurePocket);
		System.out.println(" pocket atoms: "+pocket.length);

		String[] pocketAtomNames = Arrays.stream(pocket).map(a -> a.getElement().name()).toArray(String[]::new);
		Point3d[] pocketPoints = Calc.atomsToPoints(pocket);
		Volume pocketVolume = new Volume();
		pocketVolume.setRadiusVarMult(2.0);
		pocketVolume.create(pocketPoints,pocketAtomNames,0.25);


		InvariantNorm pocketNormalization = new InvariantNorm(pocketVolume,15);
		Volume pocketVolumeRec = ZernikeMoments.reconstructVolume(pocketNormalization.getMoments(),32, 15, true, false);
		VolumeIO.write(pocketVolumeRec, "D:\\PT\\Ligands\\Combined\\pocket_original_rec.map", MapFileType.MRC);


		double radiusVar = pocketVolume.getRadiusVarReal();
		int[] indsPositive = {10,14,17};

		List<MomentTransform> tr = pocketNormalization.getConstrainedNormalizationSolution(2,2, indsPositive);

		MomentTransform pocketTransform = tr.get(0);
		Matrix3d pocketRotation = pocketTransform.rotation();
		Vector3d pocketCenter = pocketNormalization.getCenter();

		pocketVolumeRec = ZernikeMoments.reconstructVolume(new ZernikeMoments(pocketTransform.getMoments(), true),32, 15, true, false);
		VolumeIO.write(pocketVolumeRec, "D:\\PT\\Ligands\\Combined\\pocket_aligned_rec.map", MapFileType.MRC);

		// apply this to the ligand ,


		Structure structureLigand = StructureIO.getStructure(ligandFileName);
		Atom[] ligand = StructureTools.getAllAtomArray(structureLigand);
		System.out.println(" ligand  atoms: "+ligand.length);

		String[] ligandAtomNames = Arrays.stream(ligand).map(a -> a.getElement().name()).toArray(String[]::new);
		Point3d[] ligandPoints = Calc.atomsToPoints(ligand);
		Volume ligandVolume = new Volume();
		ligandVolume.create(ligandPoints, ligandAtomNames,0.25);


		ligandVolume.setRadiusVarMult(2.0);
		ligandVolume.setCenterReal(pocketVolume.getCenterReal());
//		VolumeIO.write(ligandVolume, "D:\\PT\\Ligands\\Combined\\ligand_original_shifted.map", MapFileType.MRC);

		ZernikeMoments ligandZm = new ZernikeMoments(ligandVolume, 15);

		// ligand radius?

		Volume ligandVolumeRec = ZernikeMoments.reconstructVolume(ligandZm,32, 15, true, false);
		VolumeIO.write(ligandVolumeRec, "D:\\PT\\Ligands\\Combined\\ligand_original_shifted_rec.map", MapFileType.MRC);

		ligandVolume.updateCenter();
		ligandZm = new ZernikeMoments(ligandVolume, 15);

		ligandVolumeRec = ZernikeMoments.reconstructVolume(ligandZm,32, 15, true, false);
		VolumeIO.write(ligandVolumeRec, "D:\\PT\\Ligands\\Combined\\ligand_original_rec.map", MapFileType.MRC);

		InvariantNorm ligandNormalization = new InvariantNorm(ligandZm);

		MomentTransform transform = new MomentTransform();
		transform.setRotation(pocketTransform.getA(), pocketTransform.getB());
		transform.setMoments(ligandNormalization.rotate(pocketTransform.getA(), pocketTransform.getB()));

		ZernikeMoments ligandPocketTr = new ZernikeMoments(transform.getMoments(), true);

		ligandVolumeRec = ZernikeMoments.reconstructVolume(ligandPocketTr,32, 15, true, false);
		VolumeIO.write(ligandVolumeRec, "D:\\PT\\Ligands\\Combined\\ligand_aligned_rec.map", MapFileType.MRC);



	}




}