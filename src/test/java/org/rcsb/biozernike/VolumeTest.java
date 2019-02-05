package org.rcsb.biozernike;

import org.junit.Test;
import org.rcsb.biozernike.volume.Volume;
import org.rcsb.biozernike.volume.ResidueVolumeCache;
import org.rcsb.biozernike.volume.VolumeConstants;

import javax.vecmath.Point3d;

import java.util.ArrayList;
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
}