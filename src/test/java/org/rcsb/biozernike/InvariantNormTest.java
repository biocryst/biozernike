package org.rcsb.biozernike;

import org.junit.Test;
import org.rcsb.biozernike.volume.Volume;
import org.rcsb.biozernike.volume.ResidueVolumeCache;

import javax.vecmath.Point3d;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class InvariantNormTest {

	@Test
	public void invariants() {
		int nPoints = 1000;
		double coordsLimit = 96;
		Point3d[] points = new Point3d[nPoints];

		for (int i=0;i<nPoints;i++) {
			points[i] = new Point3d(Math.random()*coordsLimit,Math.random()*coordsLimit,Math.random()*coordsLimit);
		}

		List<List<Double>> invariantsScaled = new ArrayList<>();
		for (double gridWidth: ResidueVolumeCache.GRID_WIDTHS) {
			Volume volume = new Volume();
			volume.create(points, gridWidth);

			InvariantNorm normalization = new InvariantNorm(volume,15);
			invariantsScaled.add(normalization.getInvariants(2));
		}

		// check that invariantsMap are similar whatever grid scaling was used
		List<Double> invariantsRef = invariantsScaled.get(0);
		for(int i=1;i<invariantsScaled.size();i++) {
			List<Double> invariantsTest = invariantsScaled.get(i);
			for (int j=0;j<invariantsRef.size();j++) {
				double a = invariantsRef.get(j);
				double b = invariantsTest.get(j);
				assertTrue(Math.abs(a-b)<1 || Math.abs(a-b)/(a+b+1)<0.1);
			}
		}
	}
}