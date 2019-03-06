package org.rcsb.biozernike;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Matrix3d;
import java.util.*;
import java.util.stream.Collectors;

public class RotationAlignment {
	private static final Logger logger = LoggerFactory.getLogger(RotationAlignment.class);

	public static List<Matrix3d> alignMultiple(List<InvariantNorm> invariantNorms) {
		return alignMultiple(invariantNorms, 5, false);
	}

	public static List<Matrix3d> alignMultiple(List<InvariantNorm> invariantNorms, int maxNormInd, boolean useExistingTransforms) {

		if (invariantNorms.size() < 2) {
			return new ArrayList<>();
		}

		Set<Map.Entry<Integer, Integer>> commonNormKeys;
		if (useExistingTransforms) {
			commonNormKeys = new HashSet<>(invariantNorms.get(0).getTransformsMap().keySet());
			for (int i = 1; i < invariantNorms.size(); i++) {
				commonNormKeys.retainAll(invariantNorms.get(i).getTransformsMap().keySet());
			}

		} else {
			commonNormKeys = new HashSet<>(invariantNorms.get(0).populateNormalisations(maxNormInd));
			for (int i = 1; i < invariantNorms.size(); i++) {
				commonNormKeys.retainAll(invariantNorms.get(i).populateNormalisations(maxNormInd));
			}
		}

		List<Map.Entry<Integer, Integer>> useNormKeys = commonNormKeys.stream().
				filter(k -> k.getKey() <= maxNormInd && k.getValue() <= maxNormInd).
				collect(Collectors.toList());

		return alignMultiple(invariantNorms, useNormKeys);
	}


	public static List<Matrix3d> alignMultiple(List<InvariantNorm> invariantNorms,
	                                           List<Map.Entry<Integer, Integer>> normKeysArr) {

		if (invariantNorms.size() < 2) {
			return new ArrayList<>();
		}

		List<MomentTransform> rotations = alignCombinatorial(invariantNorms, normKeysArr);

		List<Matrix3d> transforms = new ArrayList<>();
		rotations.forEach(r -> transforms.add(r.rotation()));

		return transforms;
	}


	private static List<MomentTransform> alignCombinatorial(List<InvariantNorm> invariantNorms,
	                                                        List<Map.Entry<Integer, Integer>> normKeysArr) {

		int nStructures = invariantNorms.size();

		List<List<MomentTransform>> altSolutions = new ArrayList<>();
		List<Double> altResiduals = new ArrayList<>();
		Map<Integer, Integer> invCentroidInd = new HashMap<>();

		for (Map.Entry<Integer, Integer> normKey : normKeysArr) {
			int indZero = normKey.getKey();
			int indReal = normKey.getValue();

			if (indReal % 2 == 1 && indReal == indZero) {
				continue;
			}

			int indRefStructure = 0;

			if (nStructures > 2) {
				if (!invCentroidInd.containsKey(indZero)) {
					indRefStructure = invariantCentroid(invariantNorms, indZero);
					invCentroidInd.put(indZero, indRefStructure);
				} else {
					indRefStructure = invCentroidInd.get(indZero);
				}
			}


			List<MomentTransform> solution = new ArrayList<>();
			MomentTransform rotationRef = invariantNorms.get(indRefStructure).getNormalizationSolutions(indZero, indReal).get(0);

			double solutionSumResiduals = 0;

			for (int indStructure = 0; indStructure < invariantNorms.size(); indStructure++) {

				if (indStructure == indRefStructure) {
					solution.add(rotationRef);
					continue;
				}

				List<MomentTransform> rotationsAln = invariantNorms.get(indStructure).getNormalizationSolutions(indZero, indReal);

				double minResidual = Double.POSITIVE_INFINITY;
				MomentTransform bestRotation = null;

				for (MomentTransform rotationAln : rotationsAln) {
					double residual = rotationRef.distanceTo(rotationAln);

					if (residual < minResidual) {
						minResidual = residual;
						bestRotation = rotationAln;
					}
				}

				solutionSumResiduals += minResidual;
				solution.add(bestRotation);
			}

			altResiduals.add(solutionSumResiduals);
			altSolutions.add(solution);

		}

		int indSel = argmin(altResiduals);
		logger.info("Selected key: {}, {}. Residual: {}", normKeysArr.get(indSel).getKey(), normKeysArr.get(indSel).getValue(), altResiduals.get(indSel));
		return altSolutions.get(argmin(altResiduals));
	}

	private static int invariantCentroid(List<InvariantNorm> invariantNorms, int indZero) {
		int nStructures = invariantNorms.size();
		double[] distances = new double[nStructures];
		for (int i = 0; i < invariantNorms.size() - 1; i++) {
			InvariantNorm structure1 = invariantNorms.get(i);
			for (int j = i + 1; j < invariantNorms.size(); j++) {
				InvariantNorm structure2 = invariantNorms.get(j);
				double distance = structure1.compareInvariants(structure2, indZero);
				distances[i] += distance;
				distances[j] += distance;
			}
		}

		return argmin(distances);
	}

	private static int argmin(double[] arr) {
		int indMin = 0;
		double minValue = arr[0];

		for (int i = 1; i < arr.length; i++) {
			if (arr[i] < minValue) {
				minValue = arr[i];
				indMin = i;
			}
		}
		return indMin;
	}

	private static int argmin(List<Double> arr) {
		int indMin = 0;
		double minValue = arr.get(0);

		for (int i = 1; i < arr.size(); i++) {
			if (arr.get(i) < minValue) {
				minValue = arr.get(i);
				indMin = i;
			}
		}
		return indMin;
	}

}
