package org.rcsb.biozernike.descriptor;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.*;

;
//TODO: split alignment and search
//TODO: unify arrays/lists
public class DescriptorConfig {
	private static final Logger logger = LoggerFactory.getLogger(DescriptorConfig.class);

	public int[] normOrders;
	public List<int[]> indicesZernike;
	public List<double[]> coefficientsZernike;

	public double[] coefficientsGeometry;

	public double[][] thresholdSets;
	public double referenceRadius;
	public double weightGeometry;

	public int maxOrderZernike;

	public int maxOrderZernikeAlign;
	public List<Map.Entry<Integer, Integer>> alignNormKeys;

	public DescriptorMode mode;

	private String propsFile = "";

	public DescriptorConfig(
			int maxOrderZernike,
			int[] normOrders
	) {
		this.maxOrderZernike = maxOrderZernike;
		this.normOrders = normOrders;
		this.mode = DescriptorMode.CALCULATE_RAW;
	}

	public DescriptorConfig(
			int maxOrderZernike,
			int[] normOrders,
			List<int[]> indicesZernike
) {
		this(maxOrderZernike, normOrders);
		this.indicesZernike = indicesZernike;
		this.mode = DescriptorMode.CALCULATE_PROCESSED;
	}

	public DescriptorConfig(
			int maxOrderZernike,
			int[] normOrders,
			List<int[]> indicesZernike,
			List<double[]> coefficientsZernike,
			double[] coefficientsGeometry,
			double weightGeometry,
			double referenceRadius,
			double[][] thresholdSets
) {
		this(maxOrderZernike, normOrders, indicesZernike);

		this.coefficientsZernike = coefficientsZernike;
		this.coefficientsGeometry = coefficientsGeometry;
		this.weightGeometry = weightGeometry;
		this.referenceRadius = referenceRadius;
		this.thresholdSets = thresholdSets;

		this.mode = DescriptorMode.COMPARE;
	}

	public DescriptorConfig(
			int maxOrderZernikeAlign,
			List<Map.Entry<Integer, Integer>> alignNormKeys
	) {
		this.maxOrderZernike = maxOrderZernikeAlign;
		this.maxOrderZernikeAlign = maxOrderZernikeAlign;
		this.alignNormKeys = alignNormKeys;
		this.mode = DescriptorMode.ALIGN;
	}

	public DescriptorConfig(
			int maxOrderZernike,
			int[] normOrders,
			List<int[]> indicesZernike,
			List<double[]> coefficientsZernike,
			double[] coefficientsGeometry,
			double weightGeometry,
			double referenceRadius,
			double[][] thresholdSets,
			int maxOrderZernikeAlign,
			List<Map.Entry<Integer, Integer>> alignNormKeys
	) {
		this(
				maxOrderZernike,
				normOrders,
				indicesZernike,
				coefficientsZernike,
				coefficientsGeometry,
				weightGeometry,
				referenceRadius,
				thresholdSets
		);
		this.maxOrderZernikeAlign = maxOrderZernikeAlign;
		this.alignNormKeys = alignNormKeys;
		this.mode = DescriptorMode.COMPARE_ALIGN;
	}

	public DescriptorConfig(String propsFile, DescriptorMode mode) {
		this.propsFile = propsFile;
		this.mode = mode;

		Properties props = new Properties();
		try {
			InputStream input = new FileInputStream(propsFile);
			props.load(input);
		} catch (IOException e) {
			logger.error("Invalid or missing config file {}", propsFile);
			throw new RuntimeException("Missing configuration '"+propsFile+"'. Can't continue.");
		}

		if (mode == DescriptorMode.CALCULATE_RAW || mode == DescriptorMode.CALCULATE_PROCESSED
				|| mode == DescriptorMode.COMPARE || mode == DescriptorMode.COMPARE_ALIGN) {
			normOrders = loadIntArrayField(props,"norm.orders.zernike");
			maxOrderZernike = loadIntegerField(props,"max.order.zernike");
		}

		if (mode == DescriptorMode.CALCULATE_PROCESSED || mode == DescriptorMode.COMPARE || mode == DescriptorMode.COMPARE_ALIGN) {
			indicesZernike = new ArrayList<>();
			for (int i = 0; i < normOrders.length; i++) {
				indicesZernike.add(loadIntArrayField(props, "selected.indices.zernike." + i));
			}
		}

		if (mode == DescriptorMode.COMPARE || mode == DescriptorMode.COMPARE_ALIGN) {
			coefficientsZernike = new ArrayList<>();
			for (int i=0; i<normOrders.length;i++) {
				coefficientsZernike.add(loadDoubleArrayField(props,"selected.coefficients.zernike."+i));
			}

			coefficientsGeometry = loadDoubleArrayField(props, "coefficients.geometry");
			weightGeometry = loadDoubleField(props, "weight.geometry");
			referenceRadius = loadDoubleField(props, "reference.radius");

			int numThresholdSets = loadIntegerField(props, "num.threshold.sets");
			thresholdSets = new double[numThresholdSets][2];
			for (int i=0;i<numThresholdSets;i++) {
				thresholdSets[i][0]=loadDoubleField(props, "threshold.geometry."+i);
				thresholdSets[i][1]=loadDoubleField(props, "threshold.zernike."+i);
			}
		}

		if (mode == DescriptorMode.ALIGN || mode == DescriptorMode.COMPARE_ALIGN) {
			maxOrderZernikeAlign = loadIntegerField(props, "max.order.zernike.align");
			int numNormKeys = loadIntegerField(props, "num.zernike.align.keys");
			alignNormKeys = new ArrayList<>();
			for (int i=0;i<numNormKeys;i++) {
				int indZero = loadIntegerField(props, "zernike.align.indzero."+i);
				int indReal = loadIntegerField(props, "zernike.align.indreal."+i);
				alignNormKeys.add(new AbstractMap.SimpleImmutableEntry<>(indZero, indReal));
			}
		}

		if (mode == DescriptorMode.ALIGN) {
			maxOrderZernike = maxOrderZernikeAlign;
		}

	}

	private double loadDoubleField(Properties props, String field) {
		String value = props.getProperty(field);
		double doubleValue;
		if (value == null || value.trim().equals("")) {
			logger.error("Field '{}' is not specified correctly in the config file {}", field, propsFile);
			throw new RuntimeException("Missing configuration '"+field+"'. Can't continue.");
		} else {
			logger.info("Using value '{}' for configuration field '{}'", value, field);
		}
		try {
			doubleValue = Double.parseDouble(value);
		} catch (NumberFormatException e) {
			logger.error("Could not parse double from specified value '{}' for property '{}'", value, field);
			throw new RuntimeException("Could not parse double from specified '"+field+"' property");
		}
		return doubleValue;
	}

	private int loadIntegerField(Properties props, String field) {
		String value = props.getProperty(field);
		int intValue;
		if (value == null || value.trim().equals("")) {
			logger.error("Field '{}' is not specified correctly in the config file {}", field, propsFile);
			throw new RuntimeException("Missing configuration '"+field+"'. Can't continue.");
		} else {
			logger.info("Using value '{}' for configuration field '{}'", value, field);
		}
		try {
			intValue = Integer.parseInt(value);
		} catch (NumberFormatException e) {
			logger.error("Could not parse integer from specified value '{}' for property '{}'", value, field);
			throw new RuntimeException("Could not parse integer from specified '"+field+"' property");
		}

		return intValue;
	}

	private double[] loadDoubleArrayField(Properties props, String field) {
		String value = props.getProperty(field);
		double[] doubleArrValue;
		if (value == null || value.trim().equals("")) {
			logger.error("Field '{}' is not specified correctly in the config file {}", field, propsFile);
			throw new RuntimeException("Missing configuration '"+field+"'. Can't continue.");
		} else {
			logger.info("Using value '{}' for configuration field '{}'", value, field);
		}
		String[] tokens = value.split(",\\s*");
		doubleArrValue = new double[tokens.length];
		for (int i=0; i<tokens.length; i++) {
			try {
				doubleArrValue[i] = Double.parseDouble(tokens[i]);
			} catch (NumberFormatException e) {
				logger.error("Could not parse double from specified value '{}' at index {} for property '{}'", tokens[i], i, field);
				throw new RuntimeException("Could not parse double from specified '"+field+"' property");
			}

		}
		return doubleArrValue;
	}

	private int[] loadIntArrayField(Properties props, String field) {
		String value = props.getProperty(field);
		int[] intArrValue;
		if (value == null || value.trim().equals("")) {
			logger.error("Field '{}' is not specified correctly in the config file {}", field, propsFile);
			throw new RuntimeException("Missing configuration '"+field+"'. Can't continue.");
		} else {
			logger.info("Using value '{}' for configuration field '{}'", value, field);
		}
		String[] tokens = value.split(",\\s*");
		intArrValue = new int[tokens.length];
		for (int i=0; i<tokens.length; i++) {
			try {
				intArrValue[i] = Integer.parseInt(tokens[i]);
			} catch (NumberFormatException e) {
				logger.error("Could not parse int from specified value '{}' at index {} for property '{}'", tokens[i], i, field);
				throw new RuntimeException("Could not parse int from specified '"+field+"' property");
			}
		}
		return intArrValue;
	}


}
