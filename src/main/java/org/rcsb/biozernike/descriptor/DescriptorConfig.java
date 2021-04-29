package org.rcsb.biozernike.descriptor;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.util.*;


//TODO: split alignment and search
//TODO: unify arrays/lists

/**
 * Configuration needed for BioZernike descriptors calculation and processing.
 *
 * @author Dmytro Guzenko
 */
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

	public EnumSet<DescriptorMode> mode;

	public boolean withCovEigenValsInGeom = false;

	public DescriptorConfig(
			int maxOrderZernike,
			int[] normOrders
	) {
		this.maxOrderZernike = maxOrderZernike;
		this.normOrders = normOrders;
		this.mode = EnumSet.of(DescriptorMode.CALCULATE_RAW);
	}

	public DescriptorConfig(
			int maxOrderZernike,
			int maxOrderZernikeAlign,
			int[] normOrders
	) {
		this(maxOrderZernike,normOrders);
		this.maxOrderZernikeAlign = maxOrderZernikeAlign;
		this.mode.add(DescriptorMode.ALIGN);
	}

	public DescriptorConfig(
			int maxOrderZernike,
			int[] normOrders,
			List<int[]> indicesZernike
	) {
		this(maxOrderZernike, normOrders);
		this.indicesZernike = indicesZernike;
		this.mode.add(DescriptorMode.CALCULATE_PROCESSED);
	}

	public DescriptorConfig(
			int maxOrderZernike,
			int maxOrderZernikeAlign,
			int[] normOrders,
			List<int[]> indicesZernike
	) {
		this(maxOrderZernike, maxOrderZernikeAlign, normOrders);
		this.indicesZernike = indicesZernike;
		this.mode.add(DescriptorMode.CALCULATE_PROCESSED);
	}

	public DescriptorConfig(
			int maxOrderZernike,
			int[] normOrders,
			List<int[]> indicesZernike,
			List<double[]> coefficientsZernike,
			double[] coefficientsGeometry,
			double weightGeometry,
			double referenceRadius,
			double[][] thresholdSets) {
		this(maxOrderZernike, normOrders, indicesZernike);

		this.coefficientsZernike = coefficientsZernike;
		this.coefficientsGeometry = coefficientsGeometry;
		this.weightGeometry = weightGeometry;
		this.referenceRadius = referenceRadius;
		this.thresholdSets = thresholdSets;

		this.mode.add(DescriptorMode.COMPARE);
	}

	public DescriptorConfig(
			int maxOrderZernike,
			int maxOrderZernikeAlign,
			int[] normOrders,
			List<int[]> indicesZernike,
			List<double[]> coefficientsZernike,
			double[] coefficientsGeometry,
			double weightGeometry,
			double referenceRadius,
			double[][] thresholdSets) {
		this(maxOrderZernike, maxOrderZernikeAlign, normOrders, indicesZernike);

		this.coefficientsZernike = coefficientsZernike;
		this.coefficientsGeometry = coefficientsGeometry;
		this.weightGeometry = weightGeometry;
		this.referenceRadius = referenceRadius;
		this.thresholdSets = thresholdSets;

		this.mode.add(DescriptorMode.COMPARE);
	}

	public DescriptorConfig(InputStream is, EnumSet<DescriptorMode> mode) throws IOException {
		this.mode = mode;
		ensureConsistentMode(this.mode);
		loadConfigs(is);
	}

	private void loadConfigs(InputStream is) throws IOException {

		Properties props = new Properties();
		props.load(is);

		if (mode.contains(DescriptorMode.CALCULATE_RAW)) {
			normOrders = loadIntArrayField(props,"norm.orders.zernike");
			maxOrderZernike = loadIntegerField(props,"max.order.zernike");
		}

		if (mode.contains(DescriptorMode.CALCULATE_PROCESSED)) {
			indicesZernike = new ArrayList<>();
			for (int i = 0; i < normOrders.length; i++) {
				indicesZernike.add(loadIntArrayField(props, "selected.indices.zernike." + i));
			}
		}

		if (mode.contains(DescriptorMode.COMPARE)) {
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

		if (mode.contains(DescriptorMode.ALIGN)) {
			maxOrderZernikeAlign = loadIntegerField(props, "max.order.zernike.align");
		}
	}

	private void ensureConsistentMode(EnumSet<DescriptorMode> mode) {
		if (mode.contains(DescriptorMode.COMPARE)) {
			mode.add(DescriptorMode.CALCULATE_PROCESSED);
		}
		if (mode.contains(DescriptorMode.CALCULATE_PROCESSED)) {
			mode.add(DescriptorMode.CALCULATE_RAW);
		}
		if(mode.contains(DescriptorMode.ALIGN)) {
			mode.add(DescriptorMode.CALCULATE_RAW);
		}
	}

	private double loadDoubleField(Properties props, String field) {
		String value = props.getProperty(field);
		double doubleValue;
		if (value == null || value.trim().equals("")) {
			logger.error("Field '{}' is not specified correctly in the configuration", field);
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
			logger.error("Field '{}' is not specified correctly in the configuration", field);
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
			logger.error("Field '{}' is not specified correctly in the configuration", field);
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
			logger.error("Field '{}' is not specified correctly in the configuration", field);
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
