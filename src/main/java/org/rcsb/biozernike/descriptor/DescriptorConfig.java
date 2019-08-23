package org.rcsb.biozernike.descriptor;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class DescriptorConfig {
	private static final Logger logger = LoggerFactory.getLogger(DescriptorConfig.class);

	public int[] searchIndicesZernike;
	public double[] searchCoefficientsZernike;
	public double[] searchCoefficientsGeometry;
	public double thresholdZernike;
	public double thresholdGeometry;

	public int maxOrderInvariants;
	public int maxOrderAlignment;
	public int[] normOrders;

	private String propsFile = "";
	private DescriptorConfig(){};

	public DescriptorConfig(
			int maxOrderInvariants,
			int maxOrderAlignment,
			int[] normOrders,
			int[] searchIndicesZernike,
	        double[] searchCoefficientsZernike,
	        double[] searchCoefficientsGeometry,
	        double thresholdZernike,
	        double thresholdGeometry) {
		this.maxOrderInvariants = maxOrderInvariants;
		this.maxOrderAlignment = maxOrderAlignment;
		this.normOrders = normOrders;
		this.searchIndicesZernike = searchIndicesZernike;
		this.searchCoefficientsZernike = searchCoefficientsZernike;
		this.searchCoefficientsGeometry = searchCoefficientsGeometry;
		this.thresholdZernike = thresholdZernike;
		this.thresholdGeometry = thresholdGeometry;
	}

	public DescriptorConfig(String propsFile) {
		this.propsFile = propsFile;

		Properties props = new Properties();
		try {
			InputStream input = new FileInputStream(propsFile);
			props.load(input);
		} catch (IOException e) {
			logger.error("Invalid or missing config file {}", propsFile);
			throw new RuntimeException("Missing configuration '"+propsFile+"'. Can't continue.");
		}

		thresholdZernike =  loadDoubleField(props,"threshold.zernike");
		thresholdGeometry = loadDoubleField(props,"threshold.geometry");
		searchIndicesZernike = loadIntArrayField(props, "search.indices.zernike");
		searchCoefficientsZernike = loadDoubleArrayField(props, "search.coefficients.zernike");
		searchCoefficientsGeometry = loadDoubleArrayField(props, "search.coefficients.geometry");
		normOrders = loadIntArrayField(props, "search.norm.orders.zernike");
		maxOrderInvariants = loadIntegerField(props,"search.max.order.zernike");
		maxOrderAlignment = loadIntegerField(props,"align.max.order.zernike");
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
