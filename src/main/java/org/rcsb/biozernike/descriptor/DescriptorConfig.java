package org.rcsb.biozernike.descriptor;

public class DescriptorConfig {
	public int[] searchIndicesZernike;
	public double[] searchCoefficientsZernike;
	public double[] searchCoefficientsGeometry;
	public double thresholdZernike;
	public double thresholdGeometry;

	public int maxOrderInvariants;
	public int maxOrderAlignment;

	public int[] normOrders;

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

}
