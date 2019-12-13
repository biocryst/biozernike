package org.rcsb.biozernike.descriptor;

//TODO: split ALIGN into ALIGN_CACHED and ALIGN_RT
public enum DescriptorMode {
	CALCULATE_RAW,
	CALCULATE_PROCESSED,
	COMPARE,
	ALIGN_CACHED,
	COMPARE_ALIGN
}
