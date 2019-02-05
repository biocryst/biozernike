package org.rcsb.biozernike.volume;

public class ResidueVolume {
	public double[] volume;
	public int size;

	ResidueVolume(double[] volume, int size) {
		this.volume = volume;
		this.size = size;
	}

	public double Get(int x, int y, int z) {
		return volume[(z * size + y) * size + x];
	}
}
