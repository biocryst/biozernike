package org.rcsb.biozernike.volume;

import java.util.HashMap;
import java.util.Map;

//TODO: config file
public class VolumeConstants {
	static final Map<String, Double> residueRadius = new HashMap<>();
	static final Map<String, Double> residueWeight = new HashMap<>();

	static {
		double coef = Math.sqrt(5.0 / 3.0);

		residueRadius.put("ALA", 1.963913 * coef);
		residueRadius.put("ARG", 3.374007 * coef);
		residueRadius.put("ASN", 2.695111 * coef);
		residueRadius.put("ASP", 2.525241 * coef);
		residueRadius.put("CYS", 2.413249 * coef);
		residueRadius.put("GLN", 3.088783 * coef);
		residueRadius.put("GLU", 2.883527 * coef);
		residueRadius.put("GLY", 1.841949 * coef);
		residueRadius.put("HIS", 2.652737 * coef);
		residueRadius.put("ILE", 2.575828 * coef);
		residueRadius.put("LEU", 2.736953 * coef);
		residueRadius.put("LYS", 3.177825 * coef);
		residueRadius.put("MET", 2.959014 * coef);
		residueRadius.put("MSE", 2.959014 * coef);
		residueRadius.put("PHE", 2.979213 * coef);
		residueRadius.put("PRO", 2.266054 * coef);
		residueRadius.put("SER", 2.184637 * coef);
		residueRadius.put("THR", 2.366486 * coef);
		residueRadius.put("TRP", 3.248871 * coef);
		residueRadius.put("TYR", 3.217711 * coef);
		residueRadius.put("VAL", 2.351359 * coef);


		residueRadius.put("H", 1.2);
		residueRadius.put("C", 1.7);
		residueRadius.put("N", 1.55);
		residueRadius.put("O", 1.52);
		residueRadius.put("F", 1.47);
		residueRadius.put("P", 1.8);
		residueRadius.put("CL", 1.75);
		residueRadius.put("CU", 1.4);

//		residueRadius.put("A", 4.333750 * coef);
//		residueRadius.put("T", 3.700942 * coef);
//		residueRadius.put("G", 4.443546 * coef);
//		residueRadius.put("C", 3.954067 * coef);
//		residueRadius.put("U", 3.964129 * coef);
//		residueRadius.put("I", 4.0 * coef);
//		residueRadius.put("DA", 4.333750 * coef);
//		residueRadius.put("DT", 3.700942 * coef);
//		residueRadius.put("DG", 4.443546 * coef);
//		residueRadius.put("DC", 3.954067 * coef);
//		residueRadius.put("DU", 3.964129 * coef);
//		residueRadius.put("DI", 4.0 * coef);

		residueWeight.put("ALA", 82.03854);
		residueWeight.put("ARG", 160.09176);
		residueWeight.put("ASN", 148.07768);
		residueWeight.put("ASP", 126.04834);
		residueWeight.put("CYS", 114.10454);
		residueWeight.put("GLN", 160.08868);
		residueWeight.put("GLU", 138.05934);
		residueWeight.put("GLY", 70.02754);
		residueWeight.put("HIS", 146.08502);
		residueWeight.put("ILE", 214.15954);
		residueWeight.put("LEU", 190.13754);
		residueWeight.put("LYS", 132.07828);
		residueWeight.put("MET", 138.12654);
		residueWeight.put("MSE", 138.12654);
		residueWeight.put("PHE", 154.10454);
		residueWeight.put("PRO", 106.06054);
		residueWeight.put("SER", 98.03794);
		residueWeight.put("THR", 146.08194);
		residueWeight.put("TRP", 192.13328);
		residueWeight.put("TYR", 170.10394);
		residueWeight.put("VAL", 178.12654);


		residueWeight.put("H", 1.0);
		residueWeight.put("C", 12.0);
		residueWeight.put("N", 14.0);
		residueWeight.put("O", 16.0);
		residueWeight.put("F", 19.0);
		residueWeight.put("P", 31.0);
		residueWeight.put("CL", 35.45);
		residueWeight.put("CU", 63.55);

//		residueWeight.put("A", 409.12186);
//		residueWeight.put("T", 379.11264);
//		residueWeight.put("G", 425.12126);
//		residueWeight.put("C", 385.09678);
//		residueWeight.put("U", 387.08944);
//		residueWeight.put("I", 400.0);
//		residueWeight.put("DA", 409.12186);
//		residueWeight.put("DT", 379.11264);
//		residueWeight.put("DG", 425.12126);
//		residueWeight.put("DC", 385.09678);
//		residueWeight.put("DU", 387.08944);
//		residueWeight.put("DI", 400.0);
	}

	public static double getRadius(String resName) {
		if (!residueRadius.containsKey(resName)) {
			resName = "ASP"; // a somewhat average residue
		}
		return residueRadius.get(resName);
	}

	public static double getWeight(String resName) {
		if (!residueWeight.containsKey(resName)) {
			resName = "ASP";
		}
		return residueWeight.get(resName);
	}

}

