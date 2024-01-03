package org.rcsb.biozernike.volume;

import javax.media.j3d.BoundingBox;
import javax.vecmath.Point3d;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

/**
 * Very simple volume read/write in CCP4 format, essentially copied from gmconvert tool.
 * Extremely limited file format support, intended for debugging purposes only.
 * <p>
 * See some documentation here: https://ftp.wwpdb.org/pub/emdb/doc/Map-format/current/EMDB_map_format.pdf
 * @author Dmytro Guzenko
 * @author Jose Duarte
 */
public class VolumeIO {

	/**
	 * Write volume to given file in specified format
	 * @param volume the volume
	 * @param file the output file
	 * @param fileType the file type
	 * @throws IOException if problems writing file out
	 */
	public static void write(Volume volume, String filename, MapFileType fileType) throws IOException {

		DataOutputStream os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(filename), 10485760));
		int[] dims = volume.getDimensions();
		float gridWidth = (float)volume.getGridWidth();
		BoundingBox bb = volume.getBoundingBox();
		Point3d pLower = new Point3d();
		bb.getLower(pLower);

//		float grid_width = 1;
		int i, x, y, z;

		/* NC, NR, NS */
		os.writeInt(dims[0]);
		os.writeInt(dims[1]);
		os.writeInt(dims[2]);

		/* Mode */
		i = 2;
		os.writeInt(i);

		/* NCSTART, NRSTART, NSSTART */
		//TODO: write the actual corner
//		os.writeInt((int)pLower.x);
//		os.writeInt((int)pLower.y);
//		os.writeInt((int)pLower.z);
		os.writeInt(0);
		os.writeInt(0);
		os.writeInt(0);

		/* NX, NY, NZ */
		os.writeInt(dims[0]);
		os.writeInt(dims[1]);
		os.writeInt(dims[2]);

		/* X length, Y length, Z length */
		os.writeFloat((float) dims[0]*gridWidth);
		os.writeFloat((float) dims[1]*gridWidth);
		os.writeFloat((float) dims[2]*gridWidth);

		/* Alpha, Beta, Gamma */
		os.writeFloat(90);
		os.writeFloat(90);
		os.writeFloat(90);

		/* MAPC, MAPR, MAPS */

		os.writeInt(1);
		os.writeInt(2);
		os.writeInt(3);

		double sum = 0;
		double sq_sum = 0;
		double max_val = 0;
		for (z = 0; z < dims[2]; ++z) {
			for (y = 0; y < dims[1]; ++y) {
				for (x = 0; x < dims[0]; ++x) {
					double val = volume.getValue(x, y, z);
					if (val > max_val) {
						max_val = val;
					}
					sum += val;
					sq_sum += val * val;
				}
			}
		}
		int n = dims[0] * dims[1] * dims[2];
		double mean_val = sum / n;
		double stdev = Math.sqrt(sq_sum / n - mean_val * mean_val);

		/* AMIN, AMAX, AMEAN */
		os.writeFloat(0);
		os.writeFloat((float) max_val);
		os.writeFloat((float) mean_val);

		/* ISPG, NSYMBT */
		os.writeInt(1);
		os.writeInt(0);

		/* for 'CCP4' (*.map) format **/
		if (fileType == MapFileType.CCP4) {
			/* LSKFLG */
			os.writeInt(0);

			/* SKWMAT11, SKWMAT12, ..., SKWMAT33 */
			for (x = 1; x <= 3; ++x) {
				for (y = 1; y <= 3; ++y) {
					os.writeFloat(0);
				}
			}
			/* SKWTRN1, SKWTRN2, SKWTRN2 */
			os.writeFloat(0);
			os.writeFloat(0);
			os.writeFloat(0);
			/* future use (from 38 to 52 words) */
			for (x = 38; x <= 52; ++x) {
				os.writeInt(0);
			}
		}

		/* for 'MRC' (*.mrc) format **/
		else if (fileType == MapFileType.MRC) {
			/* EXTRA */
			for (x = 25; x <= 49; ++x) {
				os.writeInt(0);
			}
			/* ORIGIN */
			os.writeFloat(0);
			os.writeFloat(0);
			os.writeFloat(0);
		}

		/* MAP, MACHST */
		os.writeBytes("MAP ");
		os.writeBytes("DA  ");

		/* ARMS */
		os.writeFloat((float) stdev);

		/* NLABL */
		os.writeInt(10);

		/* LABEL */
		String space80 = new String(new char[80]).replace('\0', ' ');

		for (x = 0; x < 10; ++x) {
			os.writeBytes(space80);
		}

		/* write voxel values */
		for (z = 0; z < dims[2]; ++z) {
			for (y = 0; y < dims[1]; ++y) {
				for (x = 0; x < dims[0]; ++x) {
//					int flat_ind = (z*volume.dims[1] + y)*volume.dims[0] + x;
					os.writeFloat((float) volume.getValue(x, y, z));
				}
			}
		}

		os.close();
	}

	/**
	 * Read volume from given file in specified format
	 * @param file the file
	 * @param fileType the volume file type
	 * @return the volume
	 * @throws IOException if problems reading file
	 */
	public static Volume read(File file, MapFileType fileType,boolean doReverse) throws IOException {
		InputStream is = new BufferedInputStream(new FileInputStream(file), 10485760);
		return read(is, fileType,doReverse);
	}

	/**
	 * Read volume from given input stream in specified format
	 * @param is the input stream
	 * @param fileType the file format that the input stream uses
	 * @return the volume
	 * @throws IOException if problems reading file
	 */
	public static Volume read(InputStream is, MapFileType fileType, boolean doReverse) throws IOException {

		float[][] skwmat = new float[3][3];
		float[] skwtrn = new float[3];
		float[] origin = new float[3];
		byte[][] label = new byte[10][81];

		DataInputStream dis = new DataInputStream(is);

		/* Read Headers */
		int nc = reversedInt(dis,doReverse);
		int nr = reversedInt(dis,doReverse);
		int ns = reversedInt(dis,doReverse);
		int mode = reversedInt(dis,doReverse);
		int ncStart = reversedInt(dis,doReverse);
		int nrStart = reversedInt(dis,doReverse);
		int nsStart = reversedInt(dis,doReverse);
		int nx = reversedInt(dis,doReverse);
		int ny = reversedInt(dis,doReverse);
		int nz = reversedInt(dis,doReverse);

		float xlength = reversedFloat(dis,doReverse);
		float ylength = reversedFloat(dis,doReverse);
		float zlength = reversedFloat(dis,doReverse);
		float alpha = reversedFloat(dis,doReverse);
		float beta = reversedFloat(dis,doReverse);
		float gamma = reversedFloat(dis,doReverse);

		int mapc = reversedInt(dis,doReverse);
		int mapr = reversedInt(dis,doReverse);
		int maps = reversedInt(dis,doReverse);
		float aMin = reversedFloat(dis,doReverse);
		float aMax = reversedFloat(dis,doReverse);
		float aMean = reversedFloat(dis,doReverse);
		int ispg = reversedInt(dis,doReverse);
		int nsymbt = reversedInt(dis,doReverse);

		/* CCP4 */
		if (fileType == MapFileType.CCP4) {
			int lskflg = reversedInt(dis,doReverse);

			for (int i=0; i<3; ++i) {
				for (int j=0; j<3; ++j) {
					skwmat[i][j] = reversedFloat(dis,doReverse);
				}
			}

			for (int i=0; i<3; ++i) {
				skwtrn[i] = reversedFloat(dis,doReverse);
			}
			for (int i=0;i<15;++i) {
				dis.readInt();
			}

			origin[0] = origin[1] = origin[2] = (float)0.0;
		}

		/* MRC */
		else if (fileType == MapFileType.MRC){
			for (int i=0; i<25; ++i) dis.readInt();
			origin[0] =  reversedFloat(dis,doReverse);
			origin[1] =  reversedFloat(dis,doReverse);
			origin[2] =  reversedFloat(dis,doReverse);
		}
		dis.readInt();
		dis.readInt();

		float aRms = reversedFloat(dis,doReverse);
		int nlabl = reversedInt(dis, doReverse);

		for (int i=0; i<10; ++i){
			for (int j=0; j<80; ++j) {
				label[i][j] = dis.readByte();
			}
		}

		/* Read NSYMBT characters */
		for (int i=0;i<nsymbt;++i) {
			dis.readByte();
		}

		double gridWidth = xlength/(float)nx;
		if (gridWidth < 0.0){
			gridWidth = ylength/(float)ny;
		}
		if (gridWidth < 0.0){
			gridWidth = zlength/(float)nz;
		}

		double[] orig_pos = new double[3];
		orig_pos[0] = ncStart * gridWidth + origin[0];
		orig_pos[1] = nrStart * gridWidth + origin[1];
		orig_pos[2] = nsStart * gridWidth + origin[2];

		int[] dims = {nc, nr, ns};
		double[] voxels  = new double[dims[0]*dims[1]*dims[2]];

		/* Read Voxel */
		for (int z=0; z<dims[2]; ++z){
			for (int y=0; y<dims[1]; ++y){
				for (int x=0; x<dims[0]; ++x){
					double val = 0;

					if (mode==0)
						val = (float)dis.readChar();
					else if (mode==1)
						val = (float)dis.readShort();
					else if (mode==2)
						val = reversedFloat(dis,doReverse);

					// [(z * dimensions[1] + y) * dimensions[0] + x];
					voxels[(z*dims[1] + y)*dims[0] + x] = val;
				}
			}
		}

		is.close();

		Volume volume = new Volume();
		volume.createFromData(dims, voxels, gridWidth);

		return volume;
	}

	private static float reversedFloat(DataInputStream is, boolean doReverse) throws IOException {
		if (doReverse) 
			return Float.intBitsToFloat(Integer.reverseBytes(Float.floatToIntBits(is.readFloat())));
		else
			return is.readFloat();
	}

	private static int reversedInt(DataInputStream is, boolean doReverse) throws IOException {
		if (doReverse)
			return Integer.reverseBytes(is.readInt());
		else
			return is.readInt();
	}

}
