package org.rcsb.biozernike.volume;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;

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
	public static void write(Volume volume, File file, MapFileType fileType) throws IOException {

		DataOutputStream os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file), 10485760));

		int[] dims = volume.getDimensions();
		float gridWidth = (float)volume.getGridWidth();
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
//					int flat_ind = (z*volume.dims[1] + y)*volume.dims[0] + x;
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
	 * @param threshold a threshold value
	 * @param multiplier a multiplier value
	 * @return the volume
	 * @throws IOException if problems reading file
	 */
	public static Volume read(File file, MapFileType fileType, double threshold, double multiplier) throws IOException {
		InputStream is = new BufferedInputStream(new FileInputStream(file), 10485760);
		return read(is, fileType, threshold, multiplier);
	}

	/**
	 * Read volume from given input stream in specified format
	 * @param is the input stream
	 * @param fileType the file format that the input stream uses
	 * @param threshold a lower threshold of density values to read. Values below this threshold will not be stored
	 * @param multiplier a multiplier value to normalise the density values
	 * @return the volume
	 * @throws IOException if problems reading file
	 */
	public static Volume read(InputStream is, MapFileType fileType, double threshold, double multiplier) throws IOException {

		int i, j, x, y, z;
		int nc,nr,ns,mode,ncStart,nrStart,nsStart,nx,ny,nz,mapc,mapr,maps,ispg,nsymbt,lskflg,nlabl;
		float xlength, ylength, zlength, alpha, beta, gamma;
		float aMin, aMax, aMean, aRms;
		float[][] skwmat = new float[3][3];
		float[] skwtrn = new float[3];
		float[] origin = new float[3];

		byte[][] label = new byte[10][81];

		DataInputStream dis = new DataInputStream(is);

		/* Read Headers */
		nc = Integer.reverseBytes(dis.readInt());
		nr = Integer.reverseBytes(dis.readInt());
		ns = Integer.reverseBytes(dis.readInt());
		mode = Integer.reverseBytes(dis.readInt());

		ncStart = Integer.reverseBytes(dis.readInt());
		nrStart = Integer.reverseBytes(dis.readInt());
		nsStart = Integer.reverseBytes(dis.readInt());
		nx = Integer.reverseBytes(dis.readInt());
		ny = Integer.reverseBytes(dis.readInt());
		nz = Integer.reverseBytes(dis.readInt());

		xlength = reversedFloat(dis);
		ylength = reversedFloat(dis);
		zlength = reversedFloat(dis);
		alpha = reversedFloat(dis);
		beta = reversedFloat(dis);
		gamma = reversedFloat(dis);

		mapc = Integer.reverseBytes(dis.readInt());
		mapr = Integer.reverseBytes(dis.readInt());
		maps = Integer.reverseBytes(dis.readInt());
		aMin = reversedFloat(dis);
		aMax = reversedFloat(dis);
		aMean = reversedFloat(dis);
		ispg = Integer.reverseBytes(dis.readInt());
		nsymbt = Integer.reverseBytes(dis.readInt());

		/* CCP4 */
		if (fileType == MapFileType.CCP4) {
			lskflg = Integer.reverseBytes(dis.readInt());

			for (i=0;i<3;++i) {
				for (j=0;j<3;++j) {
					skwmat[i][j] = reversedFloat(dis);
				}
			}

			for (i=0;i<3;++i) {
				skwtrn[i] = reversedFloat(dis);
			}
			for (i=0;i<15;++i) {
				dis.readInt();
			}


			origin[0] = origin[1] = origin[2] = (float)0.0;
		}

		/* MRC */
		else if (fileType == MapFileType.MRC){
			for (i=0;i<25;++i){ dis.readInt();}
			origin[0] =  reversedFloat(dis);
			origin[1] =  reversedFloat(dis);
			origin[2] =  reversedFloat(dis);
		}
		dis.readInt();
		dis.readInt();

		aRms = reversedFloat(dis);
		nlabl = Integer.reverseBytes(dis.readInt());

		for (i=0;i<10;++i){
			for (j=0;j<80;++j) {
				label[i][j] = dis.readByte();
			}
		}

		int[] dims = {0, 0, 0};

		dims[0] = nc;   dims[1] = nr;   dims[2] = ns;


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


		/* Read NSYMBT characters */
		for (i=0;i<nsymbt;++i) {
			dis.readByte();
		}


		int dim = dims[0];
		if (dims[1]>dim) {dim = dims[1];}
		if (dims[2]>dim) {dim = dims[2];}

		int flat_dim = dim*dim*dim;

		double[] voxels  = new double[flat_dim];
		double sumval = 0;
		int n_voxels = 0;

		/* Read Voxel */
		for (z=0;z<dims[2];++z){
			for (y=0;y<dims[1];++y){
				for (x=0;x<dims[0];++x){
					double val = 0;

					if (mode==0){ val =  (float)dis.readChar();}
					else if (mode==1){ val =  (float)dis.readShort();}
					else if (mode==2){ val =  Math.abs(reversedFloat(dis));}

					if(val <= threshold) {
						continue;
					}
					n_voxels++;
					voxels[(z*dim + y)*dim + x] = val;
					sumval+=val;
				}
			}
		}

		is.close();

		double norm_coef = n_voxels*multiplier/sumval;

		for (i=0;i<flat_dim;i++) {
			voxels[i]*= norm_coef;
		}

		Volume volume = new Volume();
		volume.createFromData(dims, voxels, gridWidth);

		return volume;
	}

	private static float reversedFloat(DataInputStream is) throws IOException {
		return Float.intBitsToFloat(Integer.reverseBytes(Float.floatToIntBits(is.readFloat())));
	}

}
