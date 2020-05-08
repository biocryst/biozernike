package org.rcsb.biozernike.volume;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

// very simple volume read/write in CCP4 format, essentially copied from gmconvert tool.
// Extremely limited file format support, intended for debugging purposes only.
public class VolumeIO {

	public static void write(Volume volume, String filename) throws IOException {

		DataOutputStream os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(filename), 10485760));

		int[] dims = volume.getDimensions();
		float gridWidth = (float)volume.getGridWidth();
//		float grid_width = 1;
		char FileType = 'C';
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
		if (FileType == 'C') {
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
		else if (FileType == 'M') {
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

	public static Volume read(String filename, double threshold, double multiplier) throws IOException {

		DataInputStream is = new DataInputStream(new BufferedInputStream(new FileInputStream(filename), 10485760));

		int L,i,j,x,y,z,malloc_ok;
		int NC,NR,NS,MODE,NCSTART,NRSTART,NSSTART,NX,NY,NZ,MAPC,MAPR,MAPS,ISPG,NSYMBT,LSKFLG,NLABL;
		float Xlength, Ylength, Zlength,Alpha,Beta,Gamma;
		float AMIN,AMAX,AMEAN,ARMS;
		float[][] SKWMAT = new float[3][3];
		float[] SKWTRN = new float[3];
		float[] ORIGIN = new float[3];

		char[] MAP = new char[5];
		char[] MACHST = new char[5];
		byte[][] LABEL = new byte[10][81];
		char[] command = new char[128];

		char  FileType = 'M'; /* 'C':ccp4 (*.map), 'M':mrc (*.mrc)  */

		NC = Integer.reverseBytes(is.readInt());

		/* printf("#NC %d Order %c\n",NC,Order); */
		/** Read Headers **/
		NR = Integer.reverseBytes(is.readInt());
		NS = Integer.reverseBytes(is.readInt());
		MODE = Integer.reverseBytes(is.readInt());

		NCSTART = Integer.reverseBytes(is.readInt());
		NRSTART = Integer.reverseBytes(is.readInt());
		NSSTART = Integer.reverseBytes(is.readInt());
		NX = Integer.reverseBytes(is.readInt());
		NY = Integer.reverseBytes(is.readInt());
		NZ = Integer.reverseBytes(is.readInt());

		Xlength = reversedFloat(is);
		Ylength = reversedFloat(is);
		Zlength = reversedFloat(is);
		Alpha = reversedFloat(is);
		Beta = reversedFloat(is);
		Gamma = reversedFloat(is);

		MAPC = Integer.reverseBytes(is.readInt());
		MAPR = Integer.reverseBytes(is.readInt());
		MAPS = Integer.reverseBytes(is.readInt());
		AMIN = reversedFloat(is);
		AMAX = reversedFloat(is);
		AMEAN = reversedFloat(is);
		ISPG = Integer.reverseBytes(is.readInt());
		NSYMBT = Integer.reverseBytes(is.readInt());

		/** CCP4 **/
		if (FileType=='C'){
			LSKFLG = Integer.reverseBytes(is.readInt());

			for (i=0;i<3;++i){
				for (j=0;j<3;++j){
					SKWMAT[i][j] = reversedFloat(is);
				}
			}

			for (i=0;i<3;++i){
				SKWTRN[i] = reversedFloat(is);
			}
			for (i=0;i<15;++i){is.readInt();}


			ORIGIN[0] = ORIGIN[1] = ORIGIN[2] = (float)0.0;
		}

		/** MRC **/
		else if (FileType=='M'){
			for (i=0;i<25;++i){ is.readInt();}
			ORIGIN[0] =  reversedFloat(is);
			ORIGIN[1] =  reversedFloat(is);
			ORIGIN[2] =  reversedFloat(is);
		}
		is.readInt();
		is.readInt();

		ARMS = reversedFloat(is);
		NLABL = Integer.reverseBytes(is.readInt());

		for (i=0;i<10;++i){
			for (j=0;j<80;++j) {
				LABEL[i][j] = is.readByte();
			}
		}

		int[] dims_ = {0, 0, 0};

		dims_[0] = NC;   dims_[1] = NR;   dims_[2] = NS;


		double gridWidth = (float)Xlength/(float)NX;
		if (gridWidth < 0.0){
			gridWidth = (float)Ylength/(float)NY;
		}
		if (gridWidth < 0.0){
			gridWidth = (float)Zlength/(float)NZ;
		}

		double[] orig_pos = new double[3];
		orig_pos[0] = NCSTART * gridWidth + ORIGIN[0];
		orig_pos[1] = NRSTART * gridWidth + ORIGIN[1];
		orig_pos[2] = NSSTART * gridWidth + ORIGIN[2];


		/** Read NSYMBT characters **/
		for (i=0;i<NSYMBT;++i){ is.readByte(); }


		int dim_ = dims_[0];
		if (dims_[1]>dim_) {dim_ = dims_[1];}
		if (dims_[2]>dim_) {dim_ = dims_[2];}

		int flat_dim = dim_*dim_*dim_;

		double[] voxels_  = new double[flat_dim];
		double sumval_ = 0;
		int n_voxels = 0;

		/** Read Voxel **/
		for (z=0;z<dims_[2];++z){
			for (y=0;y<dims_[1];++y){
				for (x=0;x<dims_[0];++x){
					double val = 0;

					if (MODE==0){ val =  (float)is.readChar();}
					else if (MODE==1){ val =  (float)is.readShort();}
					else if (MODE==2){ val =  Math.abs(reversedFloat(is));}

					if(val <= threshold) {
						continue;
					}
					n_voxels++;
					voxels_[(z*dim_ + y)*dim_ + x] = val;
					sumval_+=val;
				}
			}
		}

		is.close();

		double norm_coef = n_voxels*multiplier/sumval_;

		for (i=0;i<flat_dim;i++) {
			voxels_[i]*= norm_coef;
		}

		Volume volume = new Volume();
		volume.createFromData(dims_, voxels_, gridWidth);

		return volume;
	}

	private static float reversedFloat(DataInputStream is) throws IOException {
		return Float.intBitsToFloat(Integer.reverseBytes(Float.floatToIntBits(is.readFloat())));
	}

}
