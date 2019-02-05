package org.rcsb.biozernike.volume;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;

// very simple volume read/write in CCP4 format, essentially copied from gmconvert tool.
// Extremely limited file format support, intended for debugging purposes only.
public class VolumeIO {

	public void write(Volume volume, String filename) throws IOException {

		DataOutputStream os = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(filename), 10485760));

		int[] dims = volume.getDimensions();
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
		os.writeFloat((float) dims[0]);
		os.writeFloat((float) dims[1]);
		os.writeFloat((float) dims[2]);

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
}
