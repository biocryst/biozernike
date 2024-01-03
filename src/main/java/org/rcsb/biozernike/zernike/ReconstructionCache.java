package org.rcsb.biozernike.zernike;
import org.rcsb.biozernike.complex.Complex;

import java.io.*;
import java.util.HashMap;
import java.util.Map;

public class ReconstructionCache {
    private static Map<Integer, Map<Integer, Map<Integer, Complex[]>>> zp_arr = null;

    public  static Map<Integer, Map<Integer, Map<Integer, Complex[]>>> readZpCache() throws IOException, ClassNotFoundException {
        if (zp_arr==null) {
            System.out.println("Reading cache...");

            ObjectInputStream ois = new ObjectInputStream(new FileInputStream("D:\\PT\\SRC\\biozernike\\test.coefs"));

            Map<Integer, Map<Integer, Map<Integer, Complex[]>>> zp_vals =
                    (Map<Integer, Map<Integer, Map<Integer, Complex[]>>>) ois.readObject();
            ois.close();

            zp_arr = zp_vals;
            System.out.println("Finished reading cache");
        }
        return zp_arr;
    }

    public static void writeZpCache(Map<Integer, Map<Integer, Map<Integer, Complex[]>>> zp_vals) throws Exception {

//        ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream("test.coefs"));
        System.out.println("Writing cache...");
        ObjectOutputStream oos = new  ObjectOutputStream(new FileOutputStream("D:\\PT\\SRC\\biozernike\\test.coefs"));
        oos.writeObject(zp_vals);
        oos.close();
        System.out.println("Done writing cache...");
    }

////        DataInputStream data_in = new DataInputStream(
////                new BufferedInputStream(
////                        new FileInputStream(new File("binary_file.dat"))));
//
////        FILE *fp;
////        fp = fopen("test.coefs","wb");
//
//        int v=zp_vals.size();
//        bos.write(v);
////        fwrite(&v,sizeof(int),1,fp); // number of n-s
//
//        for (Map.Entry<Integer, Map<Integer, Map<Integer, Complex[]>>> entry : zp_vals.entrySet())
//        {
////        for(Integer in: zp_vals) {
//
//            v=entry.getKey(); // n ind
////            fwrite(&v,sizeof(int),1,fp);
//            bos.write(v);
//
//            v=entry.getValue().size(); // number of l-s
//            bos.write(v);
////            fwrite(&v,sizeof(int),1,fp);
//
////            std::cout << "N: "<< in.first<<", len(L): " << in.second.size() <<"\n";
//            for (Map.Entry<Integer, Map<Integer, Complex[]>> il : entry.getValue().entrySet()) {
////            for(auto const& il: in.second) {
////                v=il.first; // l ind
////                fwrite(&v,sizeof(int),1,fp);
//                v=il.getKey(); // n ind
//                bos.write(v);
//
////                v=il.second.size(); // number of m-s
////                fwrite(&v,sizeof(int),1,fp);
//                v=il.getValue().size(); // number of l-s
//                bos.write(v);
//
////                std::cout << "L: "<< il.first<<", len(M): " << il.second.size() <<"\n";
//                for (Map.Entry<Integer, Complex[]> im : il.getValue().entrySet()) {
////                for(auto const& im: il.second) {
//                    v=im.getKey(); // m ind
//                    bos.write(v);
////                    fwrite(&v,sizeof(int),1,fp);
////                    std::cout << "M: "<< il.first<< ", Element:"<<sizeof(ComplexT)<<"\n";
//                    // im.second; //T* dimX*dimY*dimZ
////                    fwrite(im.second,sizeof(ComplexT),volume_size,fp);
//                    bos.write(im.getValue(),);
//                }
//            }
//        }
//
//        fclose(fp);

//    }

//    public void readZpCache (ZCoefsMapT & zp, String filename, int dim){
//        int volume_size = dim * dim * dim;
//
//        FILE * fp;
//        fp = fopen(filename, "rb");
//
//        int n_N;
//
//        fread( & n_N, sizeof( int),1, fp);  // number of n-s
//
//        for (int in = 0; in < n_N; in++) {
//            int n_cur;
//            fread( & n_cur, sizeof( int),1, fp);  // n ind
//
//            int n_L;
//            fread( & n_L, sizeof( int),1, fp);  // number of l-s
//
//            for (int il = 0; il < n_L; il++) {
//                int l_cur;
//                fread( & l_cur, sizeof( int),1, fp);  // l ind
//
//                int n_M;
//                fread( & n_M, sizeof( int),1, fp);  // number of m-s
//
//                for (int im = 0; im < n_M; im++) {
//                    int m_cur;
//                    fread( & m_cur, sizeof( int),1, fp);  // m ind
//                    auto * zp_xyz = new ComplexT[volume_size] ();
//                    fread(zp_xyz, sizeof(ComplexT), volume_size, fp);  // volume
//                    zp[n_cur][l_cur][m_cur] = zp_xyz;
//                }
//            }
//        }
//
//        fclose(fp);
//    }
}
