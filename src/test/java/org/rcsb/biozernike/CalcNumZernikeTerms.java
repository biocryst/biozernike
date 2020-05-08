package org.rcsb.biozernike;

/**
 * Utility to calculate the number of 3D Zernike terms for different max orders (N).
 */
public class CalcNumZernikeTerms {

    public static void main(String[] args) {
        System.out.println("Number of 3D Zernike terms for different max orders (N)");
        System.out.printf("%5s %5s %5s %5s \n", "N", "all", "m>=0", "3DZD");
        for (int i=0; i<=20; i++) {
            int[] terms = calcNumTerms(i);
            System.out.printf("%5d %5d %5d %5d\n", i, terms[0], terms[1], terms[2]);
        }
    }

    private static int[] calcNumTerms(int maxOrder) {
        int count = 0;
        int countL = 0;
        int countMpositive = 0;
        for (int n = 0; n <= maxOrder; n++) {
            for (int l = 0; l <= n; l++) {
                if ((n-l)%2 != 0) continue;
                //System.out.print(" " + n + "," + l);
                countL++;
                for (int m = -l; m <= l; m++) {
                    count++;
                    if (m >=0 ) countMpositive++;
                }
            }
            //System.out.println();
        }

        return new int[]{count, countMpositive, countL};
    }
}
