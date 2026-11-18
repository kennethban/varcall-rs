import java.io.*;
import java.nio.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.CombinatoricsUtils;

/** Dumps the exact GATK numeric tables used by HaplotypeCaller's isActive() so the Rust port is bit-compatible. */
public class GenTables {
    public static void main(String[] a) throws Exception {
        // genotype priors used by AlleleFrequencyCalculator.calculateSingleSampleBiallelicNonRefPosterior (ploidy 2)
        final double refPseudo = 1e-3 / Math.pow(0.01, 2);
        final double snpPseudo = 1e-3 * refPseudo;
        final double[] prior = new double[3];
        for (int n = 0; n <= 2; n++) {
            prior[n] = MathUtils.logToLog10(CombinatoricsUtils.binomialCoefficientLog(2, n) + Gamma.logGamma(n + snpPseudo) + Gamma.logGamma(2 - n + refPseudo));
        }
        // band-pass kernel (BandPassActivityProfile, adaptive filter size, sigma 17, max 50)
        final double sigma = 17.0;
        final double[] full = makeKernel(50, sigma);
        final int middle = (full.length - 1) / 2;
        int filterEnd = middle;
        while (filterEnd > 0) { if (full[filterEnd - 1] < 1e-5) break; filterEnd--; }
        final int filterSize = middle - filterEnd;
        final double[] kernel = makeKernel(filterSize, sigma);
        // jacobian table, exactly as MathUtils.JacobianLogTable
        final int jsize = (int) (8.0 / 0.0001) + 1;
        final ByteBuffer bb = ByteBuffer.allocate(8 * (3 + 1 + kernel.length + 256 + 1 + jsize + 4)).order(ByteOrder.LITTLE_ENDIAN);
        for (double p : prior) bb.putDouble(p);
        bb.putDouble(filterSize);
        for (double k : kernel) bb.putDouble(k);
        for (int q = 0; q < 256; q++) bb.putDouble(q < 255 ? QualityUtils.qualToProbLog10((byte) q) : Double.NaN);  // MAX_QUAL=254
        bb.putDouble(jsize);
        for (int k = 0; k < jsize; k++) bb.putDouble(Math.log10(1.0 + Math.pow(10.0, -k * 0.0001)));
        bb.putDouble(MathUtils.LOG10_ONE_THIRD);
        bb.putDouble(Math.log10(2));
        bb.putDouble(Math.log10(1));
        bb.putDouble(1.0 / 0.0001);
        try (FileOutputStream f = new FileOutputStream(a[0])) { f.write(bb.array()); }
        System.out.printf("prior=%s filterSize=%d kernel[mid]=%.17g%n", java.util.Arrays.toString(prior), filterSize, kernel[filterSize]);
    }
    static double[] makeKernel(int filterSize, double sigma) {
        final double[] k = new double[2 * filterSize + 1];
        for (int i = 0; i < k.length; i++) k[i] = MathUtils.normalDistribution(filterSize, sigma, i);
        return MathUtils.normalizeSumToOne(k);
    }
}
