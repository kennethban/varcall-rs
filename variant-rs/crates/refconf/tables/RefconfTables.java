import java.io.*;
import java.nio.*;
import java.lang.reflect.*;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceModel;

/** Dumps the exact GATK numbers used by ReferenceConfidenceModel (GVCF mode) for the Rust port:
 *  qualToProbLog10[0..255], LOG10_ONE_THIRD, log10(2), and the diploid indel GL cache (n = 0..40). */
public class RefconfTables {
    public static void main(String[] a) throws Exception {
        final ReferenceConfidenceModel m = new ReferenceConfidenceModel(new IndexedSampleList("S"), new SAMFileHeader(), 10, 0, (byte) 30, true, false);
        final Method gi = ReferenceConfidenceModel.class.getDeclaredMethod("getIndelPLs", int.class, int.class);
        gi.setAccessible(true);
        final ByteBuffer bb = ByteBuffer.allocate(8 * (256 + 2 + 41 * 3)).order(ByteOrder.LITTLE_ENDIAN);
        for (int q = 0; q < 256; q++) bb.putDouble(q < 255 ? QualityUtils.qualToProbLog10((byte) q) : Double.NaN);
        bb.putDouble(MathUtils.LOG10_ONE_THIRD);
        bb.putDouble(Math.log10(2));
        for (int n = 0; n <= 40; n++) {
            final GenotypeLikelihoods gl = (GenotypeLikelihoods) gi.invoke(m, 2, n);
            for (double d : gl.getAsVector()) bb.putDouble(d);
            if (n < 3 || n == 40) System.out.println(n + " " + java.util.Arrays.toString(gl.getAsVector()) + " " + java.util.Arrays.toString(gl.getAsPLs()));
        }
        try (FileOutputStream f = new FileOutputStream(a[0])) { f.write(bb.array()); }
    }
}
