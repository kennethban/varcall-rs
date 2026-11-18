// Emits java.lang.Math.{exp,log,log10,pow,log1p,expm1} results for many inputs, as binary records on stdout:
//   u8 fn (0 exp, 1 log, 2 log10, 3 pow, 4 log1p, 5 expm1), u8 category, u64 x, u64 y, u64 result  (little endian, raw bits)
// Usage: java JMathGen.java [millions per function, default 10] [seed] | jmath_check
import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.SplittableRandom;

public class JMathGen {
    static OutputStream out;
    static final byte[] rec = new byte[26];

    static void put(int fn, int cat, double x, double y, double r) throws IOException {
        rec[0] = (byte) fn;
        rec[1] = (byte) cat;
        le(2, Double.doubleToRawLongBits(x));
        le(10, Double.doubleToRawLongBits(y));
        le(18, Double.doubleToRawLongBits(r));
        out.write(rec);
    }

    static void le(int o, long v) {
        for (int i = 0; i < 8; i++) rec[o + i] = (byte) (v >>> (8 * i));
    }

    static SplittableRandom rnd;

    static double anyBits() { return Double.longBitsToDouble(rnd.nextLong()); }
    static double posFinite() {  // uniform over positive finite bit patterns (incl. subnormals)
        while (true) {
            double d = Double.longBitsToDouble(rnd.nextLong() & 0x7fffffffffffffffL);
            if (d < Double.POSITIVE_INFINITY) return d;
        }
    }
    static double subnormal() { return Double.longBitsToDouble(1 + rnd.nextLong(0x000fffffffffffffL)); }
    static double uni(double a, double b) { return a + (b - a) * rnd.nextDouble(); }
    static double nearOne() {
        int k = 1 + rnd.nextInt(60);
        return 1.0 + Math.scalb(uni(-1, 1), -k);
    }

    static final double[] SPECIAL = {
        0.0, -0.0, 1.0, -1.0, 2.0, -2.0, 0.5, -0.5, 3.0, -3.0, 10.0, -10.0, 0.1, 1.5, -1.5,
        Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NaN,
        Double.longBitsToDouble(0x7ff0000000000001L), Double.longBitsToDouble(0xfff8000000000001L),
        Double.MIN_VALUE, -Double.MIN_VALUE, Double.MIN_NORMAL, -Double.MIN_NORMAL, Double.MAX_VALUE, -Double.MAX_VALUE,
        Math.nextUp(1.0), Math.nextDown(1.0), Math.nextUp(-1.0), Math.nextDown(-1.0),
        709.782712893384, 709.7827128933841, -745.1332191019411, -745.1332191019412, -708.3964185322641, -708.3964185322642,
        1023.0, 1024.0, -1022.0, -1074.0, -1075.0, 1e300, -1e300, 1e-300, 4503599627370496.0, 9007199254740992.0,
        9007199254740993.0, -9007199254740991.0, 4503599627370497.0, 4503599627370495.5, 1e17, -1e17, 1e-17, 2.5, -2.5,
    };

    public static void main(String[] a) throws IOException {
        long n = (a.length > 0 ? Long.parseLong(a[0]) : 10) * 1_000_000L;
        long seed = a.length > 1 ? Long.parseLong(a[1]) : 20261116L;
        rnd = new SplittableRandom(seed);
        out = new BufferedOutputStream(System.out, 1 << 20);

        // ---- exp ----
        for (double s : SPECIAL) put(0, 0, s, 0, Math.exp(s));
        for (int i = -1600; i <= 1600; i++) { double x = i * 0.5; put(0, 0, x, 0, Math.exp(x)); }
        for (long i = 0; i < n; i++) {
            int cat = (int) (i % 8) + 1;
            double x;
            switch (cat) {
                case 1: x = anyBits(); break;
                case 2: x = uni(-750, 750); break;
                case 3: x = uni(-1, 1); break;
                case 4: x = Math.scalb(uni(-1, 1), -rnd.nextInt(1080)); break;  // tiny, down to subnormal
                case 5: x = uni(-746, -700); break;                              // subnormal / underflow results
                case 6: x = uni(700, 710); break;                                // overflow edge
                case 7: x = uni(-50, 0); break;                                  // log-probability range
                default: x = uni(-30, 30); break;
            }
            put(0, cat, x, 0, Math.exp(x));
        }
        // ---- log, log10 ----
        for (int fn = 1; fn <= 2; fn++) {
            for (double s : SPECIAL) put(fn, 0, s, 0, fn == 1 ? Math.log(s) : Math.log10(s));
            for (int i = 1; i <= 100000; i++) put(fn, 0, i, 0, fn == 1 ? Math.log(i) : Math.log10(i));
            for (int e = -330; e <= 310; e++) { double x = Double.parseDouble("1e" + e); put(fn, 0, x, 0, fn == 1 ? Math.log(x) : Math.log10(x)); }
            for (long i = 0; i < n; i++) {
                int cat = (int) (i % 8) + 1;
                double x;
                switch (cat) {
                    case 1: x = anyBits(); break;
                    case 2: x = posFinite(); break;
                    case 3: x = nearOne(); break;
                    case 4: x = subnormal(); break;
                    case 5: x = rnd.nextDouble(); break;                           // ratios in [0,1)
                    case 6: x = uni(0.5, 2); break;
                    case 7: x = rnd.nextInt(1, 1 << 30) / (double) rnd.nextInt(1, 1 << 30); break;
                    default: x = Math.scalb(1.0 + rnd.nextDouble(), rnd.nextInt(-1100, 1030)); break;
                }
                put(fn, cat, x, 0, fn == 1 ? Math.log(x) : Math.log10(x));
            }
        }
        // ---- pow ----
        for (double s : SPECIAL) for (double t : SPECIAL) put(3, 0, s, t, Math.pow(s, t));
        for (int q = 0; q <= 2000; q++) { double y = q / -10.0; put(3, 0, 10.0, y, Math.pow(10.0, y)); }
        for (long i = 0; i < n; i++) {
            int cat = (int) (i % 12) + 1;
            double x, y;
            switch (cat) {
                case 1: x = 10.0; y = uni(-30, 0); break;                           // BQSR qualToErrorProb
                case 2: x = 10.0; y = Math.round(uni(0, 600) * 10000) / 10000.0 / -10.0; break;  // reported quals (4 dp)
                case 3: x = 10.0; y = uni(-330, 310); break;
                case 4: x = posFinite(); y = uni(-10, 10); break;
                case 5: x = uni(0, 2); y = uni(-1000, 1000); break;
                case 6: x = -posFinite(); y = rnd.nextInt(-100, 101); break;          // negative x, integer y
                case 7: x = uni(-4, 4); y = rnd.nextInt(-2000, 2001) * 0.5; break;     // integer and half-integer y
                case 8: x = posFinite(); y = rnd.nextInt(-64, 65); break;
                case 9: x = -uni(0, 2); y = Math.scalb((double) rnd.nextLong(1L << 53), rnd.nextInt(-3, 12)); break; // huge (odd/even) y
                case 10: x = anyBits(); y = anyBits(); break;
                case 11: x = nearOne(); y = uni(-1e6, 1e6); break;
                default: x = subnormal() * (rnd.nextBoolean() ? 1 : -1); y = uni(-3, 3); break;
            }
            put(3, cat, x, y, Math.pow(x, y));
        }
        // ---- log1p (StrictMath / fdlibm) ----
        final double LN_HALF = Math.log(0.5);
        for (double s : SPECIAL) put(4, 0, s, 0, Math.log1p(s));
        for (int i = -1000; i <= 1000; i++) { double x = i / 1000.0; put(4, 0, x, 0, Math.log1p(x)); }
        for (long i = 0; i < n; i++) {
            int cat = (int) (i % 10) + 1;
            double x;
            switch (cat) {
                case 1: x = anyBits(); break;
                case 2: x = posFinite(); break;
                case 3: x = uni(-1, 0); break;                                    // (-1, 0]
                case 4: { double la = uni(-750, LN_HALF); x = -Math.exp(la); break; } // log1mexp: log1p(-exp(a)), a < ln 0.5
                case 5: x = Math.scalb(uni(-1, 1), -rnd.nextInt(1, 1080)); break;  // tiny, down to subnormal
                case 6: x = uni(-0.2929, 0.41422); break;                          // k = 0 range
                case 7: x = uni(0.41, 10); break;
                case 8: x = -1.0 + Math.scalb(rnd.nextDouble(), -rnd.nextInt(1, 60)); break; // just above -1
                case 9: x = subnormal() * (rnd.nextBoolean() ? 1 : -1); break;
                default: x = Math.scalb(1.0 + rnd.nextDouble(), rnd.nextInt(-60, 1024)) * (rnd.nextInt(4) == 0 ? -1 : 1); break;
            }
            put(4, cat, x, 0, Math.log1p(x));
        }
        // ---- expm1 (StrictMath / fdlibm) ----
        for (double s : SPECIAL) put(5, 0, s, 0, Math.expm1(s));
        for (int i = -1600; i <= 1600; i++) { double x = i * 0.5; put(5, 0, x, 0, Math.expm1(x)); }
        for (long i = 0; i < n; i++) {
            int cat = (int) (i % 10) + 1;
            double x;
            switch (cat) {
                case 1: x = anyBits(); break;
                case 2: x = uni(LN_HALF, 0); break;                               // log1mexp: log(-expm1(a)), ln 0.5 <= a < 0
                case 3: x = uni(-750, 750); break;
                case 4: x = uni(-1.1, 1.1); break;                                 // k = 0, +-1 ranges
                case 5: x = Math.scalb(uni(-1, 1), -rnd.nextInt(1, 1080)); break; // tiny, down to subnormal
                case 6: x = uni(-40, 40); break;                                   // k up to 56
                case 7: x = uni(700, 710); break;                                  // overflow edge
                case 8: x = uni(-45, -35); break;                                  // -56 ln2 edge
                case 9: x = subnormal() * (rnd.nextBoolean() ? 1 : -1); break;
                default: x = uni(-30, 30); break;
            }
            put(5, cat, x, 0, Math.expm1(x));
        }
        out.flush();
    }
}
