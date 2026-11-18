import org.broadinstitute.hellbender.transformers.BQSRReadTransformer;
int[][] sets = {{10,20,30},{10,20,30,40,50},{20},{3,10,25},{2,4},{30,10,20},{7,8,9,60}};
for (int[] s : sets) for (boolean rd : new boolean[]{false,true}) {
  java.util.List<Integer> l = new java.util.ArrayList<>(); for (int x : s) l.add(x);
  byte[] m = BQSRReadTransformer.constructStaticQuantizedMapping(l, rd);
  StringBuilder sb = new StringBuilder(); for (byte b : m) sb.append(b & 0xff).append(',');
  System.out.println(java.util.Arrays.toString(s).replace(" ","") + " " + rd + " " + sb);
}
/exit
