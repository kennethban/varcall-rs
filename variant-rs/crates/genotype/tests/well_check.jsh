// jshell --class-path gatk.jar well_check.jsh > well19937c.txt : reference values for src/well.rs
import org.apache.commons.math3.random.*;
var w = new Well19937c(47382911L);
for (int i = 0; i < 5; i++) System.out.print(w.nextInt(100) + " ");
var p = new RandomDataGenerator(new Well19937c(47382911L)).nextPermutation(10, 4);
for (int x : p) System.out.print(x + " ");
var w2 = new Well19937c(47382911L);
for (int i = 0; i < 3; i++) System.out.print(w2.nextInt(1000000007) + " ");
System.out.println();
/exit
