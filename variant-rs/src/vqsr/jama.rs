//! The parts of Jama 1.1.1 (`gov.nist.math.jama`) that GATK's VQSR uses, operation for operation: square
//! matrices, the "left-looking" Crout/Doolittle LU decomposition with partial pivoting, `solve`, `inverse`
//! (`solve(identity)`), `det`, and the dot-product `times`. Java evaluates these without fused multiply-adds, as
//! Rust does here, so the results are bit-identical.

/// An n x n matrix, row-major (Jama's `double[][] A`).
#[derive(Clone, Debug)]
pub struct Matrix {
    pub n: usize,
    pub a: Vec<f64>,
}

impl Matrix {
    pub fn zeros(n: usize) -> Matrix {
        Matrix { n, a: vec![0.0; n * n] }
    }

    pub fn identity(n: usize) -> Matrix {
        let mut m = Matrix::zeros(n);
        for i in 0..n {
            m.a[i * n + i] = 1.0;
        }
        m
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.a[i * self.n + j]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.a[i * self.n + j] = v;
    }

    /// `times(double s)`: a new matrix s * A
    pub fn times_scalar(&self, s: f64) -> Matrix {
        Matrix { n: self.n, a: self.a.iter().map(|&x| s * x).collect() }
    }

    /// `timesEquals(double s)`
    pub fn times_equals(&mut self, s: f64) {
        for x in &mut self.a {
            *x = s * *x;
        }
    }

    /// `plusEquals(Matrix B)`
    pub fn plus_equals(&mut self, b: &Matrix) {
        for (x, y) in self.a.iter_mut().zip(&b.a) {
            *x += *y;
        }
    }

    pub fn transpose(&self) -> Matrix {
        let n = self.n;
        let mut t = Matrix::zeros(n);
        for i in 0..n {
            for j in 0..n {
                t.a[j * n + i] = self.a[i * n + j];
            }
        }
        t
    }

    /// `times(Matrix B)`: column j of B copied, then a dot product per row of A.
    pub fn times(&self, b: &Matrix) -> Matrix {
        let n = self.n;
        let mut c = Matrix::zeros(n);
        let mut bcol = vec![0.0; n];
        for j in 0..n {
            for k in 0..n {
                bcol[k] = b.a[k * n + j];
            }
            for i in 0..n {
                let row = &self.a[i * n..i * n + n];
                let mut s = 0.0;
                for k in 0..n {
                    s += row[k] * bcol[k];
                }
                c.a[i * n + j] = s;
            }
        }
        c
    }

    /// `det()`: the product of the LU diagonal, signed by the pivoting.
    pub fn det(&self) -> f64 {
        Lu::new(self).det()
    }

    /// `inverse()` = `solve(identity)`; Err as Jama's "Matrix is singular." RuntimeException.
    pub fn inverse(&self) -> Result<Matrix, &'static str> {
        Lu::new(self).solve(&Matrix::identity(self.n))
    }
}

/// `LUDecomposition`
struct Lu {
    n: usize,
    lu: Vec<f64>,
    piv: Vec<usize>,
    pivsign: i32,
}

impl Lu {
    fn new(m: &Matrix) -> Lu {
        let n = m.n;
        let mut lu = m.a.clone();
        let mut piv: Vec<usize> = (0..n).collect();
        let mut pivsign = 1;
        let mut colj = vec![0.0; n];
        for j in 0..n {
            for i in 0..n {
                colj[i] = lu[i * n + j];
            }
            for i in 0..n {
                let kmax = i.min(j);
                let mut s = 0.0;
                for k in 0..kmax {
                    s += lu[i * n + k] * colj[k];
                }
                colj[i] -= s;
                lu[i * n + j] = colj[i];
            }
            let mut p = j;
            for i in j + 1..n {
                if colj[i].abs() > colj[p].abs() {
                    p = i;
                }
            }
            if p != j {
                for k in 0..n {
                    lu.swap(p * n + k, j * n + k);
                }
                piv.swap(p, j);
                pivsign = -pivsign;
            }
            if lu[j * n + j] != 0.0 {
                let d = lu[j * n + j];
                for i in j + 1..n {
                    lu[i * n + j] /= d;
                }
            }
        }
        Lu { n, lu, piv, pivsign }
    }

    fn det(&self) -> f64 {
        let mut d = self.pivsign as f64;
        for j in 0..self.n {
            d *= self.lu[j * self.n + j];
        }
        d
    }

    fn solve(&self, b: &Matrix) -> Result<Matrix, &'static str> {
        let n = self.n;
        if (0..n).any(|j| self.lu[j * n + j] == 0.0) {
            return Err("Matrix is singular.");
        }
        // B(piv, :)
        let mut x = Matrix::zeros(n);
        for (i, &p) in self.piv.iter().enumerate() {
            x.a[i * n..i * n + n].copy_from_slice(&b.a[p * n..p * n + n]);
        }
        let lu = &self.lu;
        // L * Y = B(piv, :)
        for k in 0..n {
            for i in k + 1..n {
                for j in 0..n {
                    x.a[i * n + j] -= x.a[k * n + j] * lu[i * n + k];
                }
            }
        }
        // U * X = Y
        for k in (0..n).rev() {
            for j in 0..n {
                x.a[k * n + j] /= lu[k * n + k];
            }
            for i in 0..k {
                for j in 0..n {
                    x.a[i * n + j] -= x.a[k * n + j] * lu[i * n + k];
                }
            }
        }
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inverse_and_det() {
        let m = Matrix { n: 3, a: vec![4.0, 2.0, 0.6, 2.0, 2.0, 0.5, 0.6, 0.5, 3.0] };
        let inv = m.inverse().unwrap();
        let id = m.times(&inv);
        for i in 0..3 {
            for j in 0..3 {
                assert!((id.get(i, j) - if i == j { 1.0 } else { 0.0 }).abs() < 1e-12);
            }
        }
        // det of [[4,2,.6],[2,2,.5],[.6,.5,3]]
        assert!((m.det() - 11.48).abs() < 1e-9);
        assert!(Matrix::zeros(2).inverse().is_err());
        // 200 * I inverted: exactly 1/200 on the diagonal
        let e = Matrix::identity(4).times_scalar(200.0).inverse().unwrap();
        assert_eq!(e.get(2, 2), 1.0 / 200.0);
        assert_eq!(e.get(1, 2), 0.0);
    }
}
