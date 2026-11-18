//! klib ksort.h `ks_introsort` (with its combsort fallback and final insertion sort), exactly: it is not stable,
//! and bwa's outputs depend on the order it leaves equal elements in.

fn insertsort<T: Copy>(a: &mut [T], lt: &impl Fn(&T, &T) -> bool) {
    for i in 1..a.len() {
        let mut j = i;
        while j > 0 && lt(&a[j], &a[j - 1]) {
            a.swap(j, j - 1);
            j -= 1;
        }
    }
}

fn combsort<T: Copy>(a: &mut [T], lt: &impl Fn(&T, &T) -> bool) {
    const SHRINK: f64 = 1.2473309501039786540366528676643;
    let n = a.len();
    let mut gap = n;
    loop {
        if gap > 2 {
            gap = (gap as f64 / SHRINK) as usize;
            if gap == 9 || gap == 10 {
                gap = 11;
            }
        }
        let mut do_swap = false;
        let mut i = 0;
        while i + gap < n {
            let j = i + gap;
            if lt(&a[j], &a[i]) {
                a.swap(i, j);
                do_swap = true;
            }
            i += 1;
        }
        if !(do_swap || gap > 2) {
            break;
        }
    }
    if gap != 1 {
        insertsort(a, lt);
    }
}

pub fn introsort<T: Copy>(a: &mut [T], lt: impl Fn(&T, &T) -> bool) {
    let n = a.len();
    if n < 1 {
        return;
    }
    if n == 2 {
        if lt(&a[1], &a[0]) {
            a.swap(0, 1);
        }
        return;
    }
    let mut d: i32 = 2;
    while (1usize << d) < n {
        d += 1;
    }
    let mut stack: Vec<(usize, usize, i32)> = Vec::new();
    let (mut s, mut t) = (0isize, n as isize - 1);
    d <<= 1;
    loop {
        if s < t {
            d -= 1;
            if d == 0 {
                combsort(&mut a[s as usize..=t as usize], &lt);
                t = s;
                continue;
            }
            let (mut i, mut j) = (s, t);
            let mut k = i + ((j - i) >> 1) + 1;
            if lt(&a[k as usize], &a[i as usize]) {
                if lt(&a[k as usize], &a[j as usize]) {
                    k = j;
                }
            } else {
                k = if lt(&a[j as usize], &a[i as usize]) { i } else { j };
            }
            let rp = a[k as usize];
            if k != t {
                a.swap(k as usize, t as usize);
            }
            loop {
                loop {
                    i += 1;
                    if !lt(&a[i as usize], &rp) {
                        break;
                    }
                }
                loop {
                    j -= 1;
                    if !(i <= j && lt(&rp, &a[j as usize])) {
                        break;
                    }
                }
                if j <= i {
                    break;
                }
                a.swap(i as usize, j as usize);
            }
            a.swap(i as usize, t as usize);
            if i - s > t - i {
                if i - s > 16 {
                    stack.push((s as usize, (i - 1) as usize, d));
                }
                s = if t - i > 16 { i + 1 } else { t };
            } else {
                if t - i > 16 {
                    stack.push(((i + 1) as usize, t as usize, d));
                }
                t = if i - s > 16 { i - 1 } else { s };
            }
        } else {
            match stack.pop() {
                None => {
                    insertsort(a, &lt);
                    return;
                }
                Some((l, r, dd)) => {
                    s = l as isize;
                    t = r as isize;
                    d = dd;
                }
            }
        }
    }
}
