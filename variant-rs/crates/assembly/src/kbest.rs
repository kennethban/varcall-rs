//! Port of GraphBasedKBestHaplotypeFinder (with KBestHaplotypeFinder's cycle removal),
//! using an exact port of java.util.PriorityQueue so ties resolve identically.

use crate::graph::{EId, Graph, VId};
use crate::log10::java_log10;
use crate::seqgraph::{SeqEdge, SeqGraph};
use std::cmp::Ordering;

pub struct KBestPath {
    pub bases: Vec<u8>,
    pub score: f64,
}

struct PathNode {
    last: VId,
    score: f64,
    bases: Vec<u8>,
}

/// K_BEST_HAPLOTYPE_COMPARATOR: score descending, then bases (BaseUtils.BASES_COMPARATOR) descending.
#[inline]
fn compare(a: &PathNode, b: &PathNode) -> Ordering {
    // Comparator.comparingDouble(score).reversed(): Double.compare(b.score, a.score)
    match java_double_compare(b.score, a.score) {
        Ordering::Equal => {}
        o => return o,
    }
    // BASES_COMPARATOR.reversed(): compare(b, a)
    bases_compare(&b.bases, &a.bases)
}

fn java_double_compare(a: f64, b: f64) -> Ordering {
    if a < b {
        return Ordering::Less;
    }
    if a > b {
        return Ordering::Greater;
    }
    let ab = canon_bits(a);
    let bb = canon_bits(b);
    ab.cmp(&bb)
}

fn canon_bits(x: f64) -> i64 {
    if x.is_nan() {
        0x7ff8000000000000i64
    } else {
        x.to_bits() as i64
    }
}

fn bases_compare(a: &[u8], b: &[u8]) -> Ordering {
    let m = a.len().min(b.len());
    for i in 0..m {
        let c = (a[i] as i8).cmp(&(b[i] as i8));
        if c != Ordering::Equal {
            return c;
        }
    }
    a.len().cmp(&b.len())
}

/// java.util.PriorityQueue with comparator
struct JavaPQ {
    q: Vec<PathNode>,
}

impl JavaPQ {
    fn offer(&mut self, x: PathNode) {
        let mut k = self.q.len();
        self.q.push(x);
        // siftUpUsingComparator
        while k > 0 {
            let parent = (k - 1) >> 1;
            if compare(&self.q[k], &self.q[parent]) != Ordering::Less {
                break;
            }
            self.q.swap(k, parent);
            k = parent;
        }
    }

    fn poll(&mut self) -> Option<PathNode> {
        if self.q.is_empty() {
            return None;
        }
        let n = self.q.len() - 1;
        self.q.swap(0, n);
        let result = self.q.pop().unwrap();
        if n > 0 {
            // siftDownUsingComparator(0, x) where x = old last element now at index 0
            let mut k = 0usize;
            let half = n >> 1;
            while k < half {
                let mut child = 2 * k + 1;
                let right = child + 1;
                if right < n && compare(&self.q[child], &self.q[right]) == Ordering::Greater {
                    child = right;
                }
                if compare(&self.q[k], &self.q[child]) != Ordering::Greater {
                    break;
                }
                self.q.swap(k, child);
                k = child;
            }
        }
        Some(result)
    }
}

/// KBestHaplotypeFinder.removeCyclesAndVerticesThatDontLeadToSinks
fn remove_cycles(original: &Graph<Vec<u8>, SeqEdge>, source: VId, sink: VId) -> Graph<Vec<u8>, SeqEdge> {
    let n = original.vs.len();
    let mut edges_to_remove = vec![false; original.es.len()];
    let mut verts_to_remove = vec![false; n];
    let mut parent = vec![false; n];
    // recursive DFS (iterative emulation)
    // findGuiltyVerticesAndEdgesToRemoveCycles(graph, current, sinks, edgesToRemove, verticesToRemove, parentVertices)
    struct Frame {
        v: VId,
        idx: usize,
        reaches: bool,
    }
    let found_some_path;
    if source == sink {
        found_some_path = true;
    } else {
        let mut stack: Vec<Frame> = Vec::new();
        parent[source as usize] = true;
        stack.push(Frame { v: source, idx: 0, reaches: false });
        let ret;
        loop {
            let top = stack.len() - 1;
            let v = stack[top].v;
            let out = original.out_edges(v);
            if stack[top].idx < out.len() {
                let e: EId = out[stack[top].idx];
                stack[top].idx += 1;
                let child = original.tgt(e);
                if parent[child as usize] {
                    edges_to_remove[e as usize] = true;
                } else if child == sink {
                    stack[top].reaches = true;
                } else {
                    parent[child as usize] = true;
                    stack.push(Frame { v: child, idx: 0, reaches: false });
                }
            } else {
                let f = stack.pop().unwrap();
                if !f.reaches {
                    verts_to_remove[f.v as usize] = true;
                }
                match stack.last_mut() {
                    Some(p) => p.reaches = p.reaches || f.reaches,
                    None => {
                        ret = f.reaches;
                        break;
                    }
                }
            }
        }
        found_some_path = ret;
    }
    assert!(found_some_path, "could not find any path from the source vertex to the sink vertex after removing cycles");
    let mut result = original.clone();
    for e in 0..edges_to_remove.len() {
        if edges_to_remove[e] {
            result.remove_edge(e as EId);
        }
    }
    for v in 0..n {
        if verts_to_remove[v] {
            result.remove_vertex(v as VId);
        }
    }
    result
}

pub fn find_best_haplotypes(sg: &SeqGraph, source: VId, sink: VId, max_haplotypes: usize) -> Vec<KBestPath> {
    let cleaned;
    let graph: &Graph<Vec<u8>, SeqEdge> = if sg.g.has_cycles() {
        crate::prof::hit("kbest cycle removal");
        cleaned = remove_cycles(&sg.g, source, sink);
        &cleaned
    } else {
        &sg.g
    };
    let mut result: Vec<KBestPath> = Vec::new();
    let mut pq = JavaPQ { q: Vec::new() };
    pq.offer(PathNode { last: source, score: 0.0, bases: graph.v(source).clone() });
    let mut vertex_counts = vec![0usize; graph.vs.len()];
    while result.len() < max_haplotypes {
        let p = match pq.poll() {
            Some(p) => p,
            None => break,
        };
        let v = p.last;
        if v == sink {
            result.push(KBestPath { bases: p.bases, score: p.score });
        } else {
            let c = vertex_counts[v as usize];
            vertex_counts[v as usize] += 1;
            if c < max_haplotypes {
                let out = graph.out_edges(v);
                let mut total: i32 = 0;
                for &e in out {
                    total += graph.e(e).mult;
                }
                for &e in out {
                    let t = graph.tgt(e);
                    let mult = graph.e(e).mult;
                    let score = p.score + (java_log10(mult) - java_log10(total));
                    let mut bases = Vec::with_capacity(p.bases.len() + graph.v(t).len());
                    bases.extend_from_slice(&p.bases);
                    bases.extend_from_slice(graph.v(t));
                    pq.offer(PathNode { last: t, score, bases });
                }
            }
        }
    }
    result
}
