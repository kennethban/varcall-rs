//! Port of AbstractReadThreadingGraph / ReadThreadingGraph (non-linked de Bruijn graph)
//! including chain pruning (LowWeightChainPruner) and dangling head/tail recovery.

use crate::cigar::{self, Cigar};
use crate::graph::{EId, EdgeLike, Graph, VId};
use crate::kmertable::{base_code, KmerTable, PackedKey, PackedSet};
use crate::seqgraph::{SeqEdge, SeqGraph};
use crate::sw::{Overhang, SwAligner, SwParameters};

#[derive(Clone, Copy, Debug)]
pub struct RtVertex {
    /// offset of the k-mer in the graph's byte arena
    pub off: u32,
}

/// MultiSampleEdge with numPruningSamples == 1.
#[derive(Clone, Copy, Debug)]
pub struct RtEdge {
    pub mult: i32,
    pub is_ref: bool,
    /// the single retained per-sample multiplicity (priority queue of capacity 1)
    pub pq: i32,
    /// currentSingleSampleMultiplicity
    pub cur: i32,
}

impl RtEdge {
    fn new(is_ref: bool, mult: i32) -> Self {
        RtEdge { mult, is_ref, pq: mult, cur: mult }
    }
    #[inline]
    fn inc(&mut self, n: i32) {
        self.mult += n;
        self.cur += n;
    }
    #[inline]
    fn flush(&mut self) {
        // add current, then poll the lowest when size == capacity + 1
        self.pq = self.pq.max(self.cur);
        self.cur = 0;
    }
    #[inline]
    pub fn pruning_mult(&self) -> i32 {
        self.pq
    }
}

impl EdgeLike for RtEdge {
    #[inline]
    fn is_ref(&self) -> bool {
        self.is_ref
    }
    #[inline]
    fn multiplicity(&self) -> i32 {
        self.mult
    }
}

/// SequenceForKmers
#[derive(Clone, Copy, Debug)]
pub struct SeqForKmers {
    /// arena offset of sequence start (index 0 of the sequence)
    pub seq_off: u32,
    pub seq_len: u32,
    pub start: u32,
    pub stop: u32,
    pub is_ref: bool,
}

pub struct ReadThreadingGraph {
    pub k: usize,
    pub g: Graph<RtVertex, RtEdge>,
    pub arena: Vec<u8>,
    kmer_map: KmerTable,
    non_unique: NonUniqueSet,
    ref_source_off: Option<u32>,
    /// pending sequences grouped by sample, in insertion order
    pub pending: Vec<Vec<SeqForKmers>>,
    min_matching_bases_dangling: i32,
    /// pruneFactor given to dangling-branch recovery (HaplotypeCaller 2, Mutect2 0)
    pub prune_factor: i32,
    /// `--recover-all-dangling-branches` (Mutect2's mitochondria mode): follow the heaviest edges through forks
    pub recover_all: bool,
    /// VARIANT_RS_ASM_PROFILE counters: sequences threaded, k-mers threaded, followed an existing edge, hash lookups
    /// (new edge), vertices created, find_start non-unique checks
    cnt: [u64; 6],
}

pub const PRUNE_FACTOR: i32 = 2;

/// `AdaptiveChainPruner` parameters (Mutect2 defaults: 0.001, ln(10), 4 ln(10), 100).
#[derive(Clone, Copy, Debug)]
pub struct AdaptivePruning {
    pub initial_error_rate: f64,
    pub log_odds_threshold: f64,
    pub seeding_log_odds_threshold: f64,
    pub max_unpruned_variants: usize,
}

impl Default for AdaptivePruning {
    fn default() -> Self {
        AdaptivePruning {
            initial_error_rate: 0.001,
            log_odds_threshold: 1.0 * std::f64::consts::LN_10,
            seeding_log_odds_threshold: 4.0 * std::f64::consts::LN_10,
            max_unpruned_variants: 100,
        }
    }
}
pub const MIN_DANGLING_BRANCH_LENGTH: usize = 4;
const MAX_CIGAR_COMPLEXITY: usize = 3;

/// determineNonUniqueKmers for one sequence (kmers from 0 to stop-k, as in GATK).
/// Calls `f(offset)` for every k-mer occurrence that was already seen earlier in this sequence.
pub fn for_each_non_unique_kmer(arena: &[u8], s: &SeqForKmers, k: usize, seen: &mut KmerTable, mut f: impl FnMut(u32)) {
    seen.clear();
    let stop_position = s.stop as i64 - k as i64;
    let mut i: i64 = 0;
    while i <= stop_position {
        let off = s.seq_off + i as u32;
        if !seen.insert_if_absent(arena, off, 0) {
            f(off);
        }
        i += 1;
    }
}

impl ReadThreadingGraph {
    pub fn new(k: usize, arena: Vec<u8>, pending: Vec<Vec<SeqForKmers>>, min_matching_bases_dangling: i32) -> Self {
        let total: usize = pending.iter().map(|p| p.iter().map(|s| s.seq_len as usize).sum::<usize>()).sum();
        ReadThreadingGraph {
            k,
            g: Graph::new(),
            arena,
            kmer_map: KmerTable::new(k, total.max(64)),
            non_unique: NonUniqueSet::new(k),
            ref_source_off: None,
            pending,
            min_matching_bases_dangling,
            prune_factor: PRUNE_FACTOR,
            recover_all: false,
            cnt: [0; 6],
        }
    }

    #[inline]
    fn kmer_bytes(&self, off: u32) -> &[u8] {
        &self.arena[off as usize..off as usize + self.k]
    }
    #[inline]
    pub fn vseq(&self, v: VId) -> &[u8] {
        self.kmer_bytes(self.g.v(v).off)
    }
    #[inline]
    fn suffix(&self, v: VId) -> u8 {
        self.arena[self.g.v(v).off as usize + self.k - 1]
    }

    /// buildGraphIfNecessary
    pub fn build(&mut self) {
        // preprocessReads: non-unique kmers among all pending sequences
        let _p = crate::prof::scope("rt_nonunique");
        let mut nu = NonUniqueSet::new(self.k);
        for sample in &self.pending {
            for s in sample {
                nu.add_sequence(&self.arena, s);
            }
        }
        self.non_unique = nu;
        drop(_p);
        let _p = crate::prof::scope("rt_thread");
        let pending = std::mem::take(&mut self.pending);
        for sample in &pending {
            for s in sample {
                self.thread_sequence(s);
            }
            for e in 0..self.g.es.len() {
                if self.g.es[e].alive {
                    self.g.es[e].data.flush();
                }
            }
        }
        if crate::prof::enabled() {
            let names = ["rt: sequences threaded", "rt: k-mers threaded", "rt: k-mers following an existing edge", "rt: k-mers looked up (new edge)", "rt: vertices created", "rt: find_start non-unique checks"];
            for (n, c) in names.iter().zip(self.cnt) {
                crate::prof::add(n, c);
            }
            crate::prof::add("rt: graph vertices (end of threading)", self.g.vertex_count() as u64);
            crate::prof::add("rt: graph edges (end of threading)", self.g.edge_count() as u64);
        }
    }

    fn find_start(&self, s: &SeqForKmers) -> i64 {
        if s.is_ref {
            return 0;
        }
        let mut i = s.start as i64;
        while i < s.stop as i64 - self.k as i64 {
            let off = s.seq_off + i as u32;
            if !self.non_unique.contains(&self.arena, off) {
                return i;
            }
            i += 1;
        }
        -1
    }

    fn create_vertex(&mut self, off: u32) -> VId {
        self.cnt[4] += 1;
        let v = self.g.add_vertex(RtVertex { off });
        // trackKmer
        if !self.non_unique.contains(&self.arena, off) {
            self.kmer_map.insert_if_absent(&self.arena, off, v);
        }
        v
    }

    fn get_kmer_vertex(&self, off: u32, allow_ref_source: bool) -> Option<VId> {
        let kmer = self.kmer_bytes(off);
        if !allow_ref_source {
            if let Some(rs) = self.ref_source_off {
                if self.kmer_bytes(rs) == kmer {
                    return None;
                }
            }
        }
        self.kmer_map.get(&self.arena, kmer)
    }

    fn thread_sequence(&mut self, s: &SeqForKmers) {
        self.cnt[0] += 1;
        let start_pos = self.find_start(s);
        if !s.is_ref {
            let end = if start_pos == -1 { s.stop as i64 - self.k as i64 } else { start_pos + 1 };
            self.cnt[5] += (end - s.start as i64).max(0) as u64;
        }
        if start_pos == -1 {
            return;
        }
        let start_off = s.seq_off + start_pos as u32;
        let starting_vertex = match self.get_kmer_vertex(start_off, true) {
            Some(v) => v,
            None => self.create_vertex(start_off),
        };
        // increaseCountsInMatchedKmers
        let orig_off = self.g.v(starting_vertex).off;
        self.increase_counts_in_matched_kmers(starting_vertex, orig_off, self.k as i64 - 2);
        if s.is_ref {
            assert!(self.ref_source_off.is_none(), "Found two refSources");
            self.ref_source_off = Some(s.seq_off + s.start);
        }
        let mut vertex = starting_vertex;
        let last = s.stop as i64 - self.k as i64;
        let mut i = start_pos + 1;
        while i <= last {
            vertex = self.extend_chain_by_one(vertex, s.seq_off + i as u32, s.is_ref);
            i += 1;
        }
    }

    fn increase_counts_in_matched_kmers(&mut self, vertex: VId, orig_off: u32, offset: i64) {
        // iterative version of the recursion (in-degree must be 1 to recurse, so at most one edge each step)
        let mut vertex = vertex;
        let mut offset = offset;
        loop {
            if offset == -1 {
                return;
            }
            let seq_base = self.arena[orig_off as usize + offset as usize];
            let in_deg = self.g.in_degree(vertex);
            let mut next: Option<VId> = None;
            let n = self.g.in_edges(vertex).len();
            for idx in 0..n {
                let e = self.g.in_edges(vertex)[idx];
                let prev = self.g.src(e);
                if self.suffix(prev) == seq_base && in_deg == 1 {
                    self.g.e_mut(e).inc(1);
                    next = Some(prev);
                }
            }
            match next {
                Some(p) => {
                    vertex = p;
                    offset -= 1;
                }
                None => return,
            }
        }
    }

    fn extend_chain_by_one(&mut self, prev: VId, kmer_off: u32, is_ref: bool) -> VId {
        self.cnt[1] += 1;
        let next_base = self.arena[kmer_off as usize + self.k - 1];
        let n = self.g.out_degree(prev);
        for idx in 0..n {
            let e = self.g.out_edges(prev)[idx];
            let t = self.g.tgt(e);
            if self.suffix(t) == next_base {
                self.g.e_mut(e).inc(1);
                self.cnt[2] += 1;
                return t;
            }
        }
        self.cnt[3] += 1;
        let merge = self.get_kmer_vertex(kmer_off, false);
        assert!(!(is_ref && merge.is_some()), "Found a unique vertex to merge into the reference graph");
        let next = match merge {
            Some(v) => v,
            None => self.create_vertex(kmer_off),
        };
        self.g.add_edge(prev, next, RtEdge::new(is_ref, 1));
        next
    }

    fn remove_vertex(&mut self, v: VId) {
        if self.g.remove_vertex(v) {
            let off = self.g.v(v).off as usize;
            let kmer = self.arena[off..off + self.k].to_vec();
            self.kmer_map.remove(&self.arena, &kmer);
        }
    }

    /// AbstractReadThreadingGraph.removeSingletonOrphanVertices (no ref-source exception)
    fn remove_singleton_orphans(&mut self) {
        let to_remove: Vec<VId> =
            self.g.vertex_iter().filter(|&v| self.g.in_degree(v) == 0 && self.g.out_degree(v) == 0).collect();
        for v in to_remove {
            self.remove_vertex(v);
        }
    }

    /// LowWeightChainPruner(pruneFactor).pruneLowWeightChains
    pub fn prune_low_weight_chains(&mut self, prune_factor: i32) {
        let g = &self.g;
        let mut starts: std::collections::VecDeque<VId> = g.vertex_iter().filter(|&v| g.is_source(v)).collect();
        let mut seen = vec![false; g.vs.len()];
        for &s in &starts {
            seen[s as usize] = true;
        }
        let mut to_remove: Vec<EId> = Vec::new();
        let mut chain: Vec<EId> = Vec::new();
        while let Some(cs) = starts.pop_front() {
            for &oe in g.out_edges(cs) {
                chain.clear();
                chain.push(oe);
                let first = g.src(oe);
                let mut last = g.tgt(oe);
                loop {
                    let outs = g.out_edges(last);
                    if outs.len() != 1 || g.in_degree(last) > 1 || last == first {
                        break;
                    }
                    let ne = outs[0];
                    chain.push(ne);
                    last = g.tgt(ne);
                }
                if chain.iter().all(|&e| g.e(e).pruning_mult() < prune_factor && !g.e(e).is_ref) {
                    to_remove.extend_from_slice(&chain);
                }
                if !seen[last as usize] {
                    seen[last as usize] = true;
                    starts.push_back(last);
                }
            }
        }
        for e in to_remove {
            self.g.remove_edge(e);
        }
        self.remove_singleton_orphans();
    }

    /// ChainPruner.findAllChains: breadth-first from the sources; a chain ends at a vertex without exactly one
    /// out edge, with more than one in edge, or back at its first vertex. Returns (edges, first, last).
    fn find_all_chains(&self) -> Vec<(Vec<EId>, VId, VId)> {
        let g = &self.g;
        let mut starts: std::collections::VecDeque<VId> = g.vertex_iter().filter(|&v| g.is_source(v)).collect();
        let mut seen = vec![false; g.vs.len()];
        for &s in &starts {
            seen[s as usize] = true;
        }
        let mut chains = Vec::new();
        while let Some(cs) = starts.pop_front() {
            for &oe in g.out_edges(cs) {
                let mut chain = vec![oe];
                let first = g.src(oe);
                let mut last = g.tgt(oe);
                loop {
                    let outs = g.out_edges(last);
                    if outs.len() != 1 || g.in_degree(last) > 1 || last == first {
                        break;
                    }
                    let ne = outs[0];
                    chain.push(ne);
                    last = g.tgt(ne);
                }
                chains.push((chain, first, last));
                if !seen[last as usize] {
                    seen[last as usize] = true;
                    starts.push_back(last);
                }
            }
        }
        chains
    }

    fn kmer_of(&self, v: VId) -> &[u8] {
        let off = self.g.v(v).off as usize;
        &self.arena[off..off + self.k]
    }

    /// `Path.getBases()` of a chain: the first vertex's k-mer, then the last base of each following vertex.
    fn chain_bases(&self, chain: &[EId]) -> Vec<u8> {
        let mut b = self.kmer_of(self.g.src(chain[0])).to_vec();
        for &e in chain {
            b.push(*self.kmer_of(self.g.tgt(e)).last().unwrap());
        }
        b
    }

    /// `AdaptiveChainPruner.pruneLowWeightChains` (Mutect2's default assembly pruning).
    pub fn prune_adaptive(&mut self, p: &AdaptivePruning, math: &crate::m2math::M2Math) {
        let chains = self.find_all_chains();
        if !chains.is_empty() {
            let probable = self.likely_error_chains(&chains, p, math, p.initial_error_rate);
            let error_count: i32 = probable.iter().map(|&c| self.g.e(*chains[c].0.last().unwrap()).mult).sum();
            let total_bases: i32 = chains.iter().map(|c| c.0.iter().map(|&e| self.g.e(e).mult).sum::<i32>()).sum();
            let error_rate = error_count as f64 / total_bases as f64;
            let remove: Vec<usize> = self
                .likely_error_chains(&chains, p, math, error_rate)
                .into_iter()
                .filter(|&c| !chains[c].0.iter().any(|&e| self.g.e(e).is_ref))
                .collect();
            for c in remove {
                for &e in &chains[c].0 {
                    self.g.remove_edge(e);
                }
            }
        }
        self.remove_singleton_orphans();
    }

    /// `AdaptiveChainPruner.likelyErrorChains`: indices of the chains outside the grown "good" subgraph.
    fn likely_error_chains(&self, chains: &[(Vec<EId>, VId, VId)], p: &AdaptivePruning, math: &crate::m2math::M2Math, error_rate: f64) -> Vec<usize> {
        use std::cmp::Ordering;
        let g = &self.g;
        let nv = g.vs.len();
        let mult = |e: EId| g.e(e).mult as i64;
        // chainLogOdds
        let odds: Vec<(f64, f64)> = chains
            .iter()
            .map(|(edges, first, last)| {
                let left_total: i64 = g.out_edges(*first).iter().map(|&e| mult(e)).sum();
                let right_total: i64 = g.in_edges(*last).iter().map(|&e| mult(e)).sum();
                let left_mult = mult(edges[0]);
                let right_mult = mult(*edges.last().unwrap());
                let left = if g.is_source(*first) { 0.0 } else { math.log_likelihood_ratio_err(left_total - left_mult, left_mult, error_rate) };
                let right = if g.is_sink(*last) { 0.0 } else { math.log_likelihood_ratio_err(right_total - right_mult, right_mult, error_rate) };
                (left, right)
            })
            .collect();
        let mut good_in: Vec<Vec<usize>> = vec![Vec::new(); nv];
        let mut good_out: Vec<Vec<usize>> = vec![Vec::new(); nv];
        let mut seedable: Vec<usize> = vec![0; nv];
        for (c, (edges, first, last)) in chains.iter().enumerate() {
            let first_ref = g.e(edges[0]).is_ref;
            let (l, r) = odds[c];
            if r >= p.log_odds_threshold || first_ref {
                good_in[*last as usize].push(c);
            }
            if l >= p.log_odds_threshold || first_ref {
                good_out[*first as usize].push(c);
            }
            if r >= p.seeding_log_odds_threshold && l >= p.seeding_log_odds_threshold {
                seedable[*first as usize] += 1;
                seedable[*last as usize] += 1;
            }
        }
        // priority queue: higher score first, then first vertex sequence, then chain bases
        let bases: Vec<Vec<u8>> = chains.iter().map(|c| self.chain_bases(&c.0)).collect();
        let cmp = |a: &(usize, f64), b: &(usize, f64)| -> Ordering {
            (-a.1)
                .partial_cmp(&-b.1)
                .unwrap_or(Ordering::Equal)
                .then_with(|| self.kmer_of(chains[a.0].1).cmp(self.kmer_of(chains[b.0].1)))
                .then_with(|| bases[a.0].cmp(&bases[b.0]))
        };
        let mut queue: Vec<(usize, f64)> = Vec::new();
        // max-weight chain (first maximal one in chain order)
        let key = |c: usize| -> i32 { chains[c].0.iter().map(|&e| g.e(e).mult).max().unwrap_or(0) };
        let mut best = 0usize;
        for c in 1..chains.len() {
            let o = key(c)
                .cmp(&key(best))
                .then(chains[c].0.len().cmp(&chains[best].0.len()))
                .then_with(|| self.kmer_of(chains[c].1).cmp(self.kmer_of(chains[best].1)));
            if o == Ordering::Greater {
                best = c;
            }
        }
        queue.push((best, f64::INFINITY));
        let mut processed = vec![false; nv];
        for v in 0..nv {
            if seedable[v] > 2 {
                for &c in &good_out[v] {
                    queue.push((c, odds[c].0));
                }
                for &c in &good_in[v] {
                    queue.push((c, odds[c].1));
                }
                processed[v] = true;
            }
        }
        let mut good = vec![false; chains.len()];
        let mut has_good_out = vec![false; nv];
        let mut variant_count = 0usize;
        while !queue.is_empty() && variant_count <= p.max_unpruned_variants {
            // poll the minimum under `cmp`
            let mut mi = 0;
            for i in 1..queue.len() {
                if cmp(&queue[i], &queue[mi]) == Ordering::Less {
                    mi = i;
                }
            }
            let (c, _) = queue.swap_remove(mi);
            if good[c] {
                continue;
            }
            good[c] = true;
            let first = chains[c].1;
            let new_variant = has_good_out[first as usize];
            has_good_out[first as usize] = true;
            if new_variant {
                variant_count += 1;
            }
            if new_variant && variant_count > p.max_unpruned_variants {
                continue;
            }
            for v in [chains[c].1, chains[c].2] {
                if !processed[v as usize] {
                    for &c2 in &good_out[v as usize] {
                        queue.push((c2, odds[c2].0));
                    }
                    for &c2 in &good_in[v as usize] {
                        queue.push((c2, odds[c2].1));
                    }
                    processed[v as usize] = true;
                }
            }
        }
        (0..chains.len()).filter(|&c| !good[c]).collect()
    }

    /// ReadThreadingGraph.isLowQualityGraph
    pub fn is_low_quality_graph(&self) -> bool {
        self.non_unique.len() * 4 > self.kmer_map.len()
    }

    fn heaviest_incoming_edge(&self, v: VId) -> EId {
        let inc = self.g.in_edges(v);
        if inc.len() == 1 {
            return inc[0];
        }
        let mut best = inc[0];
        for &e in &inc[1..] {
            if self.g.e(e).mult > self.g.e(best).mult {
                best = e;
            }
        }
        best
    }
    fn heaviest_outgoing_edge(&self, v: VId) -> EId {
        let out = self.g.out_edges(v);
        if out.len() == 1 {
            return out[0];
        }
        let mut best = out[0];
        for &e in &out[1..] {
            if self.g.e(e).mult > self.g.e(best).mult {
                best = e;
            }
        }
        best
    }

    /// findPath: `heaviest` (recoverAll) takes the heaviest incoming / outgoing edge, else the only one
    fn find_path(
        &self,
        vertex: VId,
        prune_factor: i32,
        done: impl Fn(&Self, VId) -> bool,
        return_path: impl Fn(&Self, VId) -> bool,
        upwards: bool,
        heaviest: bool,
    ) -> Option<Vec<VId>> {
        // path is built with addFirst; we build reversed and reverse at the end
        let mut rev_path: Vec<VId> = Vec::new();
        let mut visited: Vec<VId> = Vec::new();
        let mut v = vertex;
        while !done(self, v) {
            let edge = match (heaviest, upwards) {
                (true, true) => Some(self.heaviest_incoming_edge(v)),
                (true, false) => Some(self.heaviest_outgoing_edge(v)),
                (false, true) => self.g.incoming_edge_of(v),
                (false, false) => self.g.outgoing_edge_of(v),
            };
            let edge = edge.expect("findPath: missing edge");
            if self.g.e(edge).pruning_mult() < prune_factor {
                visited.extend_from_slice(&rev_path);
                rev_path.clear();
            } else {
                rev_path.push(v);
            }
            v = if upwards { self.g.src(edge) } else { self.g.tgt(edge) };
            if rev_path.contains(&v) || visited.contains(&v) {
                // GATK prints "Dangling End recovery killed because of a loop (findPath)" to stderr
                crate::prof::hit("findPath loop killed");
                return None;
            }
        }
        rev_path.push(v);
        rev_path.reverse();
        if return_path(self, v) {
            Some(rev_path)
        } else {
            None
        }
    }

    fn reference_path_downwards(&self, start: VId, blacklist: Option<EId>) -> Vec<VId> {
        let mut path = Vec::new();
        let mut v = Some(start);
        while let Some(x) = v {
            path.push(x);
            v = self.g.next_reference_vertex(Some(x), true, blacklist);
            if path.len() > self.g.vs.len() + 1 {
                break; // guard against infinite loops (GATK would hang)
            }
        }
        path
    }

    fn reference_path_upwards(&self, start: VId) -> Vec<VId> {
        let mut path = Vec::new();
        let mut v = Some(start);
        while let Some(x) = v {
            path.push(x);
            v = self.g.prev_reference_vertex(x);
            if path.len() > self.g.vs.len() + 1 {
                break;
            }
        }
        path
    }

    fn bases_for_path(&self, path: &[VId], expand_source: bool) -> Vec<u8> {
        let mut r = Vec::with_capacity(path.len() + self.k);
        for &v in path {
            if expand_source && self.g.is_source(v) {
                let s = self.vseq(v);
                r.extend(s.iter().rev());
            } else {
                r.push(self.suffix(v));
            }
        }
        r
    }

    pub fn recover_dangling_tails(&mut self, aligner: &mut SwAligner, params: &SwParameters) {
        let n = self.g.vs.len();
        for vi in 0..n {
            let v = vi as VId;
            if !self.g.contains_vertex(v) {
                continue;
            }
            if self.g.out_degree(v) == 0 && !self.g.is_ref_sink(v) {
                self.recover_dangling_tail(v, aligner, params);
            }
        }
    }

    pub fn recover_dangling_heads(&mut self, aligner: &mut SwAligner, params: &SwParameters) {
        let heads: Vec<VId> =
            self.g.vertex_iter().filter(|&v| self.g.in_degree(v) == 0 && !self.g.is_ref_source(v)).collect();
        for v in heads {
            self.recover_dangling_head(v, aligner, params);
        }
    }

    fn recover_dangling_tail(&mut self, v: VId, aligner: &mut SwAligner, params: &SwParameters) -> i32 {
        let min_tail_path_length = MIN_DANGLING_BRANCH_LENGTH.max(1);
        // findPathUpwardsToLowestCommonAncestor(vertex, pruneFactor, giveUpAtBranch = !recoverAll)
        let has_incident_ref_edge = |s: &Self, x: VId| s.g.in_edges(x).iter().any(|&e| s.g.e(e).is_ref());
        let alt_path = if self.recover_all {
            self.find_path(
                v,
                self.prune_factor,
                |s, x| has_incident_ref_edge(s, x) || s.g.in_degree(x) == 0,
                |s, x| s.g.out_degree(x) > 1 && has_incident_ref_edge(s, x),
                true,
                true,
            )
        } else {
            self.find_path(v, self.prune_factor, |s, x| s.g.in_degree(x) != 1 || s.g.out_degree(x) >= 2, |s, x| s.g.out_degree(x) > 1, true, false)
        };
        let alt_path = match alt_path {
            Some(p) => p,
            None => return 0,
        };
        if self.g.is_ref_source(alt_path[0]) || alt_path.len() < min_tail_path_length + 1 {
            return 0;
        }
        let blacklist = self.heaviest_incoming_edge(alt_path[1]);
        let ref_path = self.reference_path_downwards(alt_path[0], Some(blacklist));
        let ref_bases = self.bases_for_path(&ref_path, false);
        let alt_bases = self.bases_for_path(&alt_path, false);
        let aln = aligner.align(&ref_bases, &alt_bases, params, Overhang::LeadingIndel);
        let cig = cigar::remove_trailing_deletions(aln.cigar);
        if !cigar_is_okay_to_merge(&cig, false, true) {
            return 0;
        }
        // mergeDanglingTail
        let last = *cig.last().unwrap();
        assert!(last.0 == b'M', "The last Cigar element must be an M");
        let last_ref_index = cigar::reference_length(&cig) - 1;
        let matching_suffix =
            (longest_suffix_match(&ref_bases, &alt_bases, last_ref_index) as i64).min(last.1 as i64);
        let fail = if self.min_matching_bases_dangling >= 0 {
            matching_suffix < self.min_matching_bases_dangling as i64
        } else {
            matching_suffix == 0
        };
        if fail {
            return 0;
        }
        let alt_index_to_merge = (cigar::read_length(&cig) - matching_suffix - 1).max(0);
        let first_is_del = cig[0].0 == b'D';
        let must_handle = first_is_del && (cig[0].1 as i64 + matching_suffix == last_ref_index + 1);
        let ref_index_to_merge = last_ref_index - matching_suffix + 1 + if must_handle { 1 } else { 0 };
        if ref_index_to_merge == 0 {
            return 0;
        }
        let a = alt_path[alt_index_to_merge as usize];
        let b = ref_path[ref_index_to_merge as usize];
        crate::prof::hit("dangling tail merged");
        self.g.add_edge(a, b, RtEdge::new(false, 1));
        1
    }

    fn recover_dangling_head(&mut self, v: VId, aligner: &mut SwAligner, params: &SwParameters) -> i32 {
        if self.g.in_degree(v) != 0 {
            // GATK would throw IllegalStateException here
            return 0;
        }
        // findPathDownwardsToHighestCommonDescendantOfReference(vertex, pruneFactor, giveUpAtBranch = !recoverAll)
        let alt_path = if self.recover_all {
            self.find_path(v, self.prune_factor, |s, x| s.g.is_reference_node(x) || s.g.out_degree(x) == 0, |s, x| s.g.is_reference_node(x), false, true)
        } else {
            self.find_path(v, self.prune_factor, |s, x| s.g.is_reference_node(x) || s.g.out_degree(x) != 1, |s, x| s.g.is_reference_node(x), false, false)
        };
        let mut alt_path = match alt_path {
            Some(p) => p,
            None => return 0,
        };
        if self.g.is_ref_sink(alt_path[0]) || alt_path.len() < MIN_DANGLING_BRANCH_LENGTH + 1 {
            return 0;
        }
        let ref_path = self.reference_path_upwards(alt_path[0]);
        let ref_bases = self.bases_for_path(&ref_path, true);
        let alt_bases = self.bases_for_path(&alt_path, true);
        let aln = aligner.align(&ref_bases, &alt_bases, params, Overhang::LeadingIndel);
        let cig = cigar::remove_trailing_deletions(aln.cigar);
        if !cigar_is_okay_to_merge(&cig, true, false) {
            return 0;
        }
        if self.min_matching_bases_dangling >= 0 {
            unimplemented!("non-legacy dangling head merging (minMatchingBasesToDanglingEndRecovery >= 0)");
        }
        // mergeDanglingHeadLegacy
        let first = cig[0];
        assert!(first.0 == b'M', "The first Cigar element must be an M");
        let indexes_to_merge = self.best_prefix_match_legacy(&ref_bases, &alt_bases, first.1 as usize);
        if indexes_to_merge <= 0 {
            return 0;
        }
        let itm = indexes_to_merge as usize;
        if itm >= ref_path.len() - 1 {
            return 0;
        }
        if itm >= alt_path.len() {
            let num = itm - alt_path.len() + 2;
            if !self.extend_dangling_path_against_reference(&mut alt_path, &ref_path, num, &cig) {
                return 0;
            }
            crate::prof::hit("dangling head path extended");
        }
        crate::prof::hit("dangling head merged");
        self.g.add_edge(ref_path[itm + 1], alt_path[itm], RtEdge::new(false, 1));
        1
    }

    fn best_prefix_match_legacy(&self, path1: &[u8], path2: &[u8], max_index: usize) -> i64 {
        let max_mismatches = 1.max(max_index / self.k);
        let mut mismatches = 0;
        let mut last_good_index: i64 = -1;
        for index in 0..max_index {
            if path1[index] != path2[index] {
                mismatches += 1;
                if mismatches > max_mismatches {
                    return -1;
                }
                last_good_index = index as i64;
            }
        }
        last_good_index
    }

    fn extend_dangling_path_against_reference(
        &mut self,
        dangling_path: &mut Vec<VId>,
        reference_path: &[VId],
        num_nodes_to_extend: usize,
        elements: &Cigar,
    ) -> bool {
        let index_of_last_dangling_node = dangling_path.len() as i64 - 1;
        let offset: i64 = elements
            .iter()
            .map(|&(op, len)| {
                (if cigar::consumes_ref(op) { len as i64 } else { 0 })
                    - (if cigar::consumes_read(op) { len as i64 } else { 0 })
            })
            .sum();
        let index_of_ref_node_to_use = index_of_last_dangling_node + offset + num_nodes_to_extend as i64;
        if index_of_ref_node_to_use >= reference_path.len() as i64 {
            return false;
        }
        let dangling_source = dangling_path.pop().unwrap();
        let ref_source_seq = self.vseq(reference_path[index_of_ref_node_to_use as usize]).to_vec();
        let mut seq_to_extend: Vec<u8> = Vec::with_capacity(num_nodes_to_extend + self.k);
        seq_to_extend.extend_from_slice(&ref_source_seq[..num_nodes_to_extend]);
        seq_to_extend.extend_from_slice(self.vseq(dangling_source));
        let source_edge = self.heaviest_outgoing_edge(dangling_source);
        let source_mult = self.g.e(source_edge).mult;
        let mut prev_v = self.g.tgt(source_edge);
        self.g.remove_edge_between(dangling_source, prev_v);
        for i in (1..=num_nodes_to_extend).rev() {
            let off = self.arena.len() as u32;
            self.arena.extend_from_slice(&seq_to_extend[i..i + self.k]);
            let new_v = self.g.add_vertex(RtVertex { off });
            let ne = self.g.add_edge(new_v, prev_v, RtEdge::new(false, 1)).unwrap();
            self.g.e_mut(ne).mult = source_mult;
            dangling_path.push(new_v);
            prev_v = new_v;
        }
        true
    }

    /// BaseGraph.removePathsNotConnectedToRef
    pub fn remove_paths_not_connected_to_ref(&mut self) {
        let src = self.g.reference_source_vertex().expect("Graph must have ref source and sink vertices");
        let sink = self.g.reference_sink_vertex().expect("Graph must have ref source and sink vertices");
        let fwd = self.g.reachable(src, false, true);
        let bwd = self.g.reachable(sink, true, false);
        let to_remove: Vec<VId> =
            self.g.vertex_iter().filter(|&v| !(fwd[v as usize] && bwd[v as usize])).collect();
        for v in to_remove {
            self.remove_vertex(v);
        }
    }

    /// BaseGraph.toSequenceGraph
    pub fn to_sequence_graph(&self) -> SeqGraph {
        let mut sg = SeqGraph::new(self.k);
        let mut map = vec![u32::MAX; self.g.vs.len()];
        for v in self.g.vertex_iter() {
            let seq = if self.g.is_source(v) { self.vseq(v).to_vec() } else { vec![self.suffix(v)] };
            map[v as usize] = sg.g.add_vertex(seq);
        }
        for e in self.g.edges() {
            let d = self.g.e(e);
            sg.g.add_edge(
                map[self.g.src(e) as usize],
                map[self.g.tgt(e) as usize],
                SeqEdge { mult: d.mult, is_ref: d.is_ref },
            );
        }
        sg
    }
}

pub fn cigar_is_okay_to_merge(c: &Cigar, require_first_m: bool, require_last_m: bool) -> bool {
    let n = c.len();
    if n == 0 || n > MAX_CIGAR_COMPLEXITY {
        return false;
    }
    if require_first_m && c[0].0 != b'M' {
        return false;
    }
    if require_last_m && c[n - 1].0 != b'M' {
        return false;
    }
    true
}

pub fn longest_suffix_match(seq: &[u8], kmer: &[u8], seq_start: i64) -> usize {
    for len in 1..=kmer.len() {
        let seq_i = seq_start - len as i64 + 1;
        let kmer_i = kmer.len() - len;
        if seq_i < 0 || seq[seq_i as usize] != kmer[kmer_i] {
            return len - 1;
        }
    }
    kmer.len()
}

enum PackedNu {
    Small(PackedSet<u64>, PackedSet<u64>),
    Large(PackedSet<u128>, PackedSet<u128>),
    None,
}

/// The set of non-unique k-mers (ReadThreadingGraph.nonUniqueKmers). K-mers made only of
/// ACGT are stored 2-bit packed (k <= 64); anything else goes through a byte-slice table.
pub struct NonUniqueSet {
    k: usize,
    packed: PackedNu,
    slow_seen: KmerTable,
    slow_nu: KmerTable,
}

fn scan_sequence<K: PackedKey>(
    arena: &[u8],
    s: &SeqForKmers,
    k: usize,
    seen: &mut PackedSet<K>,
    nu: &mut PackedSet<K>,
    slow_seen: &mut KmerTable,
    slow_nu: &mut KmerTable,
) {
    seen.clear();
    slow_seen.clear();
    let stop = s.stop as usize;
    if stop < k {
        return;
    }
    let base = s.seq_off as usize;
    let seq = &arena[base..base + stop];
    let mask = K::mask_for(k);
    let mut code = K::default();
    let mut last_bad: i64 = -1;
    for (p, &b) in seq.iter().enumerate() {
        match base_code(b) {
            Some(c) => code = code.push(c, mask),
            None => {
                last_bad = p as i64;
                code = code.push(0, mask);
            }
        }
        if p + 1 >= k {
            let start = p + 1 - k;
            if last_bad >= start as i64 {
                crate::prof::hit("non-ACGT kmer (slow path)");
                let off = (base + start) as u32;
                if !slow_seen.insert_if_absent(arena, off, 0) {
                    slow_nu.insert_if_absent(arena, off, 0);
                }
            } else if !seen.insert(code) {
                nu.insert(code);
            }
        }
    }
}

impl NonUniqueSet {
    pub fn new(k: usize) -> Self {
        let packed = if k <= 32 {
            PackedNu::Small(PackedSet::new(256), PackedSet::new(64))
        } else if k <= 64 {
            PackedNu::Large(PackedSet::new(256), PackedSet::new(64))
        } else {
            PackedNu::None
        };
        NonUniqueSet { k, packed, slow_seen: KmerTable::new(k, 256), slow_nu: KmerTable::new(k, 64) }
    }

    /// determineNonUniqueKmers(sequence) unioned into the set (k-mers from 0 to stop-k)
    pub fn add_sequence(&mut self, arena: &[u8], s: &SeqForKmers) {
        let k = self.k;
        match &mut self.packed {
            PackedNu::Small(seen, nu) => scan_sequence(arena, s, k, seen, nu, &mut self.slow_seen, &mut self.slow_nu),
            PackedNu::Large(seen, nu) => scan_sequence(arena, s, k, seen, nu, &mut self.slow_seen, &mut self.slow_nu),
            PackedNu::None => {
                let slow_nu = &mut self.slow_nu;
                for_each_non_unique_kmer(arena, s, k, &mut self.slow_seen, |off| {
                    slow_nu.insert_if_absent(arena, off, 0);
                });
            }
        }
    }

    fn pack<K: PackedKey>(kmer: &[u8], k: usize) -> Option<K> {
        let mask = K::mask_for(k);
        let mut code = K::default();
        for &b in kmer {
            code = code.push(base_code(b)?, mask);
        }
        Some(code)
    }

    pub fn contains(&self, arena: &[u8], off: u32) -> bool {
        let kmer = &arena[off as usize..off as usize + self.k];
        match &self.packed {
            PackedNu::Small(_, nu) => match Self::pack::<u64>(kmer, self.k) {
                Some(c) => nu.contains(c),
                None => self.slow_nu.contains(arena, kmer),
            },
            PackedNu::Large(_, nu) => match Self::pack::<u128>(kmer, self.k) {
                Some(c) => nu.contains(c),
                None => self.slow_nu.contains(arena, kmer),
            },
            PackedNu::None => self.slow_nu.contains(arena, kmer),
        }
    }

    pub fn len(&self) -> usize {
        let p = match &self.packed {
            PackedNu::Small(_, nu) => nu.len(),
            PackedNu::Large(_, nu) => nu.len(),
            PackedNu::None => 0,
        };
        p + self.slow_nu.len()
    }
}
