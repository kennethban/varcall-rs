//! Port of SeqGraph and its simplification transformers (zipLinearChains, MergeDiamonds,
//! MergeTails, SplitCommonSuffices, MergeCommonSuffices) plus the BaseGraph cleanup routines.

use crate::graph::{EId, EdgeLike, Graph, VId};

#[derive(Clone, Copy, Debug)]
pub struct SeqEdge {
    pub mult: i32,
    pub is_ref: bool,
}

impl EdgeLike for SeqEdge {
    #[inline]
    fn is_ref(&self) -> bool {
        self.is_ref
    }
    #[inline]
    fn multiplicity(&self) -> i32 {
        self.mult
    }
}

impl SeqEdge {
    fn add(&mut self, o: &SeqEdge) {
        self.mult += o.mult;
        self.is_ref = self.is_ref || o.is_ref;
    }
}

#[derive(Clone)]
pub struct SeqGraph {
    pub k: usize,
    pub g: Graph<Vec<u8>, SeqEdge>,
}

const MAX_REASONABLE_SIMPLIFICATION_CYCLES: usize = 100;
const MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES: usize = 10;

fn common_maximum_prefix_length(seqs: &[&[u8]]) -> usize {
    let min_len = seqs.iter().map(|s| s.len()).min().unwrap_or(0);
    for i in 0..min_len {
        let b = seqs[0][i];
        for s in &seqs[1..] {
            if s[i] != b {
                return i;
            }
        }
    }
    min_len
}

fn common_maximum_suffix_length(seqs: &[&[u8]], min_len: usize) -> usize {
    for suffix_len in 0..min_len {
        let s0 = seqs[0];
        let b = s0[s0.len() - suffix_len - 1];
        for s in &seqs[1..] {
            if s[s.len() - suffix_len - 1] != b {
                return suffix_len;
            }
        }
    }
    min_len
}

impl SeqGraph {
    pub fn new(k: usize) -> Self {
        SeqGraph { k, g: Graph::new() }
    }

    #[inline]
    pub fn seq(&self, v: VId) -> &[u8] {
        self.g.v(v)
    }

    // ----------------------------------------------------------------------------------
    // BaseGraph cleanup routines
    // ----------------------------------------------------------------------------------

    /// BaseGraph.cleanNonRefPaths
    pub fn clean_non_ref_paths(&mut self) {
        let (src, sink) = match (self.g.reference_source_vertex(), self.g.reference_sink_vertex()) {
            (Some(a), Some(b)) => (a, b),
            _ => return,
        };
        // upstream of the ref source
        let mut work: Vec<EId> = self.g.in_edges(src).to_vec();
        while let Some(e) = work.pop() {
            if !self.g.edge_alive(e) {
                continue;
            }
            if !self.g.e(e).is_ref {
                let s = self.g.src(e);
                work.extend_from_slice(self.g.in_edges(s));
                self.g.remove_edge(e);
            }
        }
        let sink2 = self.g.reference_sink_vertex().unwrap_or(sink);
        let mut work: Vec<EId> = self.g.out_edges(sink2).to_vec();
        while let Some(e) = work.pop() {
            if !self.g.edge_alive(e) {
                continue;
            }
            if !self.g.e(e).is_ref {
                let t = self.g.tgt(e);
                work.extend_from_slice(self.g.out_edges(t));
                self.g.remove_edge(e);
            }
        }
        self.remove_singleton_orphan_vertices();
    }

    /// BaseGraph.removeSingletonOrphanVertices
    pub fn remove_singleton_orphan_vertices(&mut self) {
        let to_remove: Vec<VId> = self
            .g
            .vertex_iter()
            .filter(|&v| self.g.in_degree(v) == 0 && self.g.out_degree(v) == 0 && !self.g.is_ref_source(v))
            .collect();
        for v in to_remove {
            self.g.remove_vertex(v);
        }
    }

    /// BaseGraph.removeVerticesNotConnectedToRefRegardlessOfEdgeDirection
    pub fn remove_vertices_not_connected_to_ref_regardless_of_edge_direction(&mut self) {
        let keep = match self.g.reference_source_vertex() {
            Some(r) => self.g.reachable(r, true, true),
            None => vec![false; self.g.vs.len()],
        };
        let to_remove: Vec<VId> = self.g.vertex_iter().filter(|&v| !keep[v as usize]).collect();
        for v in to_remove {
            self.g.remove_vertex(v);
        }
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
            self.g.remove_vertex(v);
        }
    }

    /// BaseGraph.getReferenceBytes(source, sink, true, true)
    pub fn reference_bytes(&self, from: VId, to: VId) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(self.seq(from));
        let mut v = self.g.next_reference_vertex(Some(from), false, None);
        while let Some(x) = v {
            if x == to {
                break;
            }
            bytes.extend_from_slice(self.seq(x));
            v = self.g.next_reference_vertex(Some(x), false, None);
        }
        if v == Some(to) {
            bytes.extend_from_slice(self.seq(to));
        }
        bytes
    }

    // ----------------------------------------------------------------------------------
    // simplification
    // ----------------------------------------------------------------------------------

    pub fn simplify_graph(&mut self) {
        self.zip_linear_chains();
        let mut prev: Option<GraphSnapshot> = None;
        let mut i = 0usize;
        loop {
            if i > MAX_REASONABLE_SIMPLIFICATION_CYCLES {
                panic!("Infinite loop detected in simplification routines for kmer graph {}", self.k);
            }
            let did = self.simplify_graph_once();
            if !did {
                break;
            }
            if i > 5 {
                crate::prof::hit("simplifyGraph >5 cycles");
                let snap = GraphSnapshot::of(self);
                if let Some(p) = &prev {
                    if p.graph_equals(&snap) {
                        crate::prof::hit("simplifyGraph graphEquals break");
                        break;
                    }
                }
                prev = Some(snap);
            }
            i += 1;
        }
    }

    fn simplify_graph_once(&mut self) -> bool {
        let mut did = false;
        did |= self.transform_until_complete(Transform::MergeDiamonds);
        did |= self.transform_until_complete(Transform::MergeTails);
        did |= self.transform_until_complete(Transform::SplitCommonSuffices);
        did |= self.transform_until_complete(Transform::MergeCommonSuffices);
        did |= self.zip_linear_chains();
        did
    }

    fn transform_until_complete(&mut self, t: Transform) -> bool {
        let mut did_at_least_one = false;
        let mut already_split: Vec<bool> = Vec::new();
        loop {
            let mut found = false;
            let n = self.g.vs.len();
            let mut vi = 0usize;
            while vi < n {
                let v = vi as VId;
                vi += 1;
                if !self.g.contains_vertex(v) {
                    continue;
                }
                found = match t {
                    Transform::MergeDiamonds => self.merge_diamonds(v),
                    Transform::MergeTails => self.merge_tails(v),
                    Transform::SplitCommonSuffices => {
                        if already_split.len() < self.g.vs.len() {
                            already_split.resize(self.g.vs.len(), false);
                        }
                        if already_split[v as usize] {
                            false
                        } else {
                            already_split[v as usize] = true;
                            self.split_common_suffix(v)
                        }
                    }
                    Transform::MergeCommonSuffices => self.shared_sequence_merge(v),
                };
                if found {
                    did_at_least_one = true;
                    break;
                }
            }
            if !found {
                break;
            }
        }
        did_at_least_one
    }

    // ---------------------------- zipLinearChains ----------------------------

    pub fn zip_linear_chains(&mut self) -> bool {
        let starts: Vec<VId> = self.g.vertex_iter().filter(|&v| self.is_linear_chain_start(v)).collect();
        if starts.is_empty() {
            return false;
        }
        let mut merged = false;
        for s in starts {
            let chain = self.trace_linear_chain(s);
            merged |= self.merge_linear_chain(&chain);
        }
        merged
    }

    fn is_linear_chain_start(&self, v: VId) -> bool {
        self.g.out_degree(v) == 1
            && (self.g.in_degree(v) != 1 || self.g.out_degree(self.g.src(self.g.in_edges(v)[0])) > 1)
    }

    fn trace_linear_chain(&self, start: VId) -> Vec<VId> {
        let mut chain = vec![start];
        let mut last_is_ref = self.g.is_reference_node(start);
        let mut last = start;
        loop {
            if self.g.out_degree(last) != 1 {
                break;
            }
            let target = self.g.tgt(self.g.out_edges(last)[0]);
            if self.g.in_degree(target) != 1 || last == target {
                break;
            }
            let target_is_ref = self.g.is_reference_node(target);
            if last_is_ref != target_is_ref {
                break;
            }
            chain.push(target);
            last = target;
            last_is_ref = target_is_ref;
        }
        chain
    }

    fn merge_linear_chain(&mut self, chain: &[VId]) -> bool {
        let first = chain[0];
        let last = *chain.last().unwrap();
        if first == last {
            return false;
        }
        let mut seq = Vec::new();
        for &v in chain {
            seq.extend_from_slice(self.seq(v));
        }
        let added = self.g.add_vertex(seq);
        let outs: Vec<EId> = self.g.out_edges(last).to_vec();
        for e in outs {
            let t = self.g.tgt(e);
            let d = *self.g.e(e);
            self.g.add_edge(added, t, d);
        }
        let ins: Vec<EId> = self.g.in_edges(first).to_vec();
        for e in ins {
            let s = self.g.src(e);
            let d = *self.g.e(e);
            self.g.add_edge(s, added, d);
        }
        for &v in chain {
            self.g.remove_vertex(v);
        }
        true
    }

    // ---------------------------- MergeDiamonds ----------------------------

    fn merge_diamonds(&mut self, top: VId) -> bool {
        let middles = self.g.outgoing_vertices(top);
        if middles.len() <= 1 {
            return false;
        }
        let mut bottom: Option<VId> = None;
        for &mi in &middles {
            if self.g.out_degree(mi) < 1 {
                return false;
            }
            if self.g.in_degree(mi) != 1 {
                return false;
            }
            for mt in self.g.outgoing_vertices(mi) {
                match bottom {
                    None => bottom = Some(mt),
                    Some(b) if b != mt => return false,
                    _ => {}
                }
            }
        }
        let bottom = bottom.unwrap();
        if self.g.in_degree(bottom) != middles.len() {
            return false;
        }
        let splitter = Splitter::new(self, &middles);
        if !(splitter.prefix.len() >= 1 || splitter.suffix.len() >= 1) {
            return false;
        }
        crate::prof::hit("MergeDiamonds applied");
        splitter.split_and_update(self, Some(top), Some(bottom));
        true
    }

    // ---------------------------- MergeTails ----------------------------

    fn merge_tails(&mut self, top: VId) -> bool {
        let tails = self.g.outgoing_vertices(top);
        if tails.len() <= 1 {
            return false;
        }
        for &t in &tails {
            if !self.g.is_sink(t) || self.g.in_degree(t) > 1 {
                return false;
            }
        }
        let splitter = Splitter::new(self, &tails);
        if splitter.suffix.len() < MIN_COMMON_SEQUENCE_TO_MERGE_SOURCE_SINK_VERTICES {
            return false;
        }
        crate::prof::hit("MergeTails applied");
        splitter.split_and_update(self, Some(top), None);
        true
    }

    // ---------------------------- CommonSuffixSplitter ----------------------------

    fn split_common_suffix(&mut self, v: VId) -> bool {
        let to_split = self.g.incoming_vertices(v);
        if to_split.len() < 2 {
            return false;
        }
        // safeToSplit
        let outgoing_of_bot = self.g.outgoing_vertices(v);
        for &m in &to_split {
            let outs = self.g.out_edges(m);
            if m == v || outs.len() != 1 || !self.g.outgoing_vertices(m).contains(&v) {
                return false;
            }
            if outgoing_of_bot.contains(&m) {
                return false;
            }
        }
        // commonSuffix
        let suffix: Vec<u8> = {
            let seqs: Vec<&[u8]> = to_split.iter().map(|&x| self.seq(x)).collect();
            let min = seqs.iter().map(|s| s.len()).min().unwrap_or(0);
            let suffix_len = common_maximum_suffix_length(&seqs, min);
            let s0 = seqs[0];
            s0[s0.len() - suffix_len..].to_vec()
        };
        if suffix.is_empty() {
            return false;
        }
        // wouldEliminateRefSource
        for &m in &to_split {
            if self.g.is_ref_source(m) {
                if self.seq(m).len() == suffix.len() {
                    return false;
                }
                break;
            }
        }
        // allVerticesAreTheCommonSuffix
        if to_split.iter().all(|&m| self.seq(m).len() == suffix.len()) {
            return false;
        }
        crate::prof::hit("SplitCommonSuffices applied");
        let mut edges_to_remove: Vec<EId> = Vec::new();
        for &mid in &to_split {
            let suffix_v = self.g.add_vertex(suffix.clone());
            let mid_len = self.seq(mid).len();
            let prefix_size = mid_len - suffix.len();
            let out = self.g.outgoing_edge_of(mid).unwrap();
            let incoming_target;
            if prefix_size == 0 {
                incoming_target = suffix_v;
            } else {
                let pseq = self.seq(mid)[..prefix_size].to_vec();
                let prefix_v = self.g.add_vertex(pseq);
                incoming_target = prefix_v;
                let out_ref = self.g.e(out).is_ref;
                self.g.add_edge(prefix_v, suffix_v, SeqEdge { mult: 1, is_ref: out_ref });
                edges_to_remove.push(out);
            }
            let out_tgt = self.g.tgt(out);
            let out_d = *self.g.e(out);
            self.g.add_edge(suffix_v, out_tgt, out_d);
            let ins: Vec<EId> = self.g.in_edges(mid).to_vec();
            for e in ins {
                let s = self.g.src(e);
                let d = *self.g.e(e);
                self.g.add_edge(s, incoming_target, d);
                edges_to_remove.push(e);
            }
        }
        for &m in &to_split {
            self.g.remove_vertex(m);
        }
        for e in edges_to_remove {
            self.g.remove_edge(e);
        }
        true
    }

    // ---------------------------- SharedSequenceMerger ----------------------------

    fn shared_sequence_merge(&mut self, v: VId) -> bool {
        let prevs = self.g.incoming_vertices(v);
        if prevs.is_empty() {
            return false;
        }
        let first = prevs[0];
        for &p in &prevs {
            if self.seq(p) != self.seq(first) {
                return false;
            }
            let pouts = self.g.outgoing_vertices(p);
            if pouts.len() != 1 {
                return false;
            }
            if pouts[0] != v {
                return false;
            }
            if self.g.in_degree(p) == 0 {
                return false;
            }
        }
        crate::prof::hit("MergeCommonSuffices applied");
        let mut seq = self.seq(first).to_vec();
        seq.extend_from_slice(self.seq(v));
        let new_v = self.g.add_vertex(seq);
        let mut edges_to_remove: Vec<EId> = Vec::new();
        for &p in &prevs {
            let ins: Vec<EId> = self.g.in_edges(p).to_vec();
            for e in ins {
                let s = self.g.src(e);
                let d = *self.g.e(e);
                if self.g.add_edge(s, new_v, d).is_none() {
                    crate::prof::hit("MergeCommonSuffices: duplicate edge dropped");
                }
                edges_to_remove.push(e);
            }
        }
        let outs: Vec<EId> = self.g.out_edges(v).to_vec();
        for e in outs {
            let t = self.g.tgt(e);
            let d = *self.g.e(e);
            self.g.add_edge(new_v, t, d);
        }
        for &p in &prevs {
            self.g.remove_vertex(p);
        }
        self.g.remove_vertex(v);
        for e in edges_to_remove {
            self.g.remove_edge(e);
        }
        true
    }
}

#[derive(Clone, Copy)]
enum Transform {
    MergeDiamonds,
    MergeTails,
    SplitCommonSuffices,
    MergeCommonSuffices,
}

/// Port of SharedVertexSequenceSplitter. The temporary "split graph" is represented by
/// the ordered outgoing edges of the prefix vertex and the ordered incoming edges of the
/// suffix vertex; the (single) prefix->suffix edge is shared between both lists.
struct Splitter {
    to_splits: Vec<VId>,
    prefix: Vec<u8>,
    suffix: Vec<u8>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SplitNode {
    Prefix,
    Suffix,
    Middle(usize),
}

impl Splitter {
    fn new(sg: &SeqGraph, to_splits: &[VId]) -> Self {
        assert!(to_splits.len() > 1);
        let seqs: Vec<&[u8]> = to_splits.iter().map(|&v| sg.seq(v)).collect();
        let min = seqs.iter().map(|s| s.len()).min().unwrap();
        let prefix_len = common_maximum_prefix_length(&seqs);
        let suffix_len = common_maximum_suffix_length(&seqs, min - prefix_len);
        let k0 = seqs[0];
        Splitter {
            to_splits: to_splits.to_vec(),
            prefix: k0[..prefix_len].to_vec(),
            suffix: k0[k0.len() - suffix_len..].to_vec(),
        }
    }

    fn split_and_update(&self, sg: &mut SeqGraph, top: Option<VId>, bot: Option<VId>) {
        // ---- split() ----
        // split-graph edges: (src, tgt, edge data)
        let mut sedges: Vec<(SplitNode, SplitNode, SeqEdge)> = Vec::new();
        let mut new_middles: Vec<Vec<u8>> = Vec::new();
        let mut edges_to_remove: Vec<EId> = Vec::new();
        let mut ps_edge: Option<usize> = None;
        for &mid in &self.to_splits {
            let to_mid = match sg.g.incoming_edge_of(mid) {
                None => SeqEdge { mult: 0, is_ref: sg.g.is_reference_node(mid) },
                Some(e) => {
                    edges_to_remove.push(e);
                    *sg.g.e(e)
                }
            };
            let from_mid = match sg.g.outgoing_edge_of(mid) {
                None => SeqEdge { mult: 0, is_ref: sg.g.is_reference_node(mid) },
                Some(e) => {
                    edges_to_remove.push(e);
                    *sg.g.e(e)
                }
            };
            let mseq = sg.seq(mid);
            let start = self.prefix.len();
            let length = mseq.len() as i64 - self.suffix.len() as i64 - self.prefix.len() as i64;
            if length > 0 {
                let rem = mseq[start..start + length as usize].to_vec();
                let idx = new_middles.len();
                new_middles.push(rem);
                sedges.push((SplitNode::Prefix, SplitNode::Middle(idx), to_mid));
                sedges.push((SplitNode::Middle(idx), SplitNode::Suffix, from_mid));
            } else {
                let mut e = to_mid;
                e.add(&from_mid);
                match ps_edge {
                    Some(i) => {
                        let d = &mut sedges[i].2;
                        d.add(&e);
                    }
                    None => {
                        ps_edge = Some(sedges.len());
                        sedges.push((SplitNode::Prefix, SplitNode::Suffix, e));
                    }
                }
            }
        }
        // ---- updateGraph(top, bot) ----
        for &m in &self.to_splits {
            sg.g.remove_vertex(m);
        }
        for e in edges_to_remove {
            sg.g.remove_edge(e);
        }
        let mut mid_ids: Vec<VId> = Vec::with_capacity(new_middles.len());
        for m in new_middles {
            mid_ids.push(sg.g.add_vertex(m));
        }
        let prefix_out: Vec<usize> = (0..sedges.len()).filter(|&i| sedges[i].0 == SplitNode::Prefix).collect();
        let suffix_in: Vec<usize> = (0..sedges.len()).filter(|&i| sedges[i].1 == SplitNode::Suffix).collect();
        let has_ps = ps_edge.is_some();
        if has_ps {
            crate::prof::hit("splitter: prefix->suffix edge");
        }
        let has_only_ps = has_ps && prefix_out.len() == 1;
        let need_prefix = !self.prefix.is_empty() || (top.is_none() && !has_only_ps);
        let need_suffix = !self.suffix.is_empty() || (bot.is_none() && !has_only_ps);
        let mut prefix_id: Option<VId> = None;
        let mut suffix_id: Option<VId> = None;
        if need_prefix {
            let p = sg.g.add_vertex(self.prefix.clone());
            prefix_id = Some(p);
            if let Some(t) = top {
                let any_ref = prefix_out.iter().any(|&i| sedges[i].2.is_ref);
                sg.g.add_edge(t, p, SeqEdge { mult: 1, is_ref: any_ref });
            }
        }
        if need_suffix {
            let s = sg.g.add_vertex(self.suffix.clone());
            suffix_id = Some(s);
            if let Some(b) = bot {
                let any_ref = suffix_in.iter().any(|&i| sedges[i].2.is_ref);
                sg.g.add_edge(s, b, SeqEdge { mult: 1, is_ref: any_ref });
            }
        }
        if !need_prefix {
            crate::prof::hit("splitter: no prefix node");
        }
        if !need_suffix {
            crate::prof::hit("splitter: no suffix node");
        }
        let top_for_connect = if need_prefix { prefix_id } else { top };
        let bot_for_connect = if need_suffix { suffix_id } else { bot };
        let resolve = |n: SplitNode| -> Option<VId> {
            match n {
                SplitNode::Prefix => prefix_id,
                SplitNode::Suffix => suffix_id,
                SplitNode::Middle(i) => Some(mid_ids[i]),
            }
        };
        let mut ps_added = false;
        if let Some(tc) = top_for_connect {
            for &i in &prefix_out {
                let (_, t, d) = sedges[i];
                if t == SplitNode::Suffix {
                    if let Some(bc) = bot_for_connect {
                        if sg.g.add_edge(tc, bc, d).is_some() {
                            ps_added = true;
                        }
                    }
                } else {
                    sg.g.add_edge(tc, resolve(t).unwrap(), d);
                }
            }
        }
        if let Some(bc) = bot_for_connect {
            for &i in &suffix_in {
                let (s, _, d) = sedges[i];
                if s == SplitNode::Prefix {
                    // the shared prefix->suffix edge object: jgrapht returns false if already in graph
                    if ps_added {
                        continue;
                    }
                    let src = resolve(s).expect("prefix vertex not in graph (GATK would throw)");
                    sg.g.add_edge(src, bc, d);
                } else {
                    sg.g.add_edge(resolve(s).unwrap(), bc, d);
                }
            }
        }
    }
}

/// Structural snapshot for BaseGraph.graphEquals
struct GraphSnapshot {
    vseqs: Vec<Vec<u8>>,
    edges: Vec<(Vec<u8>, Vec<u8>)>,
}

impl GraphSnapshot {
    fn of(sg: &SeqGraph) -> Self {
        let vseqs = sg.g.vertex_iter().map(|v| sg.seq(v).to_vec()).collect();
        let edges = sg
            .g
            .edges()
            .into_iter()
            .map(|e| (sg.seq(sg.g.src(e)).to_vec(), sg.seq(sg.g.tgt(e)).to_vec()))
            .collect();
        GraphSnapshot { vseqs, edges }
    }

    fn graph_equals(&self, o: &GraphSnapshot) -> bool {
        if self.vseqs.len() != o.vseqs.len() || self.edges.len() != o.edges.len() {
            return false;
        }
        use std::collections::HashSet;
        let v2: HashSet<&Vec<u8>> = o.vseqs.iter().collect();
        if !self.vseqs.iter().all(|s| v2.contains(s)) {
            return false;
        }
        let e2: HashSet<&(Vec<u8>, Vec<u8>)> = o.edges.iter().collect();
        if !self.edges.iter().all(|e| e2.contains(e)) {
            return false;
        }
        let e1: HashSet<&(Vec<u8>, Vec<u8>)> = self.edges.iter().collect();
        o.edges.iter().all(|e| e1.contains(e))
    }
}
