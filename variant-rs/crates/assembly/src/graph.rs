//! Arena-based directed graph that reproduces the iteration-order semantics of
//! jgrapht 1.1.0's DefaultDirectedGraph as used by GATK's BaseGraph:
//!  * vertexSet()/edgeSet() iterate in insertion order (LinkedHashMap);
//!  * incoming/outgoing edge sets are insertion-ordered lists (ArrayUnenforcedSet);
//!  * no multi-edges: addEdge(u,v,e) is a no-op returning false if an edge u->v exists;
//!  * removal preserves the relative order of the remaining elements.
//! Vertex and edge ids are assigned at insertion, so id order == insertion order.

pub type VId = u32;
pub type EId = u32;
pub const NONE: u32 = u32::MAX;

pub trait EdgeLike {
    fn is_ref(&self) -> bool;
    fn multiplicity(&self) -> i32;
}

#[derive(Clone)]
pub struct VSlot<V> {
    pub data: V,
    pub alive: bool,
    pub out: Vec<EId>,
    pub inc: Vec<EId>,
}

#[derive(Clone)]
pub struct ESlot<E> {
    pub data: E,
    pub src: VId,
    pub tgt: VId,
    pub alive: bool,
}

#[derive(Clone)]
pub struct Graph<V, E> {
    pub vs: Vec<VSlot<V>>,
    pub es: Vec<ESlot<E>>,
    n_alive_v: usize,
    n_alive_e: usize,
}

impl<V, E> Default for Graph<V, E> {
    fn default() -> Self {
        Graph { vs: Vec::new(), es: Vec::new(), n_alive_v: 0, n_alive_e: 0 }
    }
}

impl<V, E> Graph<V, E> {
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn vertex_count(&self) -> usize {
        self.n_alive_v
    }
    #[inline]
    pub fn edge_count(&self) -> usize {
        self.n_alive_e
    }

    pub fn add_vertex(&mut self, data: V) -> VId {
        let id = self.vs.len() as VId;
        self.vs.push(VSlot { data, alive: true, out: Vec::new(), inc: Vec::new() });
        self.n_alive_v += 1;
        id
    }

    #[inline]
    pub fn contains_vertex(&self, v: VId) -> bool {
        (v as usize) < self.vs.len() && self.vs[v as usize].alive
    }

    #[inline]
    pub fn v(&self, v: VId) -> &V {
        &self.vs[v as usize].data
    }
    #[inline]
    pub fn v_mut(&mut self, v: VId) -> &mut V {
        &mut self.vs[v as usize].data
    }
    #[inline]
    pub fn e(&self, e: EId) -> &E {
        &self.es[e as usize].data
    }
    #[inline]
    pub fn e_mut(&mut self, e: EId) -> &mut E {
        &mut self.es[e as usize].data
    }
    #[inline]
    pub fn src(&self, e: EId) -> VId {
        self.es[e as usize].src
    }
    #[inline]
    pub fn tgt(&self, e: EId) -> VId {
        self.es[e as usize].tgt
    }
    #[inline]
    pub fn out_edges(&self, v: VId) -> &[EId] {
        &self.vs[v as usize].out
    }
    #[inline]
    pub fn in_edges(&self, v: VId) -> &[EId] {
        &self.vs[v as usize].inc
    }
    #[inline]
    pub fn out_degree(&self, v: VId) -> usize {
        self.vs[v as usize].out.len()
    }
    #[inline]
    pub fn in_degree(&self, v: VId) -> usize {
        self.vs[v as usize].inc.len()
    }
    #[inline]
    pub fn edge_alive(&self, e: EId) -> bool {
        self.es[e as usize].alive
    }

    /// Alive vertices in insertion order (snapshot).
    pub fn vertices(&self) -> Vec<VId> {
        let mut r = Vec::with_capacity(self.n_alive_v);
        for (i, s) in self.vs.iter().enumerate() {
            if s.alive {
                r.push(i as VId);
            }
        }
        r
    }

    /// Iterator over alive vertex ids in insertion order (no snapshot).
    pub fn vertex_iter(&self) -> impl Iterator<Item = VId> + '_ {
        self.vs.iter().enumerate().filter(|(_, s)| s.alive).map(|(i, _)| i as VId)
    }

    pub fn edges(&self) -> Vec<EId> {
        let mut r = Vec::with_capacity(self.n_alive_e);
        for (i, s) in self.es.iter().enumerate() {
            if s.alive {
                r.push(i as EId);
            }
        }
        r
    }

    pub fn get_edge(&self, u: VId, v: VId) -> Option<EId> {
        self.vs[u as usize].out.iter().copied().find(|&e| self.es[e as usize].tgt == v)
    }

    /// jgrapht addEdge(u, v, e) for a fresh edge object: returns None (and does nothing)
    /// if an edge u->v already exists. Panics if a vertex is not in the graph.
    pub fn add_edge(&mut self, u: VId, v: VId, data: E) -> Option<EId> {
        assert!(self.contains_vertex(u) && self.contains_vertex(v), "no such vertex in graph");
        if self.get_edge(u, v).is_some() {
            return None;
        }
        let id = self.es.len() as EId;
        self.es.push(ESlot { data, src: u, tgt: v, alive: true });
        self.vs[u as usize].out.push(id);
        self.vs[v as usize].inc.push(id);
        self.n_alive_e += 1;
        Some(id)
    }

    pub fn remove_edge(&mut self, e: EId) -> bool {
        let slot = &mut self.es[e as usize];
        if !slot.alive {
            return false;
        }
        slot.alive = false;
        let (u, v) = (slot.src, slot.tgt);
        let out = &mut self.vs[u as usize].out;
        if let Some(p) = out.iter().position(|&x| x == e) {
            out.remove(p);
        }
        let inc = &mut self.vs[v as usize].inc;
        if let Some(p) = inc.iter().position(|&x| x == e) {
            inc.remove(p);
        }
        self.n_alive_e -= 1;
        true
    }

    pub fn remove_edge_between(&mut self, u: VId, v: VId) -> bool {
        match self.get_edge(u, v) {
            Some(e) => self.remove_edge(e),
            None => false,
        }
    }

    pub fn remove_vertex(&mut self, v: VId) -> bool {
        if !self.contains_vertex(v) {
            return false;
        }
        let inc = std::mem::take(&mut self.vs[v as usize].inc);
        for e in inc {
            // reinsert so remove_edge can clean it up consistently
            let slot = &mut self.es[e as usize];
            if slot.alive {
                slot.alive = false;
                self.n_alive_e -= 1;
                let u = slot.src;
                let out = &mut self.vs[u as usize].out;
                if let Some(p) = out.iter().position(|&x| x == e) {
                    out.remove(p);
                }
            }
        }
        let out = std::mem::take(&mut self.vs[v as usize].out);
        for e in out {
            let slot = &mut self.es[e as usize];
            if slot.alive {
                slot.alive = false;
                self.n_alive_e -= 1;
                let w = slot.tgt;
                let inc = &mut self.vs[w as usize].inc;
                if let Some(p) = inc.iter().position(|&x| x == e) {
                    inc.remove(p);
                }
            }
        }
        self.vs[v as usize].alive = false;
        self.n_alive_v -= 1;
        true
    }

    /// LinkedHashSet of targets in outgoing edge order (no multi-edges, so no duplicates
    /// except via self loops which are still unique targets).
    pub fn outgoing_vertices(&self, v: VId) -> Vec<VId> {
        let mut r: Vec<VId> = Vec::with_capacity(self.out_degree(v));
        for &e in self.out_edges(v) {
            let t = self.tgt(e);
            if !r.contains(&t) {
                r.push(t);
            }
        }
        r
    }
    pub fn incoming_vertices(&self, v: VId) -> Vec<VId> {
        let mut r: Vec<VId> = Vec::with_capacity(self.in_degree(v));
        for &e in self.in_edges(v) {
            let s = self.src(e);
            if !r.contains(&s) {
                r.push(s);
            }
        }
        r
    }

    #[inline]
    pub fn is_source(&self, v: VId) -> bool {
        self.in_degree(v) == 0
    }
    #[inline]
    pub fn is_sink(&self, v: VId) -> bool {
        self.out_degree(v) == 0
    }

    /// BaseGraph.incomingEdgeOf (single edge or None)
    pub fn incoming_edge_of(&self, v: VId) -> Option<EId> {
        let inc = self.in_edges(v);
        assert!(inc.len() <= 1, "Cannot get a single incoming edge for a vertex with multiple incoming edges");
        inc.first().copied()
    }
    pub fn outgoing_edge_of(&self, v: VId) -> Option<EId> {
        let out = self.out_edges(v);
        assert!(out.len() <= 1, "Cannot get a single outgoing edge for a vertex with multiple outgoing edges");
        out.first().copied()
    }

    /// CycleDetector.detectCycles()
    pub fn has_cycles(&self) -> bool {
        // iterative DFS with colors
        let n = self.vs.len();
        let mut color = vec![0u8; n];
        let mut stack: Vec<(VId, usize)> = Vec::new();
        for s in 0..n {
            if !self.vs[s].alive || color[s] != 0 {
                continue;
            }
            color[s] = 1;
            stack.push((s as VId, 0));
            while let Some(&mut (v, ref mut idx)) = stack.last_mut() {
                let out = &self.vs[v as usize].out;
                if *idx < out.len() {
                    let t = self.es[out[*idx] as usize].tgt as usize;
                    *idx += 1;
                    if color[t] == 1 {
                        return true;
                    }
                    if color[t] == 0 {
                        color[t] = 1;
                        stack.push((t as VId, 0));
                    }
                } else {
                    color[v as usize] = 2;
                    stack.pop();
                }
            }
        }
        false
    }

    /// Reachable vertices (BaseGraphIterator) from start following the given directions.
    pub fn reachable(&self, start: VId, follow_in: bool, follow_out: bool) -> Vec<bool> {
        let mut visited = vec![false; self.vs.len()];
        let mut todo = vec![start];
        while let Some(v) = todo.pop() {
            if visited[v as usize] {
                continue;
            }
            visited[v as usize] = true;
            if follow_in {
                for &e in self.in_edges(v) {
                    todo.push(self.src(e));
                }
            }
            if follow_out {
                for &e in self.out_edges(v) {
                    todo.push(self.tgt(e));
                }
            }
        }
        visited
    }
}

impl<V, E: EdgeLike> Graph<V, E> {
    pub fn is_reference_node(&self, v: VId) -> bool {
        if self.in_edges(v).iter().any(|&e| self.e(e).is_ref()) || self.out_edges(v).iter().any(|&e| self.e(e).is_ref())
        {
            return true;
        }
        self.vertex_count() == 1
    }

    pub fn is_ref_source(&self, v: VId) -> bool {
        if self.in_edges(v).iter().any(|&e| self.e(e).is_ref()) {
            return false;
        }
        if self.out_edges(v).iter().any(|&e| self.e(e).is_ref()) {
            return true;
        }
        self.vertex_count() == 1
    }

    pub fn is_ref_sink(&self, v: VId) -> bool {
        if self.out_edges(v).iter().any(|&e| self.e(e).is_ref()) {
            return false;
        }
        if self.in_edges(v).iter().any(|&e| self.e(e).is_ref()) {
            return true;
        }
        self.vertex_count() == 1
    }

    pub fn reference_source_vertex(&self) -> Option<VId> {
        self.vertex_iter().find(|&v| self.is_ref_source(v))
    }
    pub fn reference_sink_vertex(&self) -> Option<VId> {
        self.vertex_iter().find(|&v| self.is_ref_sink(v))
    }

    /// BaseGraph.getNextReferenceVertex(v, allowNonRefPaths, blacklistedEdge)
    pub fn next_reference_vertex(&self, v: Option<VId>, allow_non_ref: bool, blacklist: Option<EId>) -> Option<VId> {
        let v = v?;
        let out = self.out_edges(v);
        if out.is_empty() {
            return None;
        }
        for &e in out {
            if self.e(e).is_ref() {
                return Some(self.tgt(e));
            }
        }
        if !allow_non_ref {
            return None;
        }
        let mut cnt = 0;
        let mut first = NONE;
        for &e in out {
            if Some(e) == blacklist {
                continue;
            }
            if cnt == 0 {
                first = e;
            }
            cnt += 1;
            if cnt >= 2 {
                break;
            }
        }
        if cnt == 1 {
            Some(self.tgt(first))
        } else {
            None
        }
    }

    /// BaseGraph.getPrevReferenceVertex
    pub fn prev_reference_vertex(&self, v: VId) -> Option<VId> {
        self.in_edges(v).iter().map(|&e| self.src(e)).find(|&s| self.is_reference_node(s))
    }
}
