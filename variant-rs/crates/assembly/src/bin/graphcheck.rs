//! Compare the Rust SeqGraph cleanup/simplification + k-best finder against GATK results recorded
//! by scratch/GraphProbe.java, a development probe not in the repository.   usage: graphcheck <probe.jsonl>
use hc_assembly::seqgraph::{SeqEdge, SeqGraph};
use hc_assembly::{cleanup_seq_graph, kbest, AssemblyStatus};
use serde_json::Value;
use std::io::{BufRead, BufReader};

fn main() {
    let path = std::env::args().nth(1).expect("usage: graphcheck <probe.jsonl>");
    let f = BufReader::new(std::fs::File::open(path).unwrap());
    std::panic::set_hook(Box::new(|_| {}));
    let (mut n, mut bad, mut errs_both, mut status_jar) = (0, 0, 0, 0);
    for (li, line) in f.lines().enumerate() {
        let o: Value = serde_json::from_str(&line.unwrap()).unwrap();
        n += 1;
        let mut sg = SeqGraph::new(25);
        let vs: Vec<u32> = o["vertices"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| sg.g.add_vertex(v.as_str().unwrap().as_bytes().to_vec()))
            .collect();
        for e in o["edges"].as_array().unwrap() {
            let e: Vec<i64> = e.as_array().unwrap().iter().map(|x| x.as_i64().unwrap()).collect();
            sg.g.add_edge(vs[e[0] as usize], vs[e[1] as usize], SeqEdge { mult: e[2] as i32, is_ref: e[3] == 1 });
        }
        let res = std::panic::catch_unwind(move || {
            let mut sg = sg;
            sg.clean_non_ref_paths();
            let (status, sg) = cleanup_seq_graph(sg);
            let mut fv: Vec<String> = Vec::new();
            let mut idx = vec![usize::MAX; sg.g.vs.len()];
            for v in sg.g.vertex_iter() {
                idx[v as usize] = fv.len();
                fv.push(String::from_utf8_lossy(sg.seq(v)).into_owned());
            }
            let fe: Vec<(usize, usize, i32, i32)> = sg
                .g
                .edges()
                .into_iter()
                .map(|e| (idx[sg.g.src(e) as usize], idx[sg.g.tgt(e) as usize], sg.g.e(e).mult, sg.g.e(e).is_ref as i32))
                .collect();
            let mut kb: Vec<(String, f64)> = Vec::new();
            if status == AssemblyStatus::AssembledSomeVariation {
                let s = sg.g.reference_source_vertex().unwrap();
                let t = sg.g.reference_sink_vertex().unwrap();
                for p in kbest::find_best_haplotypes(&sg, s, t, 128) {
                    kb.push((String::from_utf8_lossy(&p.bases).into_owned(), p.score));
                }
            }
            (status, fv, fe, kb)
        });
        match (res, o.get("error")) {
            (Err(_), Some(_)) => {
                errs_both += 1;
                continue;
            }
            (Err(e), None) => {
                bad += 1;
                let msg = e.downcast_ref::<String>().cloned().or(e.downcast_ref::<&str>().map(|s| s.to_string()));
                println!("graph {}: Rust failed ({:?}), GATK did not", li, msg);
                continue;
            }
            (Ok(_), Some(ge)) => {
                bad += 1;
                println!("graph {}: GATK failed ({}), Rust did not", li, ge);
                continue;
            }
            (Ok((status, fv, fe, kb)), None) => {
                let st = match status {
                    AssemblyStatus::AssembledSomeVariation => "ASSEMBLED_SOME_VARIATION",
                    AssemblyStatus::JustAssembledReference => "JUST_ASSEMBLED_REFERENCE",
                    AssemblyStatus::Failed => "FAILED",
                };
                if st == "JUST_ASSEMBLED_REFERENCE" {
                    status_jar += 1;
                }
                let gfv: Vec<String> = o["final_vertices"].as_array().unwrap().iter().map(|v| v.as_str().unwrap().to_string()).collect();
                let gfe: Vec<(usize, usize, i32, i32)> = o["final_edges"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|e| {
                        let e = e.as_array().unwrap();
                        (e[0].as_u64().unwrap() as usize, e[1].as_u64().unwrap() as usize, e[2].as_i64().unwrap() as i32, e[3].as_i64().unwrap() as i32)
                    })
                    .collect();
                let gkb: Vec<(String, f64)> = o["kbest"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|x| {
                        let x = x.as_array().unwrap();
                        (x[0].as_str().unwrap().to_string(), x[1].as_str().unwrap().parse::<f64>().unwrap())
                    })
                    .collect();
                let kb_ok = gkb.len() == kb.len()
                    && gkb.iter().zip(kb.iter()).all(|(a, b)| a.0 == b.0 && (a.1.to_bits() == b.1.to_bits() || (a.1.is_nan() && b.1.is_nan())));
                if st != o["status"].as_str().unwrap() || gfv != fv || gfe != fe || !kb_ok {
                    bad += 1;
                    if bad <= 5 {
                        println!(
                            "graph {} mismatch: status {} vs {}; vertices eq {}; edges eq {}; kbest eq {} ({} vs {})",
                            li,
                            o["status"],
                            st,
                            gfv == fv,
                            gfe == fe,
                            kb_ok,
                            gkb.len(),
                            kb.len()
                        );
                    }
                }
            }
        }
    }
    println!("graphs: {} mismatches: {} (both threw: {}, JUST_ASSEMBLED_REFERENCE: {})", n, bad, errs_both, status_jar);
    print!("{}", hc_assembly::prof::hits_report());
}
