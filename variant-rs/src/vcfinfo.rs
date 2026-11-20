//! The INFO column as htsjdk reads and writes it, shared by the tools that re-encode records (VariantFiltration,
//! MakeSitesOnlyVcf, ApplyVQSR, MergeVcfs).
//! * reading (`AbstractVCFCodec.parseInfo`): `KEY=value` is a String (an empty value is '.'), `KEY=a,b` a list, a bare
//!   key a flag (or '.' when the header declares it non-flag); `FLAG=0` is dropped for a key the header declares a
//!   flag; a repeated key keeps its last value.
//! * writing (`VCFEncoder`): keys sorted, `KEY=value`, a flag as `KEY`, '.' when there are none; a key the header does
//!   not declare is an error ("We require all VCFs to have complete VCF headers by default").

use std::collections::{BTreeMap, HashMap};

#[derive(Clone, Debug, PartialEq)]
pub enum Info {
    Str(String),
    List(Vec<String>),
    Flag,
}

/// INFO ID -> whether the header declares it a Flag, from `##INFO` lines.
pub fn declared(meta: &[String]) -> HashMap<String, bool> {
    let mut m = HashMap::new();
    for l in meta {
        if let Some(r) = l.strip_prefix("##INFO=<ID=") {
            let id = r.split([',', '>']).next().unwrap_or("").to_string();
            m.entry(id).or_insert(l.contains(",Type=Flag,") || l.contains(",Type=Flag>"));
        }
    }
    m
}

pub fn parse<'a>(col: &'a str, declared: &HashMap<String, bool>) -> HashMap<&'a str, Info> {
    let mut m = HashMap::new();
    if col == "." {
        return m;
    }
    for f in col.split(';') {
        match f.split_once('=') {
            Some((k, v)) => {
                if v.contains(',') {
                    m.insert(k, Info::List(v.split(',').map(|s| s.to_string()).collect()));
                } else if v == "0" && declared.get(k) == Some(&true) {
                    continue;
                } else {
                    m.insert(k, Info::Str(if v.is_empty() { ".".into() } else { v.to_string() }));
                }
            }
            None => {
                let v = if declared.get(f) == Some(&false) { Info::Str(".".into()) } else { Info::Flag };
                m.insert(f, v);
            }
        }
    }
    m
}

/// Writes the INFO column; Err(message) for a key the header does not declare (`at` = "contig:pos").
pub fn write(info: &HashMap<&str, Info>, declared: &HashMap<String, bool>, at: &str, out: &mut Vec<u8>) -> Result<(), String> {
    if info.is_empty() {
        out.push(b'.');
        return Ok(());
    }
    let sorted: BTreeMap<&str, &Info> = info.iter().map(|(k, v)| (*k, v)).collect();
    if let Some(k) = sorted.keys().find(|k| !declared.contains_key(**k)) {
        return Err(format!(
            "Key {k} found in VariantContext field INFO at {at} but this key isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers by default."
        ));
    }
    for (n, (k, v)) in sorted.iter().enumerate() {
        if n > 0 {
            out.push(b';');
        }
        out.extend_from_slice(k.as_bytes());
        match v {
            Info::Str(s) => {
                out.push(b'=');
                out.extend_from_slice(s.as_bytes());
            }
            Info::List(l) => {
                out.push(b'=');
                out.extend_from_slice(l.join(",").as_bytes());
            }
            Info::Flag => {}
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info() {
        let meta: Vec<String> = ["##INFO=<ID=DB,Number=0,Type=Flag,Description=\"d\">", "##INFO=<ID=QD,Number=1,Type=Float,Description=\"q\">", "##INFO=<ID=AF,Number=A,Type=Float,Description=\"a\">", "##INFO=<ID=X,Number=1,Type=String,Description=\"x\">"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let d = declared(&meta);
        let i = parse("QD=;DB=0;AF=0.1,0.2;X;QD=", &d);
        let mut out = Vec::new();
        write(&i, &d, "c:1", &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "AF=0.1,0.2;QD=.;X=.");
        let i = parse("DB;Y=1", &d);
        assert!(write(&i, &d, "c:1", &mut Vec::new()).is_err());
    }
}
