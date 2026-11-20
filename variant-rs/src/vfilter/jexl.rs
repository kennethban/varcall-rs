//! The part of JEXL 2.1.1 (commons-jexl, as htsjdk 5 runs it for GATK's VariantFiltration) that filter expressions
//! use: number, string, `true` / `false` / `null` literals, variables, `!` / `not`, unary `-`, the comparisons
//! (`< <= > >= == !=` and `lt le gt ge eq ne`), `&&` / `and`, `||` / `or` and parentheses. Anything else (arithmetic,
//! method calls such as `vc.isSNP()`, regular expressions, functions) is refused when the expression is parsed.
//!
//! Evaluation follows `Interpreter` and `JexlArithmetic` with the engine htsjdk configures (strict, silent):
//! * an unknown variable, or an arithmetic error (comparing a list, a null operand), aborts the whole expression,
//!   which then counts as "missing" (`EvalErr::Missing`): a mismatch, or a match with
//!   `--missing-values-evaluate-as-failing`. `||` and `&&` short-circuit, so whether a missing value is reached at
//!   all depends on the order of the terms.
//! * comparisons: a Float / Double operand compares as doubles (a string parsed with `Double.parseDouble`, a Float
//!   literal through its `Float.toString`); otherwise an integer operand compares as longs (a string with
//!   `Long.parseLong`); otherwise as strings. A string that does not parse is a `NumberFormatException`, which
//!   GATK does not catch: the run fails (`EvalErr::Fatal`), as does a result that is not a boolean.
//! * decimal literals are Floats (`2.0`), `d` makes them Doubles; integer literals are Integers (or Longs).

use std::cmp::Ordering;

#[derive(Clone, Debug, PartialEq)]
pub enum Val {
    Null,
    Bool(bool),
    /// Integer or Long (both `isNumberable`)
    Int(i64),
    /// a Float literal, held as `toDouble` sees it (`Double.parseDouble(Float.toString(f))`)
    Float(f64),
    Double(f64),
    Str(String),
    /// a multi-valued INFO attribute (`List<String>`)
    List(Vec<String>),
}

#[derive(Debug, PartialEq)]
pub enum EvalErr {
    /// a JexlException: the expression's value is "missing"
    Missing,
    /// an exception GATK does not catch: the run fails with this message
    Fatal(String),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Cmp {
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
    Ne,
}

#[derive(Debug)]
pub enum Expr {
    Lit(Val),
    Var(String),
    Not(Box<Expr>),
    Neg(Box<Expr>),
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
    Cmp(Cmp, Box<Expr>, Box<Expr>),
}

/// The variables an expression reads.
pub trait Context {
    fn get(&self, name: &str) -> Option<Val>;
}

// ---------------------------------------------------------------- parsing

#[derive(Debug, PartialEq, Clone)]
enum Tok {
    LParen,
    RParen,
    And,
    Or,
    Not,
    Minus,
    Cmp(Cmp),
    Num(Val),
    Str(String),
    Ident(String),
    True,
    False,
    Null,
}

fn tokenize(s: &str) -> Result<Vec<Tok>, String> {
    let b = s.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    let unsupported = |what: &str| Err(format!("{what} in a filter expression is not supported (\"{s}\")"));
    while i < b.len() {
        let c = b[i];
        let two = if i + 1 < b.len() { &s[i..i + 2] } else { "" };
        if matches!(c, b' ' | b'\t' | b'\n' | b'\r' | 0x0c) {
            i += 1;
            continue;
        }
        let (tok, len) = match two {
            "&&" => (Tok::And, 2),
            "||" => (Tok::Or, 2),
            "==" => (Tok::Cmp(Cmp::Eq), 2),
            "!=" => (Tok::Cmp(Cmp::Ne), 2),
            "<=" => (Tok::Cmp(Cmp::Le), 2),
            ">=" => (Tok::Cmp(Cmp::Ge), 2),
            "=~" | "!~" => return unsupported("a regular expression match"),
            "?:" => return unsupported("?:"),
            "##" | "//" | "/*" => return unsupported("a comment"),
            _ => match c {
                b'(' => (Tok::LParen, 1),
                b')' => (Tok::RParen, 1),
                b'<' => (Tok::Cmp(Cmp::Lt), 1),
                b'>' => (Tok::Cmp(Cmp::Gt), 1),
                b'!' => (Tok::Not, 1),
                b'-' => (Tok::Minus, 1),
                b'"' | b'\'' => {
                    let end = s[i + 1..].find(c as char).map(|k| i + 1 + k).ok_or(format!("unterminated string in \"{s}\""))?;
                    let lit = &s[i + 1..end];
                    if lit.contains('\\') || lit.contains('\n') || lit.contains('\r') {
                        return unsupported("an escape in a string");
                    }
                    (Tok::Str(lit.to_string()), end + 1 - i)
                }
                b'0'..=b'9' => {
                    let (v, len) = number(&s[i..]).map_err(|e| format!("{e} (\"{s}\")"))?;
                    (Tok::Num(v), len)
                }
                b'a'..=b'z' | b'A'..=b'Z' | b'_' | b'$' | b'@' => {
                    let len = b[i..].iter().take_while(|&&c| c.is_ascii_alphanumeric() || matches!(c, b'_' | b'$' | b'@')).count();
                    let w = &s[i..i + len];
                    let t = match w {
                        "and" => Tok::And,
                        "or" => Tok::Or,
                        "not" => Tok::Not,
                        "eq" => Tok::Cmp(Cmp::Eq),
                        "ne" => Tok::Cmp(Cmp::Ne),
                        "lt" => Tok::Cmp(Cmp::Lt),
                        "le" => Tok::Cmp(Cmp::Le),
                        "gt" => Tok::Cmp(Cmp::Gt),
                        "ge" => Tok::Cmp(Cmp::Ge),
                        "true" => Tok::True,
                        "false" => Tok::False,
                        "null" => Tok::Null,
                        "if" | "else" | "for" | "foreach" | "while" | "new" | "var" | "empty" | "size" | "return" | "in" | "div" | "mod" => {
                            return unsupported(&format!("'{w}'"))
                        }
                        _ => Tok::Ident(w.to_string()),
                    };
                    (t, len)
                }
                _ => return unsupported(&format!("'{}'", c as char)),
            },
        };
        out.push(tok);
        i += len;
    }
    Ok(out)
}

/// A JEXL number literal at the start of `s`: (value, length). INTEGER_LITERAL is decimal (octal and hex are
/// refused), optionally `l`/`L`; FLOAT_LITERAL is `digits.digits[e[+-]digits]`, optionally `f`/`F` or `d`/`D`.
fn number(s: &str) -> Result<(Val, usize), String> {
    let b = s.as_bytes();
    let digits = |from: usize| b[from..].iter().take_while(|c| c.is_ascii_digit()).count();
    let int_len = digits(0);
    let mut len = int_len;
    let mut real = false;
    if len + 1 < b.len() && b[len] == b'.' && b[len + 1].is_ascii_digit() {
        real = true;
        len += 1 + digits(len + 1);
        if len < b.len() && (b[len] == b'e' || b[len] == b'E') {
            let mut k = len + 1;
            if k < b.len() && (b[k] == b'+' || b[k] == b'-') {
                k += 1;
            }
            let d = digits(k);
            if d > 0 {
                len = k + d;
            }
        }
    }
    let text = &s[..len];
    let suffix = b.get(len).copied().filter(|c| c.is_ascii_alphabetic());
    if real {
        let sig = text.split(['e', 'E']).next().unwrap().trim_start_matches(['0', '.']).replace('.', "").trim_end_matches('0').len();
        match suffix {
            Some(b'd') | Some(b'D') => Ok((Val::Double(text.parse().unwrap()), len + 1)),
            Some(b'f') | Some(b'F') | None => {
                // Float.valueOf, then toDouble's Double.parseDouble(Float.toString(f)): with at most 6 significant
                // digits the float's shortest decimal is the literal's own, so the double is the literal's
                if sig > 6 {
                    return Err(format!("the decimal literal {text} has more than 6 significant digits (a Java float); not supported"));
                }
                let f: f32 = text.parse().unwrap();
                if !f.is_finite() || (f != 0.0 && f.abs() < f32::MIN_POSITIVE) {
                    return Err(format!("the decimal literal {text} is out of float range; not supported"));
                }
                Ok((Val::Float(text.parse().unwrap()), len + suffix.map_or(0, |_| 1)))
            }
            Some(c) => Err(format!("the literal suffix '{}' is not supported", c as char)),
        }
    } else {
        if int_len > 1 && b[0] == b'0' {
            return Err(format!("the octal / hex literal {} is not supported", &s[..int_len + suffix.map_or(0, |_| 1)]));
        }
        match suffix {
            None | Some(b'l') | Some(b'L') => {
                let v: i64 = text.parse().map_err(|_| format!("the integer literal {text} does not fit a long; not supported"))?;
                Ok((Val::Int(v), len + suffix.map_or(0, |_| 1)))
            }
            Some(c) => Err(format!("the literal suffix '{}' is not supported", c as char)),
        }
    }
}

struct Parser {
    toks: Vec<Tok>,
    pos: usize,
    text: String,
}

impl Parser {
    fn peek(&self) -> Option<&Tok> {
        self.toks.get(self.pos)
    }
    fn next(&mut self) -> Option<Tok> {
        let t = self.toks.get(self.pos).cloned();
        self.pos += 1;
        t
    }
    fn err<T>(&self, what: &str) -> Result<T, String> {
        Err(format!("cannot parse the filter expression \"{}\": {what}", self.text))
    }
    fn or(&mut self) -> Result<Expr, String> {
        let mut l = self.and()?;
        while self.peek() == Some(&Tok::Or) {
            self.pos += 1;
            l = Expr::Or(Box::new(l), Box::new(self.and()?));
        }
        Ok(l)
    }
    fn and(&mut self) -> Result<Expr, String> {
        let mut l = self.equality()?;
        while self.peek() == Some(&Tok::And) {
            self.pos += 1;
            l = Expr::And(Box::new(l), Box::new(self.equality()?));
        }
        Ok(l)
    }
    // EqualityExpression: Relational ((== | !=) Relational)?
    fn equality(&mut self) -> Result<Expr, String> {
        let l = self.relational()?;
        if let Some(Tok::Cmp(op @ (Cmp::Eq | Cmp::Ne))) = self.peek().cloned() {
            self.pos += 1;
            return Ok(Expr::Cmp(op, Box::new(l), Box::new(self.relational()?)));
        }
        Ok(l)
    }
    // RelationalExpression: Unary ((< | > | <= | >=) Unary)?
    fn relational(&mut self) -> Result<Expr, String> {
        let l = self.unary()?;
        if let Some(Tok::Cmp(op @ (Cmp::Lt | Cmp::Le | Cmp::Gt | Cmp::Ge))) = self.peek().cloned() {
            self.pos += 1;
            return Ok(Expr::Cmp(op, Box::new(l), Box::new(self.unary()?)));
        }
        if self.peek() == Some(&Tok::Minus) {
            return self.err("arithmetic is not supported");
        }
        Ok(l)
    }
    fn unary(&mut self) -> Result<Expr, String> {
        match self.next() {
            Some(Tok::Minus) => Ok(Expr::Neg(Box::new(self.unary()?))),
            Some(Tok::Not) => Ok(Expr::Not(Box::new(self.unary()?))),
            Some(Tok::LParen) => {
                let e = self.or()?;
                if self.next() != Some(Tok::RParen) {
                    return self.err("')' expected");
                }
                Ok(e)
            }
            Some(Tok::Num(v)) => Ok(Expr::Lit(v)),
            Some(Tok::Str(s)) => Ok(Expr::Lit(Val::Str(s))),
            Some(Tok::True) => Ok(Expr::Lit(Val::Bool(true))),
            Some(Tok::False) => Ok(Expr::Lit(Val::Bool(false))),
            Some(Tok::Null) => Ok(Expr::Lit(Val::Null)),
            Some(Tok::Ident(n)) => {
                if matches!(n.as_str(), "vc" | "g" | "ALLELES" | "homRefCount" | "hetCount" | "homVarCount") {
                    return self.err(&format!("the variable {n} is not supported"));
                }
                Ok(Expr::Var(n))
            }
            _ => self.err("a value expected"),
        }
    }
}

/// Parses an expression as `JexlEngine.createExpression` would, refusing what is not supported.
pub fn parse(text: &str) -> Result<Expr, String> {
    let mut p = Parser { toks: tokenize(text)?, pos: 0, text: text.to_string() };
    if p.toks.is_empty() {
        return p.err("empty");
    }
    let e = p.or()?;
    if p.pos != p.toks.len() {
        return p.err("unexpected text after the expression");
    }
    Ok(e)
}

// ---------------------------------------------------------------- evaluation

fn missing<T>() -> Result<T, EvalErr> {
    Err(EvalErr::Missing)
}

fn nfe<T>(s: &str) -> Result<T, EvalErr> {
    Err(EvalErr::Fatal(format!("java.lang.NumberFormatException: For input string: \"{s}\"")))
}

/// `Double.parseDouble`: Java's floating-point literal grammar after `trim()` (decimal, `NaN`, `Infinity`, an
/// optional `f`/`d` suffix; hexadecimal is refused as a number format error here).
fn java_parse_double(s: &str) -> Result<f64, EvalErr> {
    let t = s.trim_matches(|c: char| c <= ' ');
    let (sign, body) = match t.as_bytes().first() {
        Some(b'+') => (1.0, &t[1..]),
        Some(b'-') => (-1.0, &t[1..]),
        _ => (1.0, t),
    };
    if body == "NaN" {
        return Ok(f64::NAN);
    }
    if body == "Infinity" {
        return Ok(sign * f64::INFINITY);
    }
    let num = body.strip_suffix(['f', 'F', 'd', 'D']).unwrap_or(body);
    let b = num.as_bytes();
    let mut i = 0;
    let int = b.iter().take_while(|c| c.is_ascii_digit()).count();
    i += int;
    let mut frac = 0;
    if i < b.len() && b[i] == b'.' {
        i += 1;
        frac = b[i..].iter().take_while(|c| c.is_ascii_digit()).count();
        i += frac;
    }
    if int + frac == 0 {
        return nfe(s);
    }
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        i += 1;
        if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
            i += 1;
        }
        let e = b[i..].iter().take_while(|c| c.is_ascii_digit()).count();
        if e == 0 {
            return nfe(s);
        }
        i += e;
    }
    if i != b.len() {
        return nfe(s);
    }
    Ok(sign * num.parse::<f64>().map_err(|_| EvalErr::Fatal(format!("java.lang.NumberFormatException: For input string: \"{s}\"")))?)
}

/// `Long.parseLong` (no trimming), "" as 0 (`toLong`).
fn java_parse_long(s: &str) -> Result<i64, EvalErr> {
    if s.is_empty() {
        return Ok(0);
    }
    let body = s.strip_prefix(['+', '-']).unwrap_or(s);
    if body.is_empty() || !body.bytes().all(|c| c.is_ascii_digit()) {
        return nfe(s);
    }
    s.parse::<i64>().or_else(|_| nfe(s))
}

fn to_double(v: &Val) -> Result<f64, EvalErr> {
    match v {
        Val::Double(d) | Val::Float(d) => Ok(*d),
        Val::Int(i) => Ok(*i as f64),
        Val::Bool(b) => Ok(if *b { 1.0 } else { 0.0 }),
        Val::Str(s) => {
            if s.trim_matches(|c: char| c <= ' ').is_empty() {
                Ok(f64::NAN)
            } else {
                java_parse_double(s)
            }
        }
        Val::Null | Val::List(_) => missing(),
    }
}

fn to_long(v: &Val) -> Result<i64, EvalErr> {
    match v {
        Val::Int(i) => Ok(*i),
        Val::Bool(b) => Ok(*b as i64),
        Val::Str(s) => java_parse_long(s),
        // only reached with no Float / Double operand
        Val::Double(_) | Val::Float(_) | Val::Null | Val::List(_) => missing(),
    }
}

fn to_string(v: &Val) -> String {
    match v {
        Val::Str(s) => s.clone(),
        Val::Bool(b) => b.to_string(),
        Val::Int(i) => i.to_string(),
        Val::List(l) => format!("[{}]", l.join(", ")),
        Val::Null | Val::Double(_) | Val::Float(_) => unreachable!("compared as numbers"),
    }
}

/// `toBoolean`: null is an error (strict); numbers are true unless 0 or NaN; strings unless "" or "false".
fn to_boolean(v: &Val) -> Result<bool, EvalErr> {
    Ok(match v {
        Val::Null => return missing(),
        Val::Bool(b) => *b,
        Val::Int(i) => *i != 0,
        Val::Float(d) | Val::Double(d) => !d.is_nan() && *d != 0.0,
        Val::Str(s) => !s.is_empty() && s != "false",
        Val::List(_) => false,
    })
}

/// `JexlArithmetic.compare` (both operands non-null)
fn compare(l: &Val, r: &Val, eq: bool) -> Result<Ordering, EvalErr> {
    let fp = |v: &Val| matches!(v, Val::Float(_) | Val::Double(_));
    if fp(l) || fp(r) {
        let (a, b) = (to_double(l)?, to_double(r)?);
        return Ok(match (a.is_nan(), b.is_nan()) {
            (true, true) => Ordering::Equal,
            (true, false) => Ordering::Less,
            (false, true) => Ordering::Greater,
            _ => a.partial_cmp(&b).unwrap(),
        });
    }
    if matches!(l, Val::Int(_)) || matches!(r, Val::Int(_)) {
        return Ok(to_long(l)?.cmp(&to_long(r)?));
    }
    if matches!(l, Val::Str(_)) || matches!(r, Val::Str(_)) {
        return Ok(java_compare_strings(&to_string(l), &to_string(r)));
    }
    match (l, r) {
        _ if eq => Ok(if l == r { Ordering::Equal } else { Ordering::Less }),
        (Val::Bool(a), Val::Bool(b)) => Ok(a.cmp(b)),
        (Val::Bool(_), _) | (_, Val::Bool(_)) => Err(EvalErr::Fatal("java.lang.ClassCastException in a filter expression comparison".into())),
        _ => missing(),
    }
}

/// `String.compareTo`: by UTF-16 code units.
fn java_compare_strings(a: &str, b: &str) -> Ordering {
    a.encode_utf16().cmp(b.encode_utf16())
}

fn equals(l: &Val, r: &Val) -> Result<bool, EvalErr> {
    match (l, r) {
        (Val::Null, Val::Null) => Ok(true),
        (Val::Null, _) | (_, Val::Null) => Ok(false),
        (Val::Bool(_), _) | (_, Val::Bool(_)) => Ok(to_boolean(l)? == to_boolean(r)?),
        _ => Ok(compare(l, r, true)? == Ordering::Equal),
    }
}

fn eval(e: &Expr, ctx: &dyn Context) -> Result<Val, EvalErr> {
    match e {
        Expr::Lit(v) => Ok(v.clone()),
        Expr::Var(n) => ctx.get(n).ok_or(EvalErr::Missing),
        Expr::Not(a) => {
            // toBoolean outside a try: a null here is an ArithmeticException GATK does not catch
            let v = eval(a, ctx)?;
            match to_boolean(&v) {
                Ok(b) => Ok(Val::Bool(!b)),
                Err(_) => Err(EvalErr::Fatal("java.lang.ArithmeticException: jexl.null".into())),
            }
        }
        Expr::Neg(a) => match eval(a, ctx)? {
            Val::Int(i) => Ok(Val::Int(-i)),
            Val::Float(d) => Ok(Val::Float(-d)),
            Val::Double(d) => Ok(Val::Double(-d)),
            Val::Bool(b) => Ok(Val::Bool(!b)),
            _ => missing(),
        },
        Expr::Or(a, b) => {
            if to_boolean(&eval(a, ctx)?)? {
                return Ok(Val::Bool(true));
            }
            Ok(Val::Bool(to_boolean(&eval(b, ctx)?)?))
        }
        Expr::And(a, b) => {
            if !to_boolean(&eval(a, ctx)?)? {
                return Ok(Val::Bool(false));
            }
            Ok(Val::Bool(to_boolean(&eval(b, ctx)?)?))
        }
        Expr::Cmp(op, a, b) => {
            let (l, r) = (eval(a, ctx)?, eval(b, ctx)?);
            let null = matches!(l, Val::Null) || matches!(r, Val::Null);
            let v = match op {
                Cmp::Eq => equals(&l, &r)?,
                Cmp::Ne => !equals(&l, &r)?,
                _ if null => false,
                Cmp::Lt => compare(&l, &r, false)? == Ordering::Less,
                Cmp::Gt => compare(&l, &r, false)? == Ordering::Greater,
                Cmp::Le => compare(&l, &r, false)? != Ordering::Greater,
                Cmp::Ge => compare(&l, &r, false)? != Ordering::Less,
            };
            Ok(Val::Bool(v))
        }
    }
}

/// The expression's boolean value (`(Boolean) exp.evaluate(context)`), or why it has none.
pub fn evaluate(e: &Expr, ctx: &dyn Context) -> Result<bool, EvalErr> {
    match eval(e, ctx)? {
        Val::Bool(b) => Ok(b),
        v => Err(EvalErr::Fatal(format!("java.lang.ClassCastException: the filter expression's value {v:?} is not a Boolean"))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    struct M(HashMap<&'static str, Val>);
    impl Context for M {
        fn get(&self, n: &str) -> Option<Val> {
            self.0.get(n).cloned()
        }
    }

    fn ev(e: &str, vars: &[(&'static str, Val)]) -> Result<bool, EvalErr> {
        evaluate(&parse(e).unwrap(), &M(vars.iter().cloned().collect()))
    }
    fn s(x: &str) -> Val {
        Val::Str(x.into())
    }

    #[test]
    fn jexl() {
        assert_eq!(ev("QD < 2.0", &[("QD", s("1.5"))]), Ok(true));
        assert_eq!(ev("ExcessHet > 54.69", &[("ExcessHet", s("54.69"))]), Ok(false));
        assert_eq!(ev("ExcessHet > 54.69", &[("ExcessHet", s("54.6901"))]), Ok(true));
        // a missing variable aborts the whole expression, unless short-circuited away
        assert_eq!(ev("QD < 2.0 || FS > 30.0", &[("FS", s("40"))]), Err(EvalErr::Missing));
        assert_eq!(ev("FS > 30.0 || QD < 2.0", &[("FS", s("40"))]), Ok(true));
        // an integer literal compares as longs: a decimal string is a NumberFormatException
        assert!(matches!(ev("QD > 2", &[("QD", s("10.18"))]), Err(EvalErr::Fatal(_))));
        assert_eq!(ev("DP > 10", &[("DP", s("25"))]), Ok(true));
        // a list cannot be compared
        assert_eq!(ev("AF > 0.5", &[("AF", Val::List(vec!["0.1".into(), "0.7".into()]))]), Err(EvalErr::Missing));
        // strings, word operators, negative literals, not
        assert_eq!(ev("TYPE == \"SNP\" and not (MQRankSum lt -12.5)", &[("TYPE", s("SNP")), ("MQRankSum", s("-1"))]), Ok(true));
        assert_eq!(ev("CHROM == 'chr20' && POS > 10", &[("CHROM", s("chr20")), ("POS", Val::Int(11))]), Ok(true));
        // a flag is a Boolean
        assert_eq!(ev("DB && QD < 5.0", &[("DB", Val::Bool(true)), ("QD", s("4"))]), Ok(true));
        // a non-boolean result fails the run
        assert!(matches!(ev("FILTER", &[("FILTER", s("0"))]), Err(EvalErr::Fatal(_))));
        // QUAL is a Double
        assert_eq!(ev("QUAL < 10.4139", &[("QUAL", Val::Double(10.41))]), Ok(true));
        assert!(parse("vc.isSNP()").is_err());
        assert!(parse("QD * 2 < 3").is_err());
        assert!(parse("QD - 1 < 3").is_err());
        assert!(parse("AF =~ 'x'").is_err());
        assert!(parse("QD < 0.1234567").is_err());
        assert!(parse("QD < 010").is_err());
    }
}
