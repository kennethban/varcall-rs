#!/usr/bin/env python3
"""Translate HotSpot's x86_64 libm stubs (JDK 17 macroAssembler_x86_{exp,log,log10,pow}.cpp, the
`#ifdef _LP64` half) into Rust, one statement per assembler instruction.

    for f in exp log log10 pow; do curl -sSfLO \
      https://raw.githubusercontent.com/openjdk/jdk17u/master/src/hotspot/cpu/x86/macroAssembler_x86_$f.cpp; done
    python3 gen_jmath.py <dir with the four .cpp files> > ../../crates/assembly/src/jmath.rs

The temp registers' physical names come from stubGenerator_x86_64.cpp (generate_libm*): they alias the
registers the stubs also name directly (pow: tmp4 = r11).

Each SSE2 / integer instruction becomes a call to a small helper with the instruction's exact
semantics (see PRELUDE); labels become arms of a `loop { match pc { .. } }` state machine, and the
flags register is modelled so that every jcc sees what the hardware would. The generator refuses any
instruction or operand form it does not know.
"""
import re
import sys

FUNCS = [  # (rust name, cpp file, MacroAssembler method, temp registers as stubGenerator_x86_64.cpp passes them)
    ("exp", "macroAssembler_x86_exp.cpp", "fast_exp", {"tmp": "r11"}),
    ("log", "macroAssembler_x86_log.cpp", "fast_log", {"tmp1": "r11", "tmp2": "r8"}),
    ("log10", "macroAssembler_x86_log10.cpp", "fast_log10", {}),
    ("pow", "macroAssembler_x86_pow.cpp", "fast_pow", {"tmp1": "r8", "tmp2": "r9", "tmp3": "r10", "tmp4": "r11"}),
]

GPR_ALIAS = {"eax": "rax", "rax": "rax", "ecx": "rcx", "rcx": "rcx", "edx": "rdx", "rdx": "rdx"}
GPR_EXTRA = {"r8", "r9", "r10", "r11"}
COND = {
    "equal": "fl.zf",
    "notEqual": "!fl.zf",
    "above": "(!fl.cf && !fl.zf)",
    "aboveEqual": "!fl.cf",
    "below": "fl.cf",
    "belowEqual": "(fl.cf || fl.zf)",
    "less": "(fl.sf != fl.of)",
    "greaterEqual": "(fl.sf == fl.of)",
    "greater": "(!fl.zf && fl.sf == fl.of)",
    "lessEqual": "(fl.zf || fl.sf != fl.of)",
    "negative": "fl.sf",
    "positive": "!fl.sf",
}


def die(msg):
    sys.exit("gen_jmath: " + msg)


def lp64_part(text):
    a = text.index("#ifdef _LP64")
    b = text.index("\n#else", a)
    return text[a:b]


def parse_tables(part):
    tabs = {}
    for m in re.finditer(r"juint\s+(_\w+)\[\]\s*=\s*\{(.*?)\};", part, re.S):
        vals = [int(v, 16) for v in re.findall(r"0x([0-9a-fA-F]+)UL", m.group(2))]
        tabs[m.group(1)] = vals
    return tabs


def split_args(s):
    out, depth, cur = [], 0, ""
    for ch in s:
        if ch == "," and depth == 0:
            out.append(cur.strip())
            cur = ""
            continue
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        cur += ch
    if cur.strip():
        out.append(cur.strip())
    return out


class Fn:
    def __init__(self, name, part, method, regs):
        self.name = name
        self.regs = regs  # the stub's Register parameters -> physical registers (they alias, e.g. tmp4 = r11)
        self.tables = parse_tables(part)
        body = part[part.index("void MacroAssembler::" + method):]
        self.alias = dict(re.findall(r"address\s+(\w+)\s*=\s*\(address\)(_\w+);", body))
        # one memory image per function: all its tables concatenated (byte offsets)
        self.base, off, self.image = {}, 0, []
        for t, vals in self.tables.items():
            self.base[t] = off
            self.image += vals
            off += 4 * len(vals)
        self.gprs = set()
        start = body.index("bind(start);")
        self.lines = self.clean(body[start:])

    def clean(self, code):
        out, skip_depth, in_else = [], 0, False
        lines = code.split("\n")
        i = 0
        while i < len(lines):
            raw = lines[i]
            i += 1
            line = raw.split("//")[0].strip()
            if not line:
                continue
            if line.startswith("if (VM_Version::supports_sse3())"):
                # take the SSE3 branch (both branches compute the same bits), drop the else
                while True:
                    l2 = lines[i].split("//")[0].strip()
                    i += 1
                    if l2 == "}":
                        break
                    out.append(l2)
                if lines[i].strip() == "else {":
                    i += 1
                    while lines[i].strip() != "}":
                        i += 1
                    i += 1
                continue
            if line == "}":
                break  # end of method
            out.append(line)
        return out

    # ---- operands -------------------------------------------------------------------------------
    def gpr(self, s):
        s = self.regs.get(s, s)
        if s in GPR_ALIAS:
            v = GPR_ALIAS[s]
        elif s in GPR_EXTRA:
            v = s
        else:
            return None
        self.gprs.add(v)
        return v

    def xmm(self, s):
        m = re.fullmatch(r"xmm([0-7])", s)
        return "x" + m.group(1) if m else None

    def imm(self, s):
        if s == "INT_MIN":
            return -2147483648
        if re.fullmatch(r"-?\d+", s):
            return int(s)
        if re.fullmatch(r"0x[0-9a-fA-F]+", s):
            return int(s, 16)
        return None

    def mem(self, s):
        """('stack', off) or ('mem', rust expr for the byte address)"""
        m = re.fullmatch(r"ExternalAddress\((?:(\d+)\s*\+\s*)?(\w+)\)", s)
        if m:
            off = int(m.group(1) or 0)
            t = self.alias.get(m.group(2))
            if t is None:
                die(f"{self.name}: unknown table {m.group(2)}")
            return ("mem", str(self.base[t] + off))
        m = re.fullmatch(r"Address\((.*)\)", s)
        if not m:
            return None
        a = split_args(m.group(1))
        if a[0] == "rsp":
            if len(a) != 2:
                die(f"bad stack address {s}")
            return ("stack", str(int(a[1])))
        base = self.gpr(a[0])
        if base is None:
            die(f"bad base in {s}")
        if len(a) == 2 and self.imm(a[1]) is not None:
            return ("mem", f"{base}.wrapping_add({self.imm(a[1])}i64 as u64)")
        idx = self.gpr(a[1])
        if idx is None:
            die(f"bad index in {s}")
        scale, disp = 1, 0
        if len(a) >= 3:
            scale = {"Address::times_1": 1, "Address::times_2": 2, "Address::times_4": 4, "Address::times_8": 8}[a[2]]
        if len(a) == 4:
            disp = int(a[3])
        return ("mem", f"{base}.wrapping_add({idx}.wrapping_mul({scale})).wrapping_add({disp}i64 as u64)")

    def load(self, m, bits):
        kind, e = m
        if kind == "stack":
            return f"sld{bits}(&st, {e})"
        return f"ld{bits}(&{self.memname()}, {e})"

    def memname(self):
        return "MEM_" + self.name.upper()

    # ---- instructions ---------------------------------------------------------------------------
    PACKED = {"addpd", "subpd", "mulpd", "andpd", "xorpd", "pand", "por", "pxor", "paddd", "psubd",
              "paddq", "unpcklpd", "unpckhpd", "movlhps"}
    SCALAR = {"addsd", "subsd", "mulsd", "divsd", "sqrtsd", "mulss"}
    ALU32 = {"addl": "add32", "subl": "sub32", "andl": "and32", "orl": "or32", "xorl": "xor32"}

    def insn(self, op, a):
        X, G, I, M = self.xmm, self.gpr, self.imm, self.mem
        n = len(a)
        if op in self.PACKED and n == 2 and X(a[0]):
            if X(a[1]):
                return f"{X(a[0])} = {op}({X(a[0])}, {X(a[1])});"
            if M(a[1]):
                return f"{X(a[0])} = {op}({X(a[0])}, {self.load(M(a[1]), 128)});"
        if op in self.SCALAR and n == 2 and X(a[0]):
            if X(a[1]):
                return f"{X(a[0])} = {op}({X(a[0])}, {X(a[1])});"
            if M(a[1]) and op != "mulss":
                return f"{X(a[0])} = {op}({X(a[0])}, [{self.load(M(a[1]), 64)}, 0]);"
        if op in ("psllq", "psrlq", "psrld", "psrad") and n == 2 and X(a[0]):
            if I(a[1]) is not None:
                return f"{X(a[0])} = {op}({X(a[0])}, {I(a[1])});"
            if op == "psllq" and X(a[1]):
                return f"{X(a[0])} = {op}({X(a[0])}, {X(a[1])}[0]);"
        if op == "pshufd" and n == 3 and X(a[0]) and X(a[1]) and I(a[2]) is not None:
            return f"{X(a[0])} = pshufd({X(a[1])}, {I(a[2])});"
        if op in ("rcpps", "movddup") and n == 2 and X(a[0]) and X(a[1]):
            return f"{X(a[0])} = {op}({X(a[1])});"
        if op in ("movdqu", "movapd") and n == 2 and X(a[0]):
            if X(a[1]):
                return f"{X(a[0])} = {X(a[1])};"
            if M(a[1]):
                return f"{X(a[0])} = {self.load(M(a[1]), 128)};"
        if op in ("movq", "movsd") and n == 2:
            if X(a[0]) and M(a[1]):
                return f"{X(a[0])} = [{self.load(M(a[1]), 64)}, 0];"
            if M(a[0]) and X(a[1]) and M(a[0])[0] == "stack":
                return f"sst64(&mut st, {M(a[0])[1]}, {X(a[1])}[0]);"
        if op == "movdl" and n == 2:
            if X(a[0]) and G(a[1]):
                return f"{X(a[0])} = [{G(a[1])} & 0xffff_ffff, 0];"
            if G(a[0]) and X(a[1]):
                return f"{G(a[0])} = {X(a[1])}[0] & 0xffff_ffff;"
        if op == "movdq" and n == 2:
            if X(a[0]) and G(a[1]):
                return f"{X(a[0])} = [{G(a[1])}, 0];"
            if G(a[0]) and X(a[1]):
                return f"{G(a[0])} = {X(a[1])}[0];"
        if op == "pinsrw" and n == 3 and X(a[0]) and G(a[1]) and I(a[2]) is not None:
            return f"{X(a[0])} = pinsrw({X(a[0])}, {G(a[1])}, {I(a[2])});"
        if op == "pextrw" and n == 3 and G(a[0]) and X(a[1]) and I(a[2]) is not None:
            return f"{G(a[0])} = pextrw({X(a[1])}, {I(a[2])});"
        if op == "cvtsi2sdl" and n == 2 and X(a[0]) and G(a[1]):
            return f"{X(a[0])} = cvtsi2sdl({X(a[0])}, {G(a[1])});"
        if op == "movl" and n == 2:
            if G(a[0]):
                if I(a[1]) is not None:
                    return f"{G(a[0])} = {I(a[1]) & 0xffffffff:#x};"
                if G(a[1]):
                    return f"{G(a[0])} = {G(a[1])} & 0xffff_ffff;"
                if M(a[1]):
                    return f"{G(a[0])} = {self.load(M(a[1]), 32)};"
            if M(a[0]) and M(a[0])[0] == "stack":
                if I(a[1]) is not None:
                    return f"sst32(&mut st, {M(a[0])[1]}, {I(a[1]) & 0xffffffff:#x});"
                if G(a[1]):
                    return f"sst32(&mut st, {M(a[0])[1]}, {G(a[1])});"
        if op == "movw" and n == 2 and G(a[0]) and M(a[1]):
            return f"{G(a[0])} = ({G(a[0])} & !0xffff) | {self.load(M(a[1]), 16)};"
        if op == "mov64" and n == 2 and G(a[0]) and I(a[1]) is not None:
            return f"{G(a[0])} = {I(a[1]) & (2**64 - 1):#x};"
        if op == "lea" and n == 2 and G(a[0]) and M(a[1]) and M(a[1])[0] == "mem":
            return f"{G(a[0])} = {M(a[1])[1]};"
        if op in self.ALU32 or op in ("cmpl", "testl"):
            f = self.ALU32.get(op) or {"cmpl": "sub32", "testl": "and32"}[op]
            if n == 2 and G(a[0]):
                if I(a[1]) is not None:
                    src = f"{I(a[1]) & 0xffffffff:#x}"
                elif G(a[1]):
                    src = G(a[1])
                elif M(a[1]):
                    src = self.load(M(a[1]), 32)
                else:
                    src = None
                if src is not None:
                    if op in ("cmpl", "testl"):
                        return f"{f}(&mut fl, {G(a[0])}, {src});"
                    return f"{G(a[0])} = {f}(&mut fl, {G(a[0])}, {src});"
        if op in ("addq", "subq") and a[0] == "rsp":
            return f"// {op}(rsp, {a[1]})"
        if op in ("shll", "shrl", "sarl", "shlq") and G(a[0]):
            f = {"shll": "shl32", "shrl": "shr32", "sarl": "sar32", "shlq": "shl64"}[op]
            if n == 2 and I(a[1]) is not None:
                return f"{G(a[0])} = {f}(&mut fl, {G(a[0])}, {I(a[1])});"
            if n == 1:
                self.gprs.add("rcx")
                return f"{G(a[0])} = {f}(&mut fl, {G(a[0])}, rcx as u32);"
        if op == "bsrl" and n == 2 and G(a[0]) and G(a[1]):
            return f"{G(a[0])} = bsr32(&mut fl, {G(a[0])}, {G(a[1])});"
        if op == "cmp64" and n == 2 and G(a[0]) and M(a[1]):
            return f"sub64(&mut fl, {G(a[0])}, {self.load(M(a[1]), 64)});"
        return None

    def translate(self):
        labels = []
        for l in self.lines:
            m = re.fullmatch(r"bind\((\w+)\);", l)
            if m:
                labels.append(m.group(1))
        lid = {l: i for i, l in enumerate(labels)}
        end = len(labels)
        blocks = {}
        cur = None
        body = []
        for l in self.lines:
            m = re.fullmatch(r"(\w+)\((.*)\);", l)
            if not m:
                die(f"{self.name}: cannot parse {l!r}")
            op, args = m.group(1), split_args(m.group(2))
            if op == "bind":
                if cur is not None:
                    blocks[cur] = body
                cur, body = lid[args[0]], []
                continue
            if op == "assert_different_registers":
                continue
            if op in ("jcc", "jccb"):
                c = args[0].replace("Assembler::", "")
                if c not in COND:
                    die(f"{self.name}: condition {c}")
                body.append(f"if {COND[c].strip('()')} {{ pc = {lid[args[1]]}; continue; }} // {l}")
                continue
            if op == "jmp":
                body.append(f"pc = {lid[args[0]]}; continue; // {l}")
                continue
            r = self.insn(op, args)
            if r is None:
                die(f"{self.name}: unsupported instruction {l!r}")
            body.append(r if r.startswith("//") else f"{r} // {l}")
        blocks[cur] = body
        out = []
        w = out.append
        w(f"static {self.memname()}: [u32; {len(self.image)}] = [")
        for t, vals in self.tables.items():
            w(f"    // {t} (byte offset {self.base[t]})")
            for i in range(0, len(vals), 6):
                w("    " + " ".join(f"{v:#010x}," for v in vals[i:i + 6]))
        w("];")
        w("")
        return out, labels, blocks, end

    def emit(self, doc):
        out, labels, blocks, end = self.translate()
        w = out.append
        args = "x: f64, y: f64" if self.name == "pow" else "x: f64"
        w(doc)
        w("#[allow(unused_assignments, unused_mut, unused_variables, unreachable_code, clippy::all)]")
        w(f"pub fn {self.name}({args}) -> f64 {{")
        w("    let mut x0: X = [x.to_bits(), 0];")
        w("    let mut x1: X = [" + ("y.to_bits()" if self.name == "pow" else "0") + ", 0];")
        w("    let (mut x2, mut x3, mut x4, mut x5, mut x6, mut x7): (X, X, X, X, X, X) = Default::default();")
        for g in sorted(self.gprs):
            w(f"    let mut {g}: u64 = 0;")
        w("    let mut st = [0u8; 48];")
        w("    let mut fl = Flags::default();")
        w("    let mut pc = 0usize;")
        w("    loop {")
        w("        match pc {")
        for i, lab in enumerate(labels):
            w(f"            {i} => {{ // {lab}")
            for s in blocks[i]:
                w("                " + s)
            w(f"                pc = {i + 1};")
            w("            }")
        w(f"            {end} => return f64::from_bits(x0[0]),")
        w("            _ => unreachable!(),")
        w("        }")
        w("    }")
        w("}")
        w("")
        return out


PRELUDE = r'''//! Bit-exact ports of HotSpot's x86_64 `java.lang.Math` intrinsics (JDK 17): `exp`, `log`, `log10`,
//! `pow`. GENERATED by `tests/jmath/gen_jmath.py` from openjdk/jdk17u
//! `src/hotspot/cpu/x86/macroAssembler_x86_{exp,log,log10,pow}.cpp` (Intel LIBM, `_LP64` code);
//! do not edit by hand.
//!
//! Each assembler instruction of the stub is one statement below, in the stub's order, operating on
//! modelled xmm registers (`[u64; 2]`, lane 0 = low double), 64-bit general registers, a stack frame
//! and the flags register; labels are arms of a `loop { match pc }` state machine. All double
//! arithmetic is plain IEEE add/sub/mul/div/sqrt (no FMA, no reassociation), so results are the
//! stub's bits. `rcpps` (an approximation whose bits are CPU-specific) runs the real instruction,
//! so the results match HotSpot on the same CPU, as Java does.
//! Verified against `java` with `tests/jmath/` (see run.sh there).
//!
//! `Math.log1p` / `Math.expm1` are not intrinsics (they are StrictMath, i.e. fdlibm): hand-written
//! ports in `fdlibm.rs`, re-exported here.

pub use crate::fdlibm::{expm1, log1p};

type X = [u64; 2];

#[derive(Default, Clone, Copy)]
struct Flags {
    cf: bool,
    zf: bool,
    sf: bool,
    of: bool,
}

#[inline(always)]
fn f(b: u64) -> f64 {
    f64::from_bits(b)
}

/// SSE result of a double operation: an invalid operation gives the x86 default NaN (negative
/// quiet NaN, which LLVM's constant folding would not), a NaN operand is returned quieted (first
/// source first), as the hardware does.
#[inline(always)]
fn sse(r: f64, a: u64, b: u64) -> u64 {
    if !r.is_nan() {
        r.to_bits()
    } else if f(a).is_nan() {
        a | 1 << 51
    } else if f(b).is_nan() {
        b | 1 << 51
    } else {
        0xfff8_0000_0000_0000
    }
}

// ---- packed double / integer SSE2 ----
#[inline(always)]
fn addpd(d: X, s: X) -> X {
    [sse(f(d[0]) + f(s[0]), d[0], s[0]), sse(f(d[1]) + f(s[1]), d[1], s[1])]
}
#[inline(always)]
fn subpd(d: X, s: X) -> X {
    [sse(f(d[0]) - f(s[0]), d[0], s[0]), sse(f(d[1]) - f(s[1]), d[1], s[1])]
}
#[inline(always)]
fn mulpd(d: X, s: X) -> X {
    [sse(f(d[0]) * f(s[0]), d[0], s[0]), sse(f(d[1]) * f(s[1]), d[1], s[1])]
}
#[inline(always)]
fn andpd(d: X, s: X) -> X {
    [d[0] & s[0], d[1] & s[1]]
}
#[inline(always)]
fn pand(d: X, s: X) -> X {
    andpd(d, s)
}
#[inline(always)]
fn por(d: X, s: X) -> X {
    [d[0] | s[0], d[1] | s[1]]
}
#[inline(always)]
fn xorpd(d: X, s: X) -> X {
    [d[0] ^ s[0], d[1] ^ s[1]]
}
#[inline(always)]
fn pxor(d: X, s: X) -> X {
    xorpd(d, s)
}
#[inline(always)]
fn dw(v: X) -> [u32; 4] {
    [v[0] as u32, (v[0] >> 32) as u32, v[1] as u32, (v[1] >> 32) as u32]
}
#[inline(always)]
fn from_dw(d: [u32; 4]) -> X {
    [d[0] as u64 | (d[1] as u64) << 32, d[2] as u64 | (d[3] as u64) << 32]
}
#[inline(always)]
fn paddd(d: X, s: X) -> X {
    let (a, b) = (dw(d), dw(s));
    from_dw([a[0].wrapping_add(b[0]), a[1].wrapping_add(b[1]), a[2].wrapping_add(b[2]), a[3].wrapping_add(b[3])])
}
#[inline(always)]
fn psubd(d: X, s: X) -> X {
    let (a, b) = (dw(d), dw(s));
    from_dw([a[0].wrapping_sub(b[0]), a[1].wrapping_sub(b[1]), a[2].wrapping_sub(b[2]), a[3].wrapping_sub(b[3])])
}
#[inline(always)]
fn paddq(d: X, s: X) -> X {
    [d[0].wrapping_add(s[0]), d[1].wrapping_add(s[1])]
}
#[inline(always)]
fn unpcklpd(d: X, s: X) -> X {
    [d[0], s[0]]
}
#[inline(always)]
fn unpckhpd(d: X, s: X) -> X {
    [d[1], s[1]]
}
#[allow(dead_code)] // only in the stubs' non-SSE3 alternative
#[inline(always)]
fn movlhps(d: X, s: X) -> X {
    [d[0], s[0]]
}
#[inline(always)]
fn movddup(s: X) -> X {
    [s[0], s[0]]
}
#[inline(always)]
fn psllq(d: X, n: u64) -> X {
    if n > 63 { [0, 0] } else { [d[0] << n, d[1] << n] }
}
#[inline(always)]
fn psrlq(d: X, n: u64) -> X {
    if n > 63 { [0, 0] } else { [d[0] >> n, d[1] >> n] }
}
#[inline(always)]
fn psrld(d: X, n: u32) -> X {
    let a = dw(d);
    if n > 31 { [0, 0] } else { from_dw([a[0] >> n, a[1] >> n, a[2] >> n, a[3] >> n]) }
}
#[inline(always)]
fn psrad(d: X, n: u32) -> X {
    let a = dw(d);
    let n = n.min(31);
    from_dw([
        ((a[0] as i32) >> n) as u32,
        ((a[1] as i32) >> n) as u32,
        ((a[2] as i32) >> n) as u32,
        ((a[3] as i32) >> n) as u32,
    ])
}
#[inline(always)]
fn pshufd(s: X, imm: u32) -> X {
    let a = dw(s);
    from_dw([
        a[(imm & 3) as usize],
        a[(imm >> 2 & 3) as usize],
        a[(imm >> 4 & 3) as usize],
        a[(imm >> 6 & 3) as usize],
    ])
}
#[inline(always)]
fn pinsrw(d: X, r: u64, i: u32) -> X {
    let (lane, sh) = ((i / 4) as usize, (i % 4) * 16);
    let mut v = d;
    v[lane] = (v[lane] & !(0xffffu64 << sh)) | ((r & 0xffff) << sh);
    v
}
#[inline(always)]
fn pextrw(s: X, i: u32) -> u64 {
    (s[(i / 4) as usize] >> ((i % 4) * 16)) & 0xffff
}
/// RCPPS: the CPU's own approximate reciprocal (its bits differ between vendors; HotSpot uses the
/// instruction, so do we).
#[inline(always)]
fn rcpps(s: X) -> X {
    #[cfg(target_arch = "x86_64")]
    unsafe {
        use std::arch::x86_64::*;
        let v = _mm_set_epi64x(s[1] as i64, s[0] as i64);
        let r = _mm_castps_si128(_mm_rcp_ps(_mm_castsi128_ps(v)));
        let mut o = [0u64; 2];
        _mm_storeu_si128(o.as_mut_ptr() as *mut __m128i, r);
        o
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let a = dw(s);
        let r = |x: u32| (1.0f32 / f32::from_bits(x)).to_bits();
        from_dw([r(a[0]), r(a[1]), r(a[2]), r(a[3])])
    }
}

// ---- scalar (low lane; high lane of the destination kept) ----
#[inline(always)]
fn addsd(d: X, s: X) -> X {
    [sse(f(d[0]) + f(s[0]), d[0], s[0]), d[1]]
}
#[inline(always)]
fn subsd(d: X, s: X) -> X {
    [sse(f(d[0]) - f(s[0]), d[0], s[0]), d[1]]
}
#[inline(always)]
fn mulsd(d: X, s: X) -> X {
    [sse(f(d[0]) * f(s[0]), d[0], s[0]), d[1]]
}
#[inline(always)]
fn divsd(d: X, s: X) -> X {
    [sse(f(d[0]) / f(s[0]), d[0], s[0]), d[1]]
}
#[inline(always)]
fn sqrtsd(d: X, s: X) -> X {
    [sse(f(s[0]).sqrt(), s[0], s[0]), d[1]]
}
#[inline(always)]
fn mulss(d: X, s: X) -> X {
    let (a, b) = (d[0] as u32, s[0] as u32);
    let p = f32::from_bits(a) * f32::from_bits(b);
    let r = if !p.is_nan() {
        p.to_bits()
    } else if f32::from_bits(a).is_nan() {
        a | 1 << 22
    } else if f32::from_bits(b).is_nan() {
        b | 1 << 22
    } else {
        0xffc0_0000
    };
    [(d[0] & !0xffff_ffff) | r as u64, d[1]]
}
#[inline(always)]
fn cvtsi2sdl(d: X, r: u64) -> X {
    [(r as u32 as i32 as f64).to_bits(), d[1]]
}

// ---- memory ----
#[inline(always)]
fn ld32(m: &[u32], a: u64) -> u64 {
    debug_assert!(a % 4 == 0);
    m[(a / 4) as usize] as u64
}
#[inline(always)]
fn ld64(m: &[u32], a: u64) -> u64 {
    ld32(m, a) | ld32(m, a + 4) << 32
}
#[inline(always)]
fn ld128(m: &[u32], a: u64) -> X {
    [ld64(m, a), ld64(m, a + 8)]
}
#[inline(always)]
fn sld16(st: &[u8; 48], o: usize) -> u64 {
    u16::from_le_bytes([st[o], st[o + 1]]) as u64
}
#[inline(always)]
fn sld32(st: &[u8; 48], o: usize) -> u64 {
    u32::from_le_bytes(st[o..o + 4].try_into().unwrap()) as u64
}
#[inline(always)]
fn sld64(st: &[u8; 48], o: usize) -> u64 {
    u64::from_le_bytes(st[o..o + 8].try_into().unwrap())
}
#[inline(always)]
fn sst32(st: &mut [u8; 48], o: usize, v: u64) {
    st[o..o + 4].copy_from_slice(&(v as u32).to_le_bytes());
}
#[inline(always)]
fn sst64(st: &mut [u8; 48], o: usize, v: u64) {
    st[o..o + 8].copy_from_slice(&v.to_le_bytes());
}

// ---- integer ALU (32-bit ops zero-extend into the 64-bit register) ----
#[inline(always)]
fn szf32(fl: &mut Flags, r: u32) {
    fl.zf = r == 0;
    fl.sf = (r as i32) < 0;
}
#[inline(always)]
fn add32(fl: &mut Flags, a: u64, b: u64) -> u64 {
    let (a, b) = (a as u32, b as u32);
    let (r, c) = a.overflowing_add(b);
    fl.cf = c;
    fl.of = (a as i32).overflowing_add(b as i32).1;
    szf32(fl, r);
    r as u64
}
#[inline(always)]
fn sub32(fl: &mut Flags, a: u64, b: u64) -> u64 {
    let (a, b) = (a as u32, b as u32);
    let (r, c) = a.overflowing_sub(b);
    fl.cf = c;
    fl.of = (a as i32).overflowing_sub(b as i32).1;
    szf32(fl, r);
    r as u64
}
#[inline(always)]
fn sub64(fl: &mut Flags, a: u64, b: u64) -> u64 {
    let (r, c) = a.overflowing_sub(b);
    fl.cf = c;
    fl.of = (a as i64).overflowing_sub(b as i64).1;
    fl.zf = r == 0;
    fl.sf = (r as i64) < 0;
    r
}
#[inline(always)]
fn logic32(fl: &mut Flags, r: u32) -> u64 {
    fl.cf = false;
    fl.of = false;
    szf32(fl, r);
    r as u64
}
#[inline(always)]
fn and32(fl: &mut Flags, a: u64, b: u64) -> u64 {
    logic32(fl, a as u32 & b as u32)
}
#[inline(always)]
fn or32(fl: &mut Flags, a: u64, b: u64) -> u64 {
    logic32(fl, a as u32 | b as u32)
}
#[inline(always)]
fn xor32(fl: &mut Flags, a: u64, b: u64) -> u64 {
    logic32(fl, a as u32 ^ b as u32)
}
#[inline(always)]
fn shl32(fl: &mut Flags, a: u64, n: u32) -> u64 {
    let (a, n) = (a as u32, n & 31);
    if n == 0 {
        return a as u64;
    }
    let r = a << n;
    fl.cf = (a >> (32 - n)) & 1 != 0;
    fl.of = ((r >> 31) != 0) != fl.cf;
    szf32(fl, r);
    r as u64
}
#[inline(always)]
fn shr32(fl: &mut Flags, a: u64, n: u32) -> u64 {
    let (a, n) = (a as u32, n & 31);
    if n == 0 {
        return a as u64;
    }
    let r = a >> n;
    fl.cf = (a >> (n - 1)) & 1 != 0;
    fl.of = (a >> 31) != 0;
    szf32(fl, r);
    r as u64
}
#[inline(always)]
fn sar32(fl: &mut Flags, a: u64, n: u32) -> u64 {
    let (a, n) = (a as u32 as i32, n & 31);
    if n == 0 {
        return a as u32 as u64;
    }
    let r = (a >> n) as u32;
    fl.cf = (a >> (n - 1)) & 1 != 0;
    fl.of = false;
    szf32(fl, r);
    r as u64
}
#[inline(always)]
fn shl64(fl: &mut Flags, a: u64, n: u32) -> u64 {
    let n = n & 63;
    if n == 0 {
        return a;
    }
    let r = a << n;
    fl.cf = (a >> (64 - n)) & 1 != 0;
    fl.of = ((r >> 63) != 0) != fl.cf;
    fl.zf = r == 0;
    fl.sf = (r as i64) < 0;
    r
}
/// BSR r32: zero source sets ZF and leaves the destination unchanged.
#[inline(always)]
fn bsr32(fl: &mut Flags, d: u64, s: u64) -> u64 {
    let s = s as u32;
    if s == 0 {
        fl.zf = true;
        d
    } else {
        fl.zf = false;
        (31 - s.leading_zeros()) as u64
    }
}
'''

DOCS = {
    "exp": "/// `java.lang.Math.exp` as HotSpot computes it on x86_64 (`StubRoutines::dexp`).",
    "log": "/// `java.lang.Math.log` as HotSpot computes it on x86_64 (`StubRoutines::dlog`).",
    "log10": "/// `java.lang.Math.log10` as HotSpot computes it on x86_64 (`StubRoutines::dlog10`).",
    "pow": "/// `java.lang.Math.pow` as HotSpot computes it on x86_64 (`StubRoutines::dpow`).",
}


def main():
    d = sys.argv[1]
    out = [PRELUDE]
    for name, fname, method, regs in FUNCS:
        part = lp64_part(open(f"{d}/{fname}").read())
        fn = Fn(name, part, method, regs)
        out.append("\n".join(fn.emit(DOCS[name])))
    sys.stdout.write("\n".join(out))


main()
