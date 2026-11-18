#!/usr/bin/env python3
"""gen.py: turns `gatk <Tool> --help` listings (GATK 4.7, saved as <Tool>.help.txt) into <Tool>.tsv tables that
variant-rs embeds: long name, short name (or -), value type, default value (or - when required / none listed),
whether the argument may be given several times, and its help section (required, optional, advanced, or
conditional: valid only with a given read filter or annotation). Regenerate with:
    for t in HaplotypeCaller BaseRecalibrator ApplyBQSR MarkDuplicates Mutect2 FilterMutectCalls GetPileupSummaries \\
             CalculateContamination LearnReadOrientationModel MergeMutectStats GatherPileupSummaries; do
        java -jar gatk.jar $t --help > $t.help.txt 2>&1; done; python3 gen.py *.help.txt
"""
import re
import sys

HEAD = re.compile(r"^(--[A-Za-z0-9_-]+)(?:,(-[A-Za-z0-9_-]+))?\s+<([^>]+)>(.*)$")

for path in sys.argv[1:]:
    tool = path.split("/")[-1].removesuffix(".help.txt")
    entries = []
    cur = None
    section = "optional"
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("Required Arguments"):
            section = "required"
        elif line.startswith("Optional Arguments"):
            section = "optional"
        elif line.startswith("Advanced Arguments"):
            section = "advanced"
        elif line.startswith("Conditional Arguments"):
            section = "conditional"
        m = HEAD.match(line)
        if m:
            cur = {"long": m.group(1), "short": m.group(2) or "-", "type": m.group(3), "text": m.group(4), "section": section}
            entries.append(cur)
        elif cur is not None and line.startswith(" "):
            cur["text"] += " " + line.strip()
        elif not line.strip():
            continue
        else:
            cur = None
    seen = set()
    with open(path.removesuffix(".help.txt") + ".tsv", "w") as out:
        out.write(f"# {tool}: GATK 4.7 arguments from `--help` (gen.py). long\tshort\ttype\tdefault\tmulti\tsection\n")
        for e in entries:
            if e["long"] in seen:
                continue
            seen.add(e["long"])
            text = " ".join(e["text"].split())
            d = re.search(r"Default value: (.*?)\.(?: |$)", text + " ")
            default = d.group(1).strip() if d else "-"
            multi = "1" if ("may be specified 0 or more times" in text or "at least once" in text) else "0"
            out.write(f"{e['long']}\t{e['short']}\t{e['type']}\t{default}\t{multi}\t{e['section']}\n")
    print(tool, len(seen), "arguments")
