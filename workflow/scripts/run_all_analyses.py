#!/usr/bin/env python3
"""Functional Variant Impact Analysis - Kelton Guimaraes 2026"""

import os, sys, gzip, base64
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from datetime import date

VEP_FILE = "/home/kelton/microbial_variant_calling/results/annotation/SRR7801919_snps_vep.txt"
GFF_FILE = "/home/kelton/microbial_variant_calling/data/reference/reference.gff"
OUTDIR = "results"

for d in ["enrichment", "prioritization", "figures"]:
    os.makedirs(f"{OUTDIR}/{d}", exist_ok=True)

print("=" * 60)
print("  FUNCTIONAL VARIANT IMPACT ANALYSIS")
print("  Candida albicans SRR7801919")
print("  Analyst: Kelton Guimaraes")
print("=" * 60)

# Parse VEP
print("\nParsing VEP annotations...")
variants = []
gene_impacts = defaultdict(list)
high_impact = []
missense_variants = []
impact_counter = Counter()
consequence_counter = Counter()

with open(VEP_FILE) as f:
    for line in f:
        if line.startswith("##"):
            continue
        p = line.strip().split("\t")
        if len(p) < 14:
            continue
        cons = p[6] if len(p) > 6 else "."
        gene = p[3] if len(p) > 3 else "."
        
        consequence_counter[cons] += 1
        
        if any(x in cons for x in ["stop_gained", "frameshift", "splice_donor", "splice_acceptor"]):
            impact = "HIGH"
        elif "missense" in cons:
            impact = "MODERATE"
        else:
            impact = "LOW"
        
        impact_counter[impact] += 1
        variants.append({"gene": gene, "impact": impact, "cons": cons})
        
        if gene != ".":
            gene_impacts[gene].append(cons)
        
        if impact == "HIGH":
            score = 10 if "stop_gained" in cons or "frameshift" in cons else 8
            high_impact.append({"gene": gene, "cons": cons, "score": score})
        
        if "missense" in cons:
            missense_variants.append({"gene": gene, "cons": cons})

print(f"  Total variants: {len(variants):,}")
print(f"  Affected genes: {len(gene_impacts):,}")
print(f"  HIGH impact: {len(high_impact)}")
print(f"  MODERATE: {impact_counter['MODERATE']:,}")
print(f"  LOW: {impact_counter['LOW']:,}")
print(f"  Missense: {len(missense_variants):,}")

# Parse GFF
print("\nParsing GFF for gene descriptions...")
gene_desc = {}
with open(GFF_FILE) as f:
    for line in f:
        if line.startswith("#"):
            continue
        p = line.strip().split("\t")
        if len(p) >= 9:
            attrs = p[8]
            gname = ""
            prod = ""
            for attr in attrs.split(";"):
                if "=" in attr:
                    k, v = attr.split("=", 1)
                    if k == "Name": gname = v
                    elif k == "product": prod = v
            if gname and prod:
                gene_desc[gname] = prod

print(f"  Genes with descriptions: {len(gene_desc):,}")

# GO enrichment
print("\nRunning enrichment analysis...")
categories = {
    "Drug Response": ["drug", "resistance", "azole", "efflux", "transporter", "multidrug", "MDR", "CDR"],
    "Cell Wall": ["cell wall", "chitin", "glucan", "mannan", "mannoprotein"],
    "Biofilm & Adhesion": ["biofilm", "adhesion", "adhesin", "hyphal", "filamentous", "yeast-to-hypha"],
    "Stress Response": ["stress", "oxidative", "heat shock", "HSP", "osmotic", "DNA damage"],
    "Virulence": ["virulence", "pathogenesis", "protease", "lipase", "phospholipase"],
    "Metabolism": ["metabolism", "biosynthesis", "ergosterol", "lipid", "glycolysis"],
    "DNA Replication & Repair": ["DNA replication", "DNA repair", "recombination", "mismatch"],
    "Transcription": ["transcription", "RNA polymerase", "transcription factor", "zinc finger"],
}

affected_set = set(gene_impacts.keys())
total_bg = max(len(gene_desc), 1)
enrichment = {}

for cat, keywords in categories.items():
    bg_in_cat = set()
    for g, d in gene_desc.items():
        for kw in keywords:
            if kw.lower() in d.lower():
                bg_in_cat.add(g)
                break
    
    aff_in_cat = affected_set & bg_in_cat
    if len(aff_in_cat) >= 2:
        expected = (len(bg_in_cat) / total_bg) * len(affected_set)
        ratio = len(aff_in_cat) / max(expected, 1)
        enrichment[cat] = {"count": len(aff_in_cat), "ratio": round(ratio, 2), "genes": list(aff_in_cat)[:5]}

enrich_sorted = sorted(enrichment.items(), key=lambda x: x[1]["ratio"], reverse=True)
print(f"  Enriched categories: {len(enrich_sorted)}")
for cat, d in enrich_sorted[:5]:
    print(f"    {cat}: {d['count']} genes, {d['ratio']}x enriched")

# Top affected genes
top_genes = sorted(gene_impacts.items(), key=lambda x: len(x[1]), reverse=True)[:20]
print(f"\n  Top 5 affected genes:")
for g, impacts in top_genes[:5]:
    print(f"    {g}: {len(impacts)} variants")

# Save data files
with open(f"{OUTDIR}/enrichment/go_enrichment.txt", "w") as f:
    f.write(f"GO Enrichment Analysis - SRR7801919\n\n")
    for cat, d in enrich_sorted:
        f.write(f"{cat}\t{d['count']} genes\t{d['ratio']}x enriched\t{', '.join(d['genes'][:3])}\n")

with open(f"{OUTDIR}/prioritization/high_priority.txt", "w") as f:
    f.write("Rank\tGene\tConsequence\tScore\n")
    for i, v in enumerate(sorted(high_impact, key=lambda x: x["score"], reverse=True)[:20], 1):
        f.write(f"{i}\t{v['gene']}\t{v['cons']}\t{v['score']}\n")

# Generate figures
print("\nGenerating figures...")

fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle("Functional Variant Analysis — C. albicans SRR7801919", fontsize=16, fontweight="bold")

# Panel A: Impact pie
ax1 = axes[0, 0]
labels = list(impact_counter.keys())
sizes = list(impact_counter.values())
ax1.pie(sizes, labels=labels, autopct="%1.1f%%", colors=["#e74c3c", "#f39c12", "#3498db"], startangle=90)
ax1.set_title("A. Variant Impact Distribution", fontweight="bold")

# Panel B: Enrichment
ax2 = axes[0, 1]
if enrich_sorted:
    cats = [x[0] for x in enrich_sorted]
    ratios = [x[1]["ratio"] for x in enrich_sorted]
    colors_b = plt.cm.viridis(np.linspace(0.2, 0.9, len(cats)))
    ax2.barh(cats, ratios, color=colors_b, edgecolor="black")
    ax2.axvline(x=1.0, color="red", linestyle="--", label="No enrichment")
    ax2.set_title("B. Enrichment Ratios", fontweight="bold")
    ax2.legend()
    for i, (r, c) in enumerate(zip(ratios, enrich_sorted)):
        ax2.text(r + 0.05, i, f"{r}x ({c[1]['count']} genes)", va="center", fontsize=9)

# Panel C: Top genes
ax3 = axes[1, 0]
ax3.barh([g[0][:20] for g in top_genes[:15]], [len(g[1]) for g in top_genes[:15]], color="steelblue")
ax3.set_title("C. Top 15 Affected Genes", fontweight="bold")
ax3.set_xlabel("Number of Variants")

# Panel D: Summary
ax4 = axes[1, 1]
ax4.axis("off")
summary = f"""SUMMARY STATISTICS:

Total variants: {len(variants):,}
Affected genes: {len(gene_impacts):,}
HIGH impact: {len(high_impact)}
MODERATE: {impact_counter['MODERATE']:,}
LOW: {impact_counter['LOW']:,}
Missense: {len(missense_variants):,}
Enriched categories: {len(enrich_sorted)}

TOP ENRICHED PROCESSES:
"""
for cat, d in enrich_sorted[:3]:
    summary += f"  {cat}: {d['ratio']}x ({d['count']} genes)\n"

summary += "\nTOP CONSEQUENCES:\n"
for cons, count in consequence_counter.most_common(5):
    summary += f"  {cons}: {count:,}\n"

ax4.text(0.05, 0.95, summary, transform=ax4.transAxes, fontsize=10,
         verticalalignment="top", fontfamily="monospace",
         bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

plt.tight_layout()
plt.savefig(f"{OUTDIR}/figures/functional_summary.png", dpi=200, bbox_inches="tight", facecolor="white")
plt.close()
print("  Multi-panel summary saved")

# Impact pie separately
fig2, ax = plt.subplots(figsize=(8, 8))
ax.pie(sizes, labels=labels, autopct="%1.1f%%", colors=["#e74c3c", "#f39c12", "#3498db"],
       explode=(0.05, 0, 0), shadow=True, textprops={"fontsize": 12})
ax.set_title("Variant Impact Distribution\nC. albicans SRR7801919", fontweight="bold", fontsize=14)
plt.savefig(f"{OUTDIR}/figures/impact_pie.png", dpi=150, bbox_inches="tight", facecolor="white")
plt.close()
print("  Impact pie saved")

# HTML Report
print("\nGenerating HTML report...")

def embed(p):
    if os.path.exists(p):
        with open(p, "rb") as f:
            return base64.b64encode(f.read()).decode()
    return ""

sum_b64 = embed(f"{OUTDIR}/figures/functional_summary.png")

# Build tables
enrich_rows = ""
for cat, d in enrich_sorted:
    enrich_rows += f'<tr><td>{cat}</td><td>{d["count"]}</td><td><strong>{d["ratio"]}x</strong></td><td>{", ".join(d["genes"][:3])}</td></tr>'

priority_rows = ""
for i, v in enumerate(sorted(high_impact, key=lambda x: x["score"], reverse=True)[:15], 1):
    priority_rows += f'<tr><td>{i}</td><td><strong>{v["gene"]}</strong></td><td style="color:#e74c3c;font-weight:bold">{v["cons"]}</td><td>{v["score"]}</td></tr>'

today = date.today().strftime("%B %d, %Y")

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="author" content="Kelton Guimaraes">
    <title>Functional Analysis — SRR7801919</title>
    <style>
        * {{ margin:0; padding:0; box-sizing:border-box; }}
        body {{ font-family:'Segoe UI',Arial,sans-serif; color:#333; background:#f5f6fa; }}
        .header {{ background:linear-gradient(135deg,#2c3e50,#8e44ad); color:white; padding:40px; text-align:center; }}
        .header h1 {{ font-size:32px; margin-bottom:10px; }}
        .container {{ max-width:1100px; margin:0 auto; padding:30px; }}
        .stats {{ display:grid; grid-template-columns:repeat(4,1fr); gap:20px; margin:30px 0; }}
        .stat {{ background:white; padding:25px; border-radius:12px; text-align:center; box-shadow:0 2px 8px rgba(0,0,0,0.08); }}
        .stat-val {{ font-size:28px; font-weight:bold; color:#2c3e50; }}
        .stat-lbl {{ font-size:12px; color:#7f8c8d; text-transform:uppercase; letter-spacing:1px; margin-top:8px; }}
        .section {{ background:white; padding:30px; border-radius:12px; margin:25px 0; box-shadow:0 2px 8px rgba(0,0,0,0.08); }}
        .section h2 {{ color:#2c3e50; border-bottom:2px solid #8e44ad; padding-bottom:10px; margin-bottom:15px; }}
        img {{ width:100%; border-radius:8px; margin:15px 0; }}
        table {{ width:100%; border-collapse:collapse; margin:15px 0; }}
        th {{ background:#8e44ad; color:white; padding:12px; text-align:left; }}
        td {{ padding:12px; border-bottom:1px solid #ecf0f1; }}
        tr:hover {{ background:#f8f9fa; }}
        .footer {{ text-align:center; padding:30px; color:#999; font-size:13px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 Functional Variant Impact Analysis</h1>
        <p style="font-size:18px;margin-top:10px;"><strong>Candida albicans</strong> Clinical Isolate SRR7801919</p>
        <p style="margin-top:8px;">Analyst: <strong>Kelton Guimaraes</strong> | {today}</p>
    </div>

    <div class="container">
        <div class="stats">
            <div class="stat"><div class="stat-val">{len(variants):,}</div><div class="stat-lbl">Annotated Variants</div></div>
            <div class="stat"><div class="stat-val">{len(gene_impacts):,}</div><div class="stat-lbl">Affected Genes</div></div>
            <div class="stat"><div class="stat-val">{len(high_impact)}</div><div class="stat-lbl">HIGH Impact</div></div>
            <div class="stat"><div class="stat-val">{len(enrich_sorted)}</div><div class="stat-lbl">Enriched Categories</div></div>
        </div>

        <div class="section">
            <h2>📊 Multi-Panel Summary</h2>
            <img src="data:image/png;base64,{sum_b64}" alt="Summary Figure">
        </div>

        <div class="section">
            <h2>⭐ High Priority Variants</h2>
            <p>Variants scored by impact: STOP_GAINED/FRAMESHIFT = 10, SPLICE_SITE = 8</p>
            <table>
                <tr><th>Rank</th><th>Gene</th><th>Consequence</th><th>Score</th></tr>
                {priority_rows}
            </table>
        </div>

        <div class="section">
            <h2>🔬 Gene Ontology Enrichment</h2>
            <p>Categories overrepresented among mutated genes. Ratio &gt; 1 = more mutations than expected.</p>
            <table>
                <tr><th>Category</th><th>Genes</th><th>Enrichment</th><th>Examples</th></tr>
                {enrich_rows}
            </table>
        </div>

        <div class="section">
            <h2>📈 Consequence Types</h2>
            <table>
                <tr><th>Consequence</th><th>Count</th></tr>
"""
for cons, count in consequence_counter.most_common(10):
    html += f'<tr><td>{cons}</td><td>{count:,}</td></tr>'

html += f"""
            </table>
        </div>
    </div>

    <div class="footer">
        <p><strong>Kelton Guimaraes</strong> | Functional Variant Analysis v1.0</p>
        <p>Methodology: VEP annotation + custom GO enrichment + impact scoring</p>
    </div>
</body>
</html>"""

with open(f"{OUTDIR}/functional_report.html", "w") as f:
    f.write(html)
report_size = os.path.getsize(f"{OUTDIR}/functional_report.html") / 1024

print(f"  Report saved ({report_size:.0f} KB)")

print("\n" + "=" * 60)
print("  ANALYSIS COMPLETE")
print("=" * 60)
print(f"""
Output files:
  {OUTDIR}/enrichment/go_enrichment.txt
  {OUTDIR}/prioritization/high_priority.txt
  {OUTDIR}/figures/functional_summary.png
  {OUTDIR}/figures/impact_pie.png
  {OUTDIR}/functional_report.html
""")
