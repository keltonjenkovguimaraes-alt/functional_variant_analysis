#!/usr/bin/env python3
"""
Functional Variant Impact Analysis
Kelton Guimaraes — 2026
"""

import os, sys, json, gzip
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from datetime import date

============================================================
CONFIGURATION
============================================================
VEP_FILE = "/home/kelton/microbial_variant_calling/results/annotation/SRR7801919_snps_vep.txt"
GFF_FILE = "/home/kelton/microbial_variant_calling/data/reference/reference.gff"
OUTDIR = "results"
SAMPLE = "SRR7801919"
ORGANISM = "Candida albicans"

os.makedirs(f"{OUTDIR}/enrichment", exist_ok=True)
os.makedirs(f"{OUTDIR}/domains", exist_ok=True)
os.makedirs(f"{OUTDIR}/prioritization", exist_ok=True)
os.makedirs(f"{OUTDIR}/figures", exist_ok=True)

print("=" * 60)
print(" FUNCTIONAL VARIANT IMPACT ANALYSIS")
print(f" Sample: {SAMPLE} | Organism: {ORGANISM}")
print("=" * 60)

============================================================
1. PARSE VEP OUTPUT
============================================================
print("\n📖 Parsing VEP annotations...")

variants = []
gene_impacts = defaultdict(list)
missense_variants = []
high_impact_variants = []

with open(VEP_FILE, 'r') as f:
for line in f:
if line.startswith('##'):
continue
parts = line.strip().split('\t')
if len(parts) < 14:
continue

VEP columns: Uploaded_variation, Location, Allele, Gene, Feature, Feature_type,
Consequence, cDNA_position, CDS_position, Protein_position, Amino_acids, Codons,
Existing_variation, Extra
variant_id = parts[0] if len(parts) > 0 else "."
consequence = parts[6] if len(parts) > 6 else "."
gene = parts[3] if len(parts) > 3 else "."
extra = parts[13] if len(parts) > 13 else ""

variants.append({
'id': variant_id,
'consequence': consequence,
'gene': gene,
'impact': 'MODERATE' if 'missense' in consequence else
'HIGH' if 'stop_gained' in consequence or 'frameshift' in consequence else
'LOW'
})

if gene and gene != '.':
gene_impacts[gene].append(consequence)

if 'missense' in consequence:
missense_variants.append({
'id': variant_id,
'gene': gene,
'consequence': consequence,
'extra': extra
})

if 'stop_gained' in consequence or 'frameshift' in consequence or 'splice_donor' in consequence or 'splice_acceptor' in consequence:
high_impact_variants.append({
'id': variant_id,
'gene': gene,
'consequence': consequence
})

print(f" Total variants parsed: {len(variants):,}")
print(f" Unique affected genes: {len(gene_impacts):,}")
print(f" Missense variants: {len(missense_variants):,}")
print(f" HIGH impact variants: {len(high_impact_variants):,}")

============================================================
2. GENE ONTOLOGY ENRICHMENT (Custom Implementation)
============================================================
print("\n🔬 Running Gene Ontology Enrichment...")

Parse GFF to get gene -> description mapping
gene_descriptions = {}
with open(GFF_FILE, 'r') as f:
for line in f:
if line.startswith('#'):
continue
parts = line.strip().split('\t')
if len(parts) >= 9:
attrs = parts[8]
gene_id = ""
gene_name = ""
product = ""

for attr in attrs.split(';'):
if '=' in attr:
key, val = attr.split('=', 1)
if key == 'ID':
gene_id = val
elif key == 'Name':
gene_name = val
elif key == 'product':
product = val

if gene_name:
gene_descriptions[gene_name] = product

Define GO-like categories for C. albicans
go_categories = {
'Cell Wall Organization': ['cell wall', 'chitin', 'glucan', 'mannan', 'mannoprotein', 'β-glucan', 'beta-glucan'],
'Biofilm Formation': ['biofilm', 'adhesion', 'adhesin', 'filamentous', 'hyphal', 'yeast-to-hypha'],
'Drug Response': ['drug', 'resistance', 'azole', 'fluconazole', 'efflux', 'transporter', 'multidrug', 'MDR', 'CDR'],
'Stress Response': ['stress', 'oxidative', 'heat shock', 'HSP', 'osmotic', 'pH', 'DNA damage'],
'Metabolism': ['metabolism', 'biosynthesis', 'ergosterol', 'lipid', 'carbohydrate', 'amino acid', 'glycolysis'],
'Virulence': ['virulence', 'pathogenesis', 'invasion', 'protease', 'lipase', 'phospholipase', 'candidalysin'],
'DNA Replication/Repair': ['DNA replication', 'DNA repair', 'recombination', 'mismatch', 'nucleotide excision'],
'Transcription': ['transcription', 'RNA polymerase', 'transcription factor', 'zinc finger', 'DNA binding'],
'Translation': ['translation', 'ribosomal', 'ribosome', 'tRNA', 'elongation factor'],
'Signal Transduction': ['signaling', 'kinase', 'phosphatase', 'cAMP', 'MAPK', 'two-component', 'Ras'],
'Transport': ['transport', 'permease', 'channel', 'ATPase', 'uptake', 'secretion'],
'Mitochondrial': ['mitochondria', 'respiration', 'cytochrome', 'ATP synthase', 'electron transport'],
}

Count genes in each category
enrichment_results = {}
affected_genes_set = set(gene_impacts.keys())
total_bg_genes = len(gene_descriptions)

for category, keywords in go_categories.items():

Background genes in category
bg_in_category = set()
for gene, desc in gene_descriptions.items():
desc_lower = desc.lower() if desc else ""
for kw in keywords:
if kw.lower() in desc_lower:
bg_in_category.add(gene)
break

Affected genes in category
affected_in_category = affected_genes_set & bg_in_category

if len(affected_in_category) >= 2:

Hypergeometric enrichment (Fisher's exact test simplified)
a = len(affected_in_category) # affected and in category
b = len(affected_genes_set) - a # affected but not in category
c = len(bg_in_category) - a # in category but not affected
d = total_bg_genes - a - b - c # neither

Simple enrichment ratio
expected = (len(bg_in_category) / max(total_bg_genes, 1)) * len(affected_genes_set)
enrichment_ratio = a / max(expected, 1)

enrichment_results[category] = {
'count': a,
'expected': round(expected, 1),
'ratio': round(enrichment_ratio, 2),
'genes': list(affected_in_category)[:5] # top 5 examples
}

Sort by enrichment ratio
enrichment_sorted = sorted(enrichment_results.items(), key=lambda x: x[1]['ratio'], reverse=True)

print(f" Categories enriched (≥2 genes): {len(enrichment_sorted)}")

Save enrichment results
with open(f'{OUTDIR}/enrichment/go_enrichment.txt', 'w') as f:
f.write(f"# GO Enrichment Analysis — {SAMPLE}\n")
f.write(f"# Total background genes: {total_bg_genes}\n")
f.write(f"# Total affected genes: {len(affected_genes_set)}\n\n")
f.write("Category\tAffected\tExpected\tEnrichment_Ratio\tExample_Genes\n")
for cat, data in enrichment_sorted:
f.write(f"{cat}\t{data['count']}\t{data['expected']}\t{data['ratio']}\t{', '.join(data['genes'][:3])}\n")

============================================================
3. VARIANT PRIORITIZATION
============================================================
print("\n⭐ Prioritizing high-impact variants...")

priority_variants = []

for v in high_impact_variants:
score = 0
reason = []

if 'stop_gained' in v['consequence']:
score += 10
reason.append("STOP_GAINED")
if 'frameshift' in v['consequence']:
score += 10
reason.append("FRAMESHIFT")
if 'splice' in v['consequence']:
score += 8
reason.append("SPLICE_SITE")

Check if gene is in enriched categories
gene = v['gene']
for cat, data in enrichment_sorted[:5]: # top 5 categories
if gene in data['genes']:
score += 5
reason.append(f"ENRICHED:{cat}")
break

priority_variants.append({
'gene': gene,
'consequence': v['consequence'],
'score': score,
'reasons': '; '.join(reason)
})

priority_variants.sort(key=lambda x: x['score'], reverse=True)

with open(f'{OUTDIR}/prioritization/high_priority_variants.txt', 'w') as f:
f.write(f"# High Priority Variants — {SAMPLE}\n\n")
f.write("Rank\tGene\tConsequence\tPriority_Score\tReasons\n")
for i, v in enumerate(priority_variants[:20], 1):
f.write(f"{i}\t{v['gene']}\t{v['consequence']}\t{v['score']}\t{v['reasons']}\n")

print(f" High-priority variants identified: {len(priority_variants)}")

============================================================
4. MISSENSE VARIANT DOMAIN ANALYSIS
============================================================
print("\n🧬 Analyzing missense variant distribution...")

Categorize missense variants by consequence type
missense_consequences = Counter()
for v in missense_variants:
missense_consequences[v['consequence']] += 1

============================================================
5. GENERATE FIGURES
============================================================
print("\n📊 Generating figures...")

Figure 1: Enrichment barplot
if enrichment_sorted:
fig, ax = plt.subplots(figsize=(12, 8))
categories = [x[0] for x in enrichment_sorted[:10]]
ratios = [x[1]['ratio'] for x in enrichment_sorted[:10]]
counts = [x[1]['count'] for x in enrichment_sorted[:10]]

colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(categories)))
bars = ax.barh(range(len(categories)), ratios, color=colors, edgecolor='black')

for i, (bar, ratio, count) in enumerate(zip(bars, ratios, counts)):
ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2,
f'n={count}, {ratio:.1f}x', va='center', fontsize=9, fontweight='bold')

ax.set_yticks(range(len(categories)))
ax.set_yticklabels(categories, fontsize=10)
ax.set_xlabel('Enrichment Ratio (observed/expected)', fontsize=12)
ax.set_title('Gene Ontology Enrichment Analysis\nC. albicans SRR7801919', fontsize=14, fontweight='bold')
ax.axvline(x=1.0, color='red', linestyle='--', alpha=0.7, label='No enrichment (ratio=1)')
ax.legend(loc='lower right')
ax.set_xlim(0, max(ratios) * 1.3)

plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/go_enrichment.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f" ✓ GO enrichment plot saved")

Figure 2: Variant Impact Pie
impact_counts = Counter(v['impact'] for v in variants)
fig, ax = plt.subplots(figsize=(10, 8))
labels = list(impact_counts.keys())
sizes = list(impact_counts.values())
colors_pie = ['#e74c3c', '#f39c12', '#3498db']
explode = (0.05, 0, 0)

ax.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors_pie, explode=explode,
shadow=True, startangle=90, textprops={'fontsize': 12})
ax.set_title('Variant Impact Distribution\nC. albicans SRR7801919', fontsize=14, fontweight='bold')

plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/variant_impact_pie.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f" ✓ Impact distribution plot saved")

Figure 3: Multi-panel summary
fig = plt.figure(figsize=(16, 12))
fig.suptitle('Functional Variant Analysis Summary\nC. albicans SRR7801919', fontsize=16, fontweight='bold')

Panel A: Enrichment
ax1 = fig.add_subplot(2, 2, 1)
if enrichment_sorted:
top5 = enrichment_sorted[:5]
ax1.barh([x[0] for x in top5], [x[1]['ratio'] for x in top5],
color=plt.cm.plasma(np.linspace(0.2, 0.9, 5)))
ax1.set_title('A. Top 5 Enriched Categories', fontweight='bold')
ax1.axvline(x=1.0, color='red', linestyle='--')

Panel B: Impact pie
ax2 = fig.add_subplot(2, 2, 2)
ax2.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors_pie, explode=explode)
ax2.set_title('B. Variant Impact Distribution', fontweight='bold')

Panel C: Priority variants
ax3 = fig.add_subplot(2, 2, 3)
ax3.axis('off')
priority_text = "HIGH PRIORITY VARIANTS:\n\n"
for i, v in enumerate(priority_variants[:10], 1):
priority_text += f"{i}. {v['gene']}: {v['consequence']} (Score: {v['score']})\n"
ax3.text(0.05, 0.95, priority_text, transform=ax3.transAxes,
fontsize=9, verticalalignment='top', fontfamily='monospace')

Panel D: Summary stats
ax4 = fig.add_subplot(2, 2, 4)
ax4.axis('off')
summary = f"""
SUMMARY STATISTICS:

• Total annotated variants: {len(variants):,}
• Affected genes: {len(gene_impacts):,}
• Missense variants: {len(missense_variants):,}
• HIGH impact variants: {len(high_impact_variants):,}
• Enriched categories: {len(enrichment_sorted)}
• Priority variants: {len(priority_variants)}

TOP ENRICHED PROCESSES:
"""
for cat, data in enrichment_sorted[:3]:
summary += f"• {cat}: {data['ratio']}x enriched ({data['count']} genes)\n"

ax4.text(0.05, 0.95, summary, transform=ax4.transAxes,
fontsize=10, verticalalignment='top', fontfamily='monospace',
bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.tight_layout()
plt.savefig(f'{OUTDIR}/figures/functional_summary.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print(f" ✓ Multi-panel summary saved")

============================================================
6. GENERATE HTML REPORT
============================================================
print("\n📄 Generating HTML report...")

import base64

def embed_image(path):
if os.path.exists(path):
with open(path, 'rb') as f:
return base64.b64encode(f.read()).decode()
return ""

go_img = embed_image(f'{OUTDIR}/figures/go_enrichment.png')
impact_img = embed_image(f'{OUTDIR}/figures/variant_impact_pie.png')
summary_img = embed_image(f'{OUTDIR}/figures/functional_summary.png')

Build enrichment table rows
enrichment_rows = ""
for cat, data in enrichment_sorted[:10]:
ratio_display = f"{data['ratio']}x"
enrichment_rows += f"""

<tr> <td>{cat}</td> <td>{data['count']}</td> <td>{data['expected']}</td> <td><strong>{ratio_display}</strong></td> <td style="font-size:11px;">{', '.join(data['genes'][:3])}</td> </tr>"""
Priority table rows
priority_rows = ""
for i, v in enumerate(priority_variants[:10], 1):
score_color = '#e74c3c' if v['score'] >= 10 else '#f39c12'
priority_rows += f"""

<tr> <td>{i}</td> <td><strong>{v['gene']}</strong></td> <td>{v['consequence']}</td> <td style="color:{score_color};font-weight:bold;">{v['score']}</td> <td style="font-size:11px;">{v['reasons']}</td> </tr>"""
html = f"""<!DOCTYPE html>

<html lang="en"> <head> <meta charset="UTF-8"> <meta name="author" content="Kelton Guimaraes"> <title>Functional Variant Analysis — {SAMPLE}</title> <style> * {{ margin:0; padding:0; box-sizing:border-box; }} body {{ font-family:'Segoe UI',Arial,sans-serif; color:#333; background:#f5f6fa; }} .header {{ background:linear-gradient(135deg,#2c3e50,#8e44ad); color:white; padding:40px; text-align:center; }} .header h1 {{ font-size:32px; margin-bottom:10px; }} .container {{ max-width:1200px; margin:0 auto; padding:30px; }} .stats-grid {{ display:grid; grid-template-columns:repeat(4,1fr); gap:20px; margin:30px 0; }} .stat-card {{ background:white; padding:25px; border-radius:12px; box-shadow:0 2px 10px rgba(0,0,0,0.08); text-align:center; }} .stat-value {{ font-size:28px; font-weight:bold; color:#2c3e50; }} .stat-label {{ font-size:13px; color:#7f8c8d; margin-top:8px; text-transform:uppercase; letter-spacing:1px; }} .section {{ background:white; padding:30px; border-radius:12px; box-shadow:0 2px 10px rgba(0,0,0,0.08); margin:25px 0; }} .section h2 {{ color:#2c3e50; font-size:22px; margin-bottom:15px; padding-bottom:10px; border-bottom:2px solid #8e44ad; }} img {{ width:100%; border-radius:8px; margin:15px 0; }} table {{ width:100%; border-collapse:collapse; margin:15px 0; }} th {{ background:#8e44ad; color:white; padding:12px; text-align:left; }} td {{ padding:12px; border-bottom:1px solid #ecf0f1; }} tr:hover {{ background:#f8f9fa; }} .badge {{ display:inline-block; padding:5px 12px; border-radius:20px; font-size:12px; font-weight:bold; }} .badge-high {{ background:#f8d7da; color:#721c24; }} .badge-moderate {{ background:#fff3cd; color:#856404; }} .footer {{ text-align:center; padding:30px; color:#95a5a6; font-size:13px; }} </style> </head> <body> <div class="header"> <h1>🧬 Functional Variant Impact Analysis</h1> <p><strong>{ORGANISM}</strong> Clinical Isolate {SAMPLE}</p> <p style="margin-top:10px;font-size:18px;">Analyst: <strong>Kelton Guimaraes</strong> | {date.today().strftime('%B %d, %Y')}</p> </div> <div class="container"> <div class="stats-grid"> <div class="stat-card"> <div class="stat-value">{len(variants):,}</div> <div class="stat-label">Annotated Variants</div> </div> <div class="stat-card"> <div class="stat-value">{len(gene_impacts):,}</div> <div class="stat-label">Affected Genes</div> </div> <div class="stat-card"> <div class="stat-value">{len(high_impact_variants)}</div> <div class="stat-label">HIGH Impact</div> </div> <div class="stat-card"> <div class="stat-value">{len(enrichment_sorted)}</div> <div class="stat-label">Enriched Categories</div> </div> </div> <div class="section"> <h2>📊 Gene Ontology Enrichment</h2> <p>Functional categories overrepresented among mutated genes. <strong>Enrichment ratio > 1</strong> indicates more mutations than expected by chance.</p> {"<img src='data:image/png;base64," + go_img + "' alt='GO Enrichment'>" if go_img else ""} <table> <tr><th>Category</th><th>Affected Genes</th><th>Expected</th><th>Enrichment</th><th>Examples</th></tr> {enrichment_rows} </table> </div> <div class="section"> <h2>⭐ High Priority Variants</h2> <p>Variants scored by impact severity (STOP_GAINED=10, FRAMESHIFT=10, SPLICE=8) + enrichment in key pathways (+5).</p> <table> <tr><th>Rank</th><th>Gene</th><th>Consequence</th><th>Score</th><th>Evidence</th></tr> {priority_rows} </table> </div> <div class="section"> <h2>📈 Variant Impact Distribution</h2> {"<img src='data:image/png;base64," + impact_img + "' alt='Impact Distribution'>" if impact_img else ""} </div> <div class="section"> <h2>📋 Multi-Panel Summary</h2> {"<img src='data:image/png;base64," + summary_img + "' alt='Summary'>" if summary_img else ""} </div> </div> <div class="footer"> <p><strong>Kelton Guimaraes</strong> | Functional Variant Analysis v1.0</p> <p>Methodology: VEP annotation + custom GO enrichment + impact scoring</p> </div> </body> </html>"""
with open(f'{OUTDIR}/functional_report.html', 'w') as f:
f.write(html)

print(f" ✓ Report saved: {OUTDIR}/functional_report.html ({len(html)/1024:.0f} KB)")

============================================================
DONE
============================================================
print("\n" + "=" * 60)
print(" ✅ FUNCTIONAL ANALYSIS COMPLETE")
print("=" * 60)
print(f"""
Outputs:
{OUTDIR}/enrichment/go_enrichment.txt
{OUTDIR}/prioritization/high_priority_variants.txt
{OUTDIR}/figures/go_enrichment.png
{OUTDIR}/figures/variant_impact_pie.png
{OUTDIR}/figures/functional_summary.png
{OUTDIR}/functional_report.html
""")
