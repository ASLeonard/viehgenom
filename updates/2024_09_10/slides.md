---
title: BPC update
author: Alex
institute: ETH Zürich
date: 2024/09/10
output:
  beamer_presentation:
    theme: Boadilla
    keep_tex: true
---

# Proposed agenda

 - Contributed assemblies and QC update (Lloyd & Alex)
 - Pangenome analysis topics
   - bovine diversity/distribution/introgressions etc
   - optimal size of pangenome to minimise reference bias
   - other suggestions?
 - Publication strategy (monolith versus companion papers)
 - Any other business

---

# Contributed assemblies

Firstly, a big thanks to **everyone**!!

. . .

> Over 200 assemblies contributed (with some duplicates)

. . .

Fill in missing metadata at some point (not sure what at the moment)

 - *Bos taurus* subspecies
 - "heritage"/normal/"composite" breed

---

# Contributed assemblies

![reference bias in alignment](2024_09_10/breed_breakdown.svg){ width=100% }

---

# Assembly QC

![reference bias in alignment](2024_09_10/N50.svg){ width=100% }

---

# Assembly QC

Polished CLR/ONT assemblies look okay.

Non-cattle *bovina* have similar divergence w.r.t. cattle.

![reference bias in alignment](2024_09_10/SNPs.svg){ width=100% }

---

# Assembly QC

"Assembly sex" depends on F1 sex && phasing method

![reference bias in alignment](2024_09_10/sex_chromosomes.svg){ width=100% }

---

# QC thresholds

Generally be forgiving.

. . .

| N50 threshold | Assemblies | Loss                                        |
|---------------|------------|---------------------------------------------|
| 0 Mb          | 68         | None                                        |
| 1 Mb          | 64         | Limousine, Aubrac, Abondance, Montbeliarde  |
| 10 Mb         | 62         | Vechur, Sunandini                           |

. . .

No direct estimate of assembly correctness, but most follow expected SNP pattern.

---

# Reality

Substantially heterozygous set of data/methods/etc.

Likely uneven quality, phasing, centromeres/telomeres, etc.

---

# Proposed pangenome analyses

We have a large set of assemblies (**and breeds**), what *can* we say?

. . .

BUSCO genes that are

 - always missing/duplicated
 - appear in different locations

---

# Proposed pangenome analyses

How much bias is reduced with a pangenome using:

 - one of each taurine/indicine/sanga?
 - one of each representative breed?
 - every assembly?

. . .

Is any bias introduced by using too complex of a graph (e.g., non-cattle assemblies)?

---

# Proposed pangenome analyses

Pangenome building:

:::incremental
 - per chromosome (valid for *most* bovina species)
 - re-scaffold reference-guided assemblies to new T2T assembly
 - `minigraph` and `pggb` pangenomes
   - `cactus` interesting from Benedict Paten
 - `pangene` type graph?
:::

---

# Proposed pangenome analyses

Genotype variants (small and structural) in populations?

. . .

Other suggestions?

---

# Publication strategy

Monolith versus companion papers

. . .

Proposed timeline

. . .


Non-compete with specific assembly analyses.

. . .

Authorship --  gathering names/affiliations for sample providers, sequencing, analysis, etc

---

# Data sharing

Sharing all raw/processed data or assemblies only?

. . .

Responsibility for ensuring data is public?

. . .

github/central location for links?

> https://github.com/ASLeonard/viehgenom

> https://cmacphillamy.github.io/BPC_demo/the-bovine-pan-genome-consortium.html

---

# Any other business

More *formal* meeting at PAG?

. . .

```
_____________________
/ Any other business \
\  to bring up now?  /
---------------------
       \   ^__^
        \  (oo)\_______
           (__)\       )\/\
               ||----w |
               ||     ||
```
