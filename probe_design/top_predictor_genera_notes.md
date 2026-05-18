# Genus composition of top RF predictor families
## AB_lower_contact vs C_elevated binary classification

Generated from QIIME2/SILVA taxonomy (taxonomy_export/taxonomy.tsv) and
genus-level count table (genus_export/genus_table.tsv).
Full data: top_predictor_genera.csv

---

## Summary table

| RF rank | Family | n genera | Dominant genus | % reads | Signal type |
|---|---|---|---|---|---|
| 1 | Deinococcaceae | **1** | *Deinococcus* | 100% | single genus |
| 2 | Intrasporangiaceae | 6 | unclassified + *Ornithinimicrobium* | 51 + 47% | 2 genera |
| 3 | Hungateiclostridiaceae | 14 | *Fastidiosipila* | 57% | distributed |
| 4 | Christensenellaceae | **3** | *Christensenellaceae* R-7 group | 98% | single genus |
| 6 | Eggerthellaceae | 5 | uncultured | 82% | 1–2 genera |
| 7 | UCG-010 | **1** | UCG-010 | 100% | single genus |
| 8 | Micrococcaceae | 13 | *Micrococcus* | 48% | distributed |
| 9 | Peptostreptococcales-Tissierellales | 18 | *Tissierella* | 76% | 1–2 genera |
| 10 | Dietziaceae | **1** | *Dietzia* | 100% | single genus |

---

## Per-family notes

### Deinococcaceae — RF rank 1, C_elevated enriched, env-shared
- **1 genus: *Deinococcus* (100% of reads and ASVs)**
- The family signal is entirely driven by *Deinococcus*.
- *Deinococcus* is famous for extreme resistance to UV radiation, desiccation,
  and ionising radiation (multiple DNA repair mechanisms, carotenoid pigments).
- Enrichment on elevated, sun-exposed sites (Dorsum, Pectoral, Forehead, Neck)
  is ecologically consistent: these sites receive maximal solar UV exposure,
  especially in working horses outdoors.
- Confirmed by Strompfová & Štempelová (2024): highest on equine back (26.5%).
- **Clean probe target**: DEIN-V3 probe targets *Deinococcus* specifically.

### Intrasporangiaceae — RF rank 2, AB_lower_contact enriched, animal-only
- **6 genera; 2 dominate: unclassified (51% reads) + *Ornithinimicrobium* (47%)**
- 98% of reads covered by these two; remaining 4 genera < 2% combined.
- The "unclassified" fraction at genus level likely represents *Intrasporangium*
  or related genera not yet resolved in SILVA 138 — but all belong to
  Intrasporangiaceae s.s. and are phylogenetically tight.
- *Ornithinimicrobium* is an actinomycete originally isolated from human skin;
  also found in soil. Presence on lower-contact sites (pastern, muzzle)
  suggests soil/grass transfer.
- Animal-only: absent from environment swab — may represent host-associated
  actinomycetes selectively colonising skin.
- **Probe note**: INTR-V3 was designed against the full family pool and should
  cover both the unclassified and *Ornithinimicrobium* fractions; verify by
  BLAST against *Ornithinimicrobium* type strain 16S (NCBI).

### Hungateiclostridiaceae — RF rank 3, AB_lower_contact enriched, animal-only
- **14 genera; signal is distributed but *Fastidiosipila* dominates reads**
- *Fastidiosipila* (57% reads, only 17% of ASVs): few but highly abundant strains.
- *Saccharofermentans* (22% reads, 26% ASVs): most diverse genus by ASV count
  but lower per-ASV abundance than *Fastidiosipila*.
- *Ruminiclostridium* (8% reads, 16% ASVs): third contributor.
- All three are strictly anaerobic, fermentative gut bacteria (class Clostridia).
  Their presence on lower-contact skin (muzzle near soil, pastern) is consistent
  with transient faecal/soil contamination or persistent microsite colonisation
  in moist areas.
- Animal-only: not detected in environment swab.
- **Probe note**: HUNG-V3 was designed to exclude Christensenellaceae; validate
  specifically against *Fastidiosipila* spp. 16S sequences.

### Christensenellaceae — RF rank 4, AB_lower_contact enriched, env-shared
- **3 genera; essentially 1 genus: *Christensenellaceae* R-7 group (98% reads)**
- The R-7 group is the dominant lineage in Christensenellaceae globally — it is
  the most commonly detected clade in human gut microbiome studies and is
  strongly heritable (associated with host BMI in GWAS).
- Env-shared: the family was also detected in the environment swab. This may
  reflect its ubiquity in cattle/horse faecal environments (grass/soil
  contamination from livestock).
- Enrichment at lower-contact sites (pastern, muzzle) is consistent with
  direct ground/faecal contact.
- **Probe note**: CHRI-V3 targets *Christensenellaceae* R-7 group specifically;
  specificity vs. Hungateiclostridiaceae confirmed in-study (0% cross-match).

### Eggerthellaceae — RF rank 6, AB_lower_contact enriched, animal-only
- **5 genera; "uncultured" fraction dominates (82% reads, 61% ASVs)**
- *Enterorhabdus* (12% reads): a cultured genus in this family, found in mammalian
  intestines and occasionally on skin.
- High uncultured fraction suggests novel Eggerthellaceae lineages not yet in
  culture collections — common for strict anaerobes.
- Animal-only: absent from environment swab.
- Likely gut-origin commensal colonising skin at lower-contact sites.

### UCG-010 — RF rank 7, AB_lower_contact enriched, animal-only
- **1 genus: UCG-010 (100% of reads and ASVs)**
- UCG-010 (Uncultured Clostridiales Genus) is an entirely uncultured lineage
  in the order Oscillospirales; known from ruminant gut microbiome studies.
- Single genus, single family — the entire signal is one phylogenetic lineage.
- 188 ASVs (high diversity within this single genus) suggest well-established
  colonisation rather than transient contamination.
- Animal-only and completely gut-associated in published references.

### Micrococcaceae — RF rank 8, C_elevated enriched, env-shared
- **13 genera; 3 drive 86% of reads**
- *Micrococcus* (48% reads, 8% ASVs): high per-ASV abundance, few but dominant strains.
- *Glutamicibacter* (29% reads, 4% ASVs): formerly *Arthrobacter* sp., soil actinomycete.
- *Arthrobacter* (9% reads, 6% ASVs): ubiquitous soil/plant-associated actinomycete.
- Env-shared: all three dominant genera are well-known environmental actinomycetes.
  Their enrichment on elevated/harness-contact sites (Dorsum, Neck, Pectoral)
  may reflect deposition from the environment (dust, hay, harness contact).
- Note: *Micrococcus* is a classic skin-associated genus in mammals (including horses);
  its enrichment on elevated sites may reflect different skin physiology
  (drier, less occluded microenvironment) compared to lower-contact sites.

### Peptostreptococcales-Tissierellales — RF rank 9, AB_lower_contact enriched, animal-only
- **18 genera; *Tissierella* dominates (76% reads, 26% ASVs)**
- *Tissierella* is a strictly anaerobic gut bacterium, found in human/animal intestines
  and occasionally in wound infections.
- W5053 (12% reads): an uncultured genus in this order.
- *Peptoniphilus* (3%) and *Helcococcus* (3%): both anaerobes found in skin
  wounds and animal gut.
- Distributed signal (18 genera) but majority driven by one genus.
- Animal-only, anaerobic, gut-associated profile consistent with faecal/soil
  transfer to lower-contact skin sites.

### Dietziaceae — RF rank 10, AB_lower_contact enriched, animal-only
- **1 genus: *Dietzia* (100% of reads and ASVs)**
- *Dietzia* is a lipophilic actinomycete closely related to *Rhodococcus*,
  found in soil, marine sediments, and occasionally animal skin.
- 24 ASVs suggest moderate within-genus diversity.
- Animal-only: not detected in the environment swab despite being soil-associated
  in the literature. This may reflect a horse-specific skin niche or very low
  environmental abundance below the detection threshold.
- Lower-contact enrichment differs from typical *Dietzia* ecology; may indicate
  skin sebum utilisation (Dietziaceae are oleophilic).

---

## Key messages for the paper

1. **Four families are effectively single-genus**: Deinococcaceae (*Deinococcus*),
   Christensenellaceae (R-7 group), UCG-010, Dietziaceae (*Dietzia*).
   The probe sequences for these families target one known genus — no ambiguity.

2. **Two families have two co-dominant genera** (Intrasporangiaceae, Eggerthellaceae).
   Probes must cover both; verify INTR-V3 specificity across *Ornithinimicrobium*
   reference sequences.

3. **Three families have distributed multi-genus signals** (Hungateiclostridiaceae,
   Micrococcaceae, Peptostreptococcales-Tissierellales). For these, a single probe
   captures the dominant genus but misses 13–44% of family reads. Probe pools
   (2–3 oligos) would improve sensitivity.

4. **ASV count ≠ read abundance** (clearest in Hungateiclostridiaceae):
   *Saccharofermentans* has the most ASVs (26%) but *Fastidiosipila* contributes
   the most reads (57%). Report both metrics when characterising family composition.

5. **Uncultured and unclassified fractions are substantial** in Eggerthellaceae (82%)
   and UCG-010 (100%). These represent genuine microbiome diversity not yet
   accessible to culture-based methods — a point for the Discussion.
