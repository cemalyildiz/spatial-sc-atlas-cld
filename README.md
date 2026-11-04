# Spatial &amp; Single-Cell Atlas

A browsable, filterable catalog of **spatial transcriptomics** and **single-cell RNA-seq**
datasets for cancer and other diseases, harvested straight from
[NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) and
[EBI ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress).

**→ [cemalyildiz.github.io/spatial-sc-atlas-cld](https://cemalyildiz.github.io/spatial-sc-atlas-cld/)**

Built for the question "which public dataset should I actually use?" rather than for reading
abstracts. The filter that usually matters first is **imaging-based vs sequencing-based**, and the
records worth the most are the **paired** ones — a single deposit carrying both spatial and
single-cell data, which is what you need for deconvolution, label transfer and spot annotation.

## Filters

| Filter | Values |
| --- | --- |
| Modality | Paired (spatial + single-cell) · Spatial only · Single-cell only |
| Production method | Imaging-based · Sequencing-based · Not stated (spatial, but the record never names its platform) |
| Platform | Xenium, MERFISH/MERSCOPE, CosMx, seqFISH, STARmap, in situ sequencing, Visium, Visium HD, GeoMx DSP, Slide-seq, Stereo-seq, Slide-tags, DBiT-seq, 10x Chromium, … |
| Disease group | Cancer · Other disease · Unspecified |
| Disease / cancer type | NSCLC, SCLC, mesothelioma, breast, colorectal, glioma/GBM, meningioma, pancreatic, melanoma, …, plus non-cancer categories (fibrosis, neurodegenerative, autoimmune, cardiovascular, …) |
| Subtype (multi-value) | LUAD, LUSC, LCNEC, TNBC, HER2-positive, luminal, ccRCC, CRPC, NEPC, ESCC, PDAC, HCC, DLBCL, B-ALL, T-ALL, CLL, MIBC, HGSOC, HNSCC, DIPG, plus molecular labels — EGFR-mutant, ALK-rearranged, KRAS-mutant, IDH-mutant, MSI-high, HPV-positive |
| Study type | Biological study · Method / tool development (the data are usually a demonstration set; titled `TOOL: …`, `Benchmarking …`, or the abstract says "we develop a method/tool/…") |
| Tissue / organ | Lung, brain/CNS, breast, colon/intestine, liver, kidney, … |
| Organism | Human · Mouse · Human &amp; Mouse |
| Source database | GEO · ArrayExpress |
| Access | Open · Controlled |
| Release year | Range |
| Minimum cells | ≥ 10k / 50k / 100k / 500k |

Filter state lives in the URL, so a filtered view is a shareable link. Results export as CSV or
JSON, and the accessions of the whole filtered set copy to the clipboard in one click — paste them
straight into a download script.

## Which fields are trustworthy

| Verbatim from the source database | Derived here, automatically |
| --- | --- |
| accession · title · description · organism · sample count · release date · PubMed ID · ArrayExpress study type | modality · platform (and `platform_source`) · imaging vs sequencing · disease category · subtype · study type · tissue · cell count |

Derived fields are pattern-matched from each record's title and description by the vocabularies at
the top of [`scripts/harvest.py`](scripts/harvest.py). They are good enough to filter and browse
with, and wrong often enough that you should open the source record — one click from every card —
before a dataset enters an analysis.

A `Paired` badge means one of two things: the record's own text mentions both a spatial platform and
single-cell data, or the record belongs to a study GEO split across several series and a sibling
carries the other half. The second kind shows a `+N series` chip and lists the sibling accessions.
GEO's API does not expose SuperSeries links (`relations` is empty for subseries), so siblings are
found by grouping records that share the study abstract verbatim. ArrayExpress companion studies
(e.g. E-MTAB-13526 scRNA-seq + E-MTAB-13530 Visium) do not share an abstract, so they are linked when
their titles match once assay words are removed and their accession numbers are close. Both are
heuristics, so check before relying on them.

GEO subseries whose titles differ only in the trailing `[assay]` tag are linked as well, for studies
whose subseries write their own abstract (GSE292469 snRNA-seq + GSE342824 Xenium). The tag also
decides what each subseries holds: `[scRNA-seq]` → single-cell (paired only through a sibling),
`[ST]`/`[Visium]`/`[Xenium]` → spatial, and a tag that is neither (bulk RNA-seq, ATAC including
single-cell ATAC, TCR/BCR/VDJ repertoire, ChIP, CUT&Tag, methylation, RRBS/WGBS, aCGH) → left out, as
is an untagged record whose title is about ATAC only. A plain `[RNA-seq]` inside a multi-series study is
the bulk companion and is dropped — unless every other part of the study is a non-expression assay
(GSE118389 `[RNA-Seq]` beside GSE118303 `[WES]`), where it is the only expression data and stays.
Groups are formed before anything is dropped, so removing a bulk subseries never splits a study.

The platform of a spatial record comes from the strongest evidence available: its own tag or title,
then the **file names it deposited** (`tissue_positions.csv`/`scalefactors_json` → Visium,
`transcripts.parquet` → Xenium, `.dcc`/`.pkc` → GeoMx, `exprMat_file` → CosMx, `cell_by_gene` →
MERSCOPE, `.gef` → Stereo-seq, `bead_locations` → Slide-seq, `square_008um` → Visium HD), then its
processing/protocol text (Space Ranger, AtoMx, Vizgen), and only then the abstract. The record's
`platform_source` field (shown as *Platform evidence*) says which one decided. File names are read
from the GEO series and sample records and the BioStudies file list, fetched once per accession and
cached in `data/enrichment.json`; bump `ENRICH_VERSION` in `harvest.py` after changing the signatures
to fetch everything again. A record whose evidence names nothing, or conflicting platforms, is
`Spatial (unspecified)` with production method *Not stated* rather than counted as sequencing-based.

The cancer named in the *title* wins over the first one found in the abstract, except that a title
saying only "lung cancer" does not override NSCLC/SCLC from the abstract. A title that says it spans
many cancer types (pan-cancer, "30 cancer types") is filed as pan-cancer. When the text names no
organ, the tissue is taken from the disease (glioma → brain, leukaemia → bone marrow / blood, …).
Organ-specific subtypes are kept only on a matching disease: a breast-cancer study that cites lung
adenocarcinoma does not carry `LUAD`. Molecular labels (EGFR-, KRAS-mutant, HER2-positive, IDH, MSI,
HPV) occur in several cancers and are kept whatever the disease.

Release dates are ISO (`YYYY-MM-DD`) for both sources. For ArrayExpress, organism, sample count and
study type come from the BioStudies study record; the study type also settles modality (a study typed
"RNA-seq of coding RNA from single cells" is single-cell even when its title names its Visium partner),
and a study typed only as ATAC/ChIP/methylation is left out. Cell counts exist only where a record
states one, so the *Minimum cells* filter hides every record without a stated count.

Improving a label means editing a regex in `harvest.py` and re-running the harvest — not hand-editing
`data/datasets.json`, which is regenerated on every run.

One rule when adding subtype patterns: short abbreviations collide. `ILC` is innate lymphoid cell
far more often than invasive lobular carcinoma, `IDC` is also interdigitating dendritic cell, `FL`
and `MCL` hit gene and cell-type names. Those are matched by full phrase only. Check a candidate
abbreviation against the existing catalog before adding it bare.

## How it stays current

`.github/workflows/refresh.yml` re-runs the harvest every Monday at 04:17 UTC and commits the
refreshed catalog **straight to `main`**, so the site is always current without anyone approving a
pull request. Nothing is hand-edited in between: `scripts/harvest.py` regenerates `data/` from
scratch on every run, and refuses to overwrite the catalog if the harvest returns fewer than 200
records, so a failed or throttled source leaves the previous catalog in place.

```
Actions → Refresh catalog → Run workflow
  mode: direct  commit straight to main          (default, used by the schedule)
  mode: pr      open a pull request for review   (use when changing the classification rules)
```

Use `mode: pr` after editing the vocabularies in `scripts/harvest.py`, when the point is to read the
label diff before it goes live.

Optionally add a repository secret `NCBI_API_KEY`
([get one here](https://www.ncbi.nlm.nih.gov/account/settings/)) to raise the NCBI rate limit from
3 to 10 requests per second. The harvest works without it, just more slowly.

## Running the harvest locally

```bash
python3 scripts/harvest.py     # stdlib only, no dependencies
python3 scripts/harvest.py --reclassify   # re-apply label rules to the existing catalog, no network
python3 scripts/harvest.py --reclassify --fetch   # same, but fetch missing platform evidence first
python3 -m http.server 8000    # then open http://localhost:8000
```

The harvest refuses to overwrite `data/datasets.json` if it comes back with fewer than 200 records,
so a half-finished network run cannot wipe the catalog.

## Layout

```
index.html            single-page app, no build step
assets/style.css      design tokens, light + dark
assets/app.js         filtering, comparison, charts — no dependencies
scripts/harvest.py    GEO + ArrayExpress harvest and classification
data/datasets.json    the catalog (generated; do not hand-edit)
data/meta.json        counts and the generation timestamp
```

## Citing

This site is an index, not a source. Cite the original study and its accession. Check each record's
terms of use before redistributing data — a few carry controlled access even where the metadata is
public.

## License

Code: MIT. Dataset metadata belongs to the depositing authors and their source repositories.
