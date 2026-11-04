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
| Production method | Imaging-based · Sequencing-based |
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
| accession · title · description · organism · sample count · release date · PubMed ID | modality · platform · imaging vs sequencing · disease category · subtype · tissue · cell count |

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

A record whose *title* says it spans many cancer types (pan-cancer, "30 cancer types") is filed as
pan-cancer even when its abstract names one cancer first.

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
