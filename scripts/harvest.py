#!/usr/bin/env python3
"""
Harvest spatial transcriptomics and single-cell RNA-seq dataset records from
NCBI GEO (E-utilities) and EBI ArrayExpress (BioStudies API), classify them,
and write data/datasets.json + data/meta.json.

Only stdlib is used so the GitHub Action needs no dependency install.

Provenance note
---------------
Accession, title, description, organism, sample count, release date and PubMed
ID come verbatim from the source databases. Everything else (modality,
platform, platform class, disease category, tissue, cell count) is DERIVED by
the rules in this file from the title + description text, and is therefore a
best-effort label, not a curator's judgement. The site labels these fields as
auto-classified.
"""

import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BIOSTUDIES = "https://www.ebi.ac.uk/biostudies/api/v1/arrayexpress/search"
NCBI_KEY = os.environ.get("NCBI_API_KEY", "").strip()
UA = "spatial-sc-atlas-cld/1.0 (https://github.com/cemalyildiz/spatial-sc-atlas-cld)"

# Polite pacing: 10 req/s with an API key, 3 req/s without.
PAUSE = 0.12 if NCBI_KEY else 0.40


# --------------------------------------------------------------------------
# HTTP helpers
# --------------------------------------------------------------------------
def _request(url, data=None, tries=4):
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(
                url,
                data=data.encode() if isinstance(data, str) else data,
                headers={
                    "User-Agent": UA,
                    "Accept": "application/json",
                    **(
                        {"Content-Type": "application/x-www-form-urlencoded"}
                        if data
                        else {}
                    ),
                },
            )
            with urllib.request.urlopen(req, timeout=90) as resp:
                return json.loads(resp.read().decode("utf-8", "replace"))
        except Exception as exc:  # noqa: BLE001 - network flakiness is expected
            last = exc
            time.sleep(2 * (attempt + 1))
    print(f"  ! giving up on {url[:90]}: {last}", file=sys.stderr)
    return None


def _request_text(url, tries=3):
    """GET a text page; "" for a 404 (nothing deposited), None when unreachable."""
    last = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=90) as resp:
                return resp.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return ""
            last = exc
        except Exception as exc:  # noqa: BLE001
            last = exc
        time.sleep(3 * (attempt + 1))
    print(f"  ! giving up on {url[:90]}: {last}", file=sys.stderr)
    return None


def esearch(term, retmax=6000):
    params = {
        "db": "gds",
        "retmode": "json",
        "retmax": str(retmax),
        "term": term,
    }
    if NCBI_KEY:
        params["api_key"] = NCBI_KEY
    j = _request(EUTILS + "esearch.fcgi?" + urllib.parse.urlencode(params))
    time.sleep(PAUSE)
    if not j:
        return []
    return j.get("esearchresult", {}).get("idlist", []) or []


def esummary(uids, batch=250):
    out = []
    for i in range(0, len(uids), batch):
        chunk = uids[i : i + batch]
        body = {"db": "gds", "retmode": "json", "id": ",".join(chunk)}
        if NCBI_KEY:
            body["api_key"] = NCBI_KEY
        j = _request(EUTILS + "esummary.fcgi", data=urllib.parse.urlencode(body))
        time.sleep(PAUSE)
        if not j or "result" not in j:
            continue
        res = j["result"]
        for uid in res.get("uids", []):
            d = res.get(uid)
            if not d:
                continue
            out.append(
                {
                    "acc": d.get("accession", ""),
                    "title": d.get("title", "") or "",
                    "summary": d.get("summary", "") or "",
                    "taxon": d.get("taxon", "") or "",
                    "gdstype": d.get("gdstype", "") or "",
                    "n": d.get("n_samples", 0) or 0,
                    "pdat": d.get("pdat", "") or "",
                    "pmid": (d.get("pubmedids") or [""])[0],
                    "ftp": d.get("ftplink", "") or "",
                }
            )
        print(f"  esummary {min(i + batch, len(uids))}/{len(uids)}", flush=True)
    return out


def ae_search(query, max_pages=30):
    out = {}
    for page in range(1, max_pages + 1):
        url = (
            BIOSTUDIES
            + "?"
            + urllib.parse.urlencode(
                {"query": query, "pageSize": "100", "page": str(page)}
            )
        )
        j = _request(url)
        time.sleep(0.25)
        if not j:
            break
        hits = j.get("hits") or []
        if not hits:
            break
        for h in hits:
            acc = h.get("accession", "")
            if acc:
                out[acc] = {
                    "acc": acc,
                    "title": h.get("title", "") or "",
                    "content": (h.get("content", "") or "")[:900],
                    "rdate": h.get("release_date", "") or "",
                    "files": h.get("files", 0) or 0,
                }
        print(f"  arrayexpress page {page}: {len(out)} total", flush=True)
        if len(hits) < 100:
            break
    return list(out.values())


# --------------------------------------------------------------------------
# Classification vocabularies
# --------------------------------------------------------------------------
def _c(patterns):
    return [(name, re.compile(rx, re.I)) for name, rx in patterns]


IMAGING = _c(
    [
        ("Xenium", r"\bxenium\b"),
        ("MERFISH / MERSCOPE", r"\bmerfish\b|\bmerscope\b|\bvizgen\b"),
        ("CosMx", r"\bcosmx\b|spatial molecular imager"),
        ("seqFISH", r"\bseqfish\b"),
        ("STARmap", r"\bstarmap\b"),
        ("In situ sequencing", r"\bin situ sequencing\b|\bhybiss\b|\bcartana\b"),
        ("IMC / MIBI", r"\bmibi-?tof\b|imaging mass cytometry"),
        ("Molecular Cartography", r"molecular cartography|resolve bioscience"),
    ]
)

SEQUENCING = _c(
    [
        ("Visium HD", r"visium\s*hd"),
        ("Visium", r"\bvisium\b"),
        ("GeoMx DSP", r"\bgeomx\b|digital spatial profil"),
        ("Slide-seq", r"\bslide-?seq\b"),
        ("Stereo-seq", r"\bstereo-?seq\b|\bstereoseq\b"),
        ("Slide-tags", r"\bslide-?tags\b"),
        ("DBiT-seq", r"\bdbit-?seq\b"),
        ("Curio Seeker", r"\bcurio seeker\b"),
        ("Open-ST", r"\bopen-?st\b"),
        ("Tomo-seq", r"\btomo-?seq\b"),
        # Plenty of studies never write "spatial transcriptomics" out — they say
        # "ST" and nothing else, which is how a spatial series ends up looking
        # like a single-cell one. A bare \bST\b is far too noisy (ST-EPN =
        # supratentorial ependymoma, ST segment, ST2), so ST is only accepted
        # next to a word that makes it unambiguous.
        (
            "Spatial (unspecified)",
            r"spatial transcriptom|spatially resolved|spatial gene expression|spatial rna"
            r"|\[st\]|\(st\)"
            r"|\bst[- ](data|slides?|spots?|sections?|samples?|experiments?|profil\w+|analys\w+|arrays?|methods?|platforms?)\b"
            r"|\b(scrna-?seq|snrna-?seq|single[- ]cell)\s+and\s+st\b"
            r"|\bst\s+and\s+(scrna|snrna|single[- ]cell)",
        ),
    ]
)

SC_RE = re.compile(
    r"single[- ]cell rna|scrna-?seq|snrna-?seq|single[- ]nucleus rna|"
    r"single[- ]cell transcriptom|10x genomics chromium|\bchromium\b|cite-?seq|"
    r"single[- ]cell sequencing|smart-?seq2?|single[- ]cell atlas",
    re.I,
)

SC_PLATFORM = _c(
    [
        # "10x Genomics Visium" / "10x Genomics Xenium" name the spatial product, not Chromium
        ("10x Chromium", r"\bchromium\b|10x genomics(?![\s\W]*(visium|xenium))"),
        ("Smart-seq2", r"smart-?seq2?"),
        ("snRNA-seq", r"snrna-?seq|single[- ]nucleus"),
        ("BD Rhapsody", r"bd rhapsody"),
        ("Drop-seq", r"drop-?seq"),
        ("CITE-seq", r"cite-?seq"),
    ]
)

# Cancer subtypes. A record can carry several (LUAD *and* EGFR-mutant), so this
# is matched as a multi-value field rather than a single winner.
#
# Every pattern here must be unambiguous IN THIS CORPUS. Short abbreviations are
# the trap: ILC is "innate lymphoid cell" far more often than "invasive lobular
# carcinoma", FL and MCL collide with gene and cell-type names, IDC is also
# "interdigitating dendritic cell", BCC is also "B-cell ...". Those are matched
# by their full phrase only. Add a bare abbreviation here only after checking it
# against the catalog for collisions.
SUBTYPES = _c(
    [
        # lung
        ("LUAD", r"\bluad\b|lung adenocarcinoma|adenocarcinoma of (the )?lungs?|pulmonary adenocarcinoma"),
        ("LUSC", r"\blusc\b|lung squamous|squamous cell (carcinoma )?of (the )?lungs?|pulmonary squamous"),
        ("LCNEC / large cell", r"\blcnec\b|large[- ]cell (neuroendocrine )?carcinoma"),
        ("Adenosquamous", r"adenosquamous"),
        ("EGFR-mutant", r"\begfr[- ](mutant|mutated|driven|altered)|egfr mutation"),
        ("ALK-rearranged", r"\balk[- ](rearrang|fusion|positive|translocat)"),
        ("KRAS-mutant", r"\bkras[- ]?(g12|g13|mutant|mutated|mutation)"),
        # breast
        ("TNBC", r"\btnbc\b|triple[- ]negative breast"),
        ("HER2-positive", r"her2[- ]?(positive|enriched|amplified)|erbb2[- ]amplif"),
        ("Luminal", r"luminal [ab]\b|luminal subtype|luminal breast"),
        ("Invasive ductal (IDC)", r"invasive ductal carcinoma"),
        ("Invasive lobular (ILC)", r"invasive lobular carcinoma"),
        ("DCIS", r"\bdcis\b|ductal carcinoma in situ"),
        # kidney
        ("ccRCC", r"\bccrcc\b|clear[- ]cell renal"),
        ("pRCC", r"\bprcc\b|papillary renal"),
        ("chRCC", r"\bchrcc\b|chromophobe renal"),
        # prostate
        ("CRPC", r"\bm?crpc\b|castration[- ]resistant"),
        ("NEPC", r"\bnepc\b|neuroendocrine prostate"),
        # brain
        ("IDH-mutant", r"\bidh[12]?[- ]?(mutant|mutated)"),
        ("IDH-wildtype", r"\bidh[12]?[- ]?(wild[- ]?type|wt)\b"),
        ("DIPG / diffuse midline", r"\bdipg\b|diffuse (intrinsic pontine|midline) glioma"),
        # blood
        ("AML", r"\baml\b|acute myeloid leuk"),
        ("B-ALL", r"\bb-all\b|b[- ]cell acute lymphoblastic"),
        ("T-ALL", r"\bt-all\b|t[- ]cell acute lymphoblastic"),
        ("CLL", r"\bcll\b|chronic lymphocytic leuk"),
        ("CML", r"\bcml\b|chronic myeloid leuk"),
        ("MDS", r"\bmds\b|myelodysplastic"),
        ("DLBCL", r"\bdlbcl\b|diffuse large b[- ]cell lymphoma"),
        ("Follicular lymphoma", r"follicular lymphoma"),
        ("Mantle cell lymphoma", r"mantle cell lymphoma"),
        ("Hodgkin", r"hodgkin"),
        # GI
        ("ESCC", r"\bescc\b|(o)?esophageal squamous"),
        ("Esophageal adenocarcinoma", r"(o)?esophageal adenocarcinoma"),
        ("iCCA", r"\bicca\b|intrahepatic cholangiocarcinoma"),
        ("eCCA", r"\becca\b|extrahepatic cholangiocarcinoma|perihilar cholangio"),
        ("MSI-high", r"\bmsi-?h\b|microsatellite instab"),
        ("MSS", r"microsatellite stable"),
        ("Signet ring", r"signet[- ]ring"),
        # pancreas / liver
        ("PDAC", r"\bpdac\b|pancreatic ductal adenocarcinoma"),
        ("HCC", r"\bhcc\b|hepatocellular carcinoma"),
        # bladder
        ("MIBC", r"\bmibc\b|muscle[- ]invasive bladder"),
        ("NMIBC", r"\bnmibc\b|non[- ]muscle[- ]invasive bladder"),
        # gynaecological
        ("HGSOC", r"\bhgsoc\b|high[- ]grade serous"),
        ("Ovarian clear cell", r"ovarian clear cell"),
        # head & neck / skin
        ("HNSCC", r"\bhnscc\b|head and neck squamous"),
        ("OSCC", r"\boscc\b|oral squamous cell"),
        ("HPV-positive", r"hpv[- ]?(positive|associated|driven)"),
        ("Basal cell carcinoma", r"basal cell carcinoma"),
        ("cSCC", r"\bcscc\b|cutaneous squamous cell"),
    ]
)

CANCER = _c(
    [
        ("NSCLC", r"non[-\u2010\u2013\s]?small[-\u2010\u2013\s]cell lung|\bnsclc\b|lung adenocarcinoma|\bluad\b|lung squamous|\blusc\b"
                  r"|adenocarcinoma of (the )?lungs?|pulmonary adenocarcinoma|\blcnec\b|bronchioloalveolar"
                  r"|lepidic|adenosquamous carcinoma of (the )?lung|large[- ]cell carcinoma of (the )?lung"),
        # the lookbehinds keep "non-small cell lung cancer" (any dash or a space)
        # from also reading as small cell lung cancer
        ("SCLC", r"(?<!non-)(?<!non\u2013)(?<!non\u2010)(?<!non )small[- ]cell lung cancer|\bsclc\b"),
        # Mesothelioma is a pleural/peritoneal cancer, not a lung cancer; it used to
        # fall into "Lung cancer (other)" and from there into lung searches.
        ("Mesothelioma", r"mesothelioma"),
        ("Lung cancer (other)", r"lung (cancer|tumou?r|carcinoma|metasta|neoplas)|pulmonary (cancer|carcinoma)"),
        ("Breast cancer", r"breast (cancer|tumou?r|carcinoma)|triple[- ]negative breast|\btnbc\b"
                          r"|ductal carcinoma in situ|\bdcis\b|invasive (ductal|lobular) carcinoma"
                          r"|luminal [ab] (breast|subtype)"),
        ("Colorectal cancer", r"colorectal|\bcrc\b|colon (cancer|adenocarcinoma|tumou?r)|rectal cancer"),
        ("Glioma / GBM", r"glioblastoma|\bgbm\b|\bglioma\b|astrocytoma|oligodendroglioma|medulloblastoma"
                         r"|ependymoma|\bdipg\b|diffuse (intrinsic pontine|midline) glioma"),
        ("Meningioma", r"meningioma"),
        ("Pancreatic cancer", r"pancreatic (cancer|ductal|adenocarcinoma|tumou?r)|\bpdac\b"),
        ("Prostate cancer", r"prostate (cancer|tumou?r|adenocarcinoma)|\bm?crpc\b|\bnepc\b|castration[- ]resistant"),
        ("Melanoma", r"melanoma"),
        ("Liver / HCC", r"hepatocellular|\bhcc\b|liver (cancer|tumou?r)|cholangiocarcinoma|\bicca\b|\becca\b"),
        ("Gastric cancer", r"gastric (cancer|adenocarcinoma|tumou?r)|stomach cancer"),
        ("Ovarian cancer", r"ovarian (cancer|carcinoma|tumou?r|clear cell)|\bhgsoc\b"),
        ("Kidney / RCC", r"renal cell carcinoma|\brcc\b|\bccrcc\b|\bprcc\b|\bchrcc\b|kidney (cancer|tumou?r)|wilms"),
        ("Head & neck cancer", r"head and neck|\bhnscc\b|oral squamous|\boscc\b|nasopharyngeal|laryngeal cancer"),
        ("Lymphoma", r"lymphoma|\bdlbcl\b"),
        ("Leukemia", r"leukemia|leukaemia|\baml\b|\bcll\b|\bcml\b|\bapl\b|\bb-all\b|\bt-all\b|myelodysplas"),
        ("Multiple myeloma", r"multiple myeloma|plasma cell myeloma"),
        ("Esophageal cancer", r"(o)?esophageal (cancer|carcinoma|adenocarcinoma|squamous)|\bescc\b"),
        ("Bladder cancer", r"bladder cancer|urothelial (carcinoma|cancer)|\bmibc\b|\bnmibc\b"),
        ("Cervical cancer", r"cervical (cancer|carcinoma)"),
        ("Endometrial cancer", r"endometrial (cancer|carcinoma)|uterine (cancer|carcinoma)"),
        ("Sarcoma", r"sarcoma|\bgist\b"),
        ("Neuroblastoma", r"neuroblastoma"),
        ("Thyroid cancer", r"thyroid (cancer|carcinoma)"),
        ("Skin cancer (non-melanoma)", r"basal cell carcinoma|cutaneous squamous|\bcscc\b"),
        ("Cancer (other / pan-cancer)", r"\bcancers?\b|\btumou?rs?\b|carcinoma|malignan|metasta|neoplas|oncogen|oncolog"),
    ]
)

OTHER_DISEASE = _c(
    [
        ("Neurodegenerative", r"alzheimer|parkinson|\bals\b|amyotrophic lateral|huntington|multiple sclerosis|dementia|neurodegener|tauopath"),
        ("Fibrosis", r"\bfibrosis\b|fibrotic|\bipf\b|idiopathic pulmonary fibrosis|cirrhosis"),
        ("Autoimmune / inflammatory", r"rheumatoid|\blupus\b|psoriasis|crohn|ulcerative colitis|inflammatory bowel|\bibd\b|autoimmun|atopic dermatitis|\bsle\b"),
        ("Cardiovascular", r"myocardial infarction|heart failure|cardiomyopath|atheroscleros|cardiac (injury|disease|remodel)|aortic aneurysm"),
        ("Infectious disease", r"covid|sars-cov-2|tuberculosis|influenza|\bhiv\b|\bsepsis\b|malaria|helicobacter|viral infection"),
        ("Metabolic", r"diabet|obesity|\bnafld\b|\bnash\b|steatohepatitis|metabolic syndrome|insulin resistance"),
        ("Kidney disease", r"chronic kidney disease|renal fibrosis|nephropath|glomerulo"),
        ("Respiratory (non-cancer)", r"\bcopd\b|\basthma\b|cystic fibrosis|\bards\b|emphysema|bronchopulmonary"),
        ("Injury / regeneration", r"\bwound\b|regenerat|\binjury\b|\bischemi"),
        ("Development / embryo", r"\bembryo|organogenes|\bfetal\b|\bfoetal\b|gestation|developmental atlas"),
        ("Ageing", r"\bageing\b|\baging\b|senescen"),
    ]
)

HEALTHY = re.compile(
    r"healthy (donor|tissue|control|human|adult)|normal (tissue|human|donor)|"
    r"reference (atlas|map)|cell atlas",
    re.I,
)

TISSUE = _c(
    [
        ("Lung", r"\blungs?\b|pulmonary|bronch|alveol|\bairway|\bnsclc\b|\bsclc\b|\bluad\b|\blusc\b"),
        ("Brain / CNS", r"\bbrain\b|cortex|cortical|hippocamp|cerebell|\bneuron|\bcns\b|spinal cord|striatum|\bglia"),
        ("Breast", r"\bbreast\b|mammary"),
        ("Colon / intestine", r"\bcolon\b|colorect|intestin|\bileum\b|\brectal\b|duoden|\bcolitis\b"),
        ("Pancreas", r"pancrea|\bislet"),
        ("Liver", r"\bliver\b|hepatic|hepatocyt|\bhcc\b"),
        ("Kidney", r"\bkidney\b|\brenal\b|nephron|glomerul"),
        ("Skin", r"\bskin\b|epiderm|dermal|cutaneous|keratinocyt"),
        ("Prostate", r"prostat"),
        ("Stomach", r"\bstomach\b|\bgastric\b"),
        ("Ovary / uterus / placenta", r"\bovar|uter|endometri|\bcervix\b|\bcervical\b|placenta|\bdecidua"),
        ("Lymph node / spleen", r"lymph node|\bspleen\b|\btonsil|lymphoid"),
        ("Bone marrow / blood", r"bone marrow|\bblood\b|\bpbmc\b|h(a)?ematopoie|\bmyeloid\b"),
        ("Heart", r"\bheart\b|cardiac|myocard|ventric|cardiomyocyt"),
        ("Bladder / urinary", r"\bbladder\b|urothel|\burinary\b"),
        ("Muscle", r"skeletal muscle|\bmyofib|\bmuscles?\b"),
        ("Eye / retina", r"\bretina|\beyes?\b|\bcornea"),
        ("Bone / cartilage", r"\bbones?\b|cartilage|osteo|chondrocyt"),
        ("Testis", r"\btest(is|es|icular)\b|\bsperm"),
        ("Thymus", r"thymus|thymic"),
        ("Adrenal", r"adrenal"),
        ("Thyroid", r"thyroid"),
        ("Esophagus", r"esophag|oesophag"),
        ("Head & neck", r"\boral\b|\btongue\b|salivary|nasopharyn|laryn|\bhnscc\b"),
        ("Adipose", r"adipose"),
        ("Embryo / fetal", r"\bembryo|\bfetal\b|\bfoetal\b"),
    ]
)

# Phrases that contain a cancer keyword but say nothing about the study being
# about cancer. They are blanked out before disease matching, otherwise every
# TNF / TNFR immunology study lands in the pan-cancer bucket.
DISEASE_NOISE = re.compile(
    r"tumou?r necrosis factor|\btnf-?(alpha|a|r\d?)\b|\btnfrsf\d*\b|"
    r"tumou?r necrosis|anti-?tnf",
    re.I,
)

# A study whose TITLE says it spans many cancer types is pan-cancer, even when
# its abstract names lung adenocarcinoma first (KCNQ3 pan-cancer analysis,
# the 30-cancer-type tumour/normal atlas). Title only: abstracts mention
# "multiple cancer types" in passing far too often to be trusted.
PAN_CANCER_TITLE = re.compile(
    r"pan-?cancer|\b(\d+|multiple|several|diverse|various|twenty|thirty)\s+(different\s+)?(human\s+)?"
    r"(cancer|tumou?r) types|across (\w+\s+){0,3}(cancer|tumou?r) types|thousand tumou?rs",
    re.I,
)

# Method / tool development studies. Their data are usually a demonstration
# set (NEST, DGAT ...), so they are labelled rather than dropped, and the site
# can hide them. Title "ACRONYM: ... <tool word>" or "Benchmarking ...", or an
# abstract that says "we develop/present a <computational> method/tool/...".
_TOOL_WORDS = (
    r"(network|method|framework|tool|toolkit|algorithm|pipeline|learning|inferr?\w*|predict\w*"
    r"|deconvol\w*|benchmark\w*|computational|software|package|embedding|transformer|graph"
    r"|imputation|integration of)"
)
METHOD_TITLE = re.compile(
    r"^\s*[A-Za-z][\w\-\.\+]{1,24}\s*:\s.*\b" + _TOOL_WORDS
    + r"|^\s*(benchmark\w*|a (computational|deep[- ]learning|statistical) (method|framework|approach))",
    re.I,
)
METHOD_BODY = re.compile(
    r"\bwe (develop|developed|present|introduce|propose|proposed|describe|built|build)\w*\s+(a |an )?"
    r"(new |novel |general |scalable |unified |open[- ]source )*"
    r"(computational |statistical |deep[- ]learning |machine[- ]learning |graph[- ]based |bayesian "
    r"|probabilistic |ai[- ]driven |generative )*"
    r"(method|tool|toolkit|framework|algorithm|pipeline|software|package)\b",
    re.I,
)


def study_type(title, body):
    if METHOD_TITLE.search(title or "") or METHOD_BODY.search(body or ""):
        return "Method / tool development"
    return "Biological study"


# GEO subseries usually end their title with the assay: "... [scRNA-seq]",
# "... [Xenium in situ]", "... [bulkRNA-seq_Tregs]". The abstract is shared by
# every subseries and mentions all assays, so it cannot tell them apart; the
# tag can. A single-cell subseries of a spatial study is single-cell (and paired
# only through its sibling), and a bulk / ChIP / ATAC subseries is neither.
TITLE_TAG = re.compile(r"[\[\(]([^\[\]\(\)]{2,60})[\]\)]\s*$")
TAG_SPATIAL = re.compile(
    r"\bst\b|spatial|visium|xenium|cosmx|merfish|merscope|geomx|stereo-?seq|slide-?seq|seqfish"
    r"|starmap|\bcodex\b|\bimc\b|in situ|dbit|curio|molecular cartography",
    re.I,
)
TAG_SC = re.compile(
    r"sc-?rna|sn-?rna|single[- _]?(cell|nucle)|cite-?seq|chromium|multiome|sc-?atac|smart-?seq"
    r"|drop-?seq|rhapsody|sc-?tcr|sc-?bcr|\bvdj\b|10x\s*(3|5)",
    re.I,
)
TAG_OTHER = re.compile(
    r"bulk|chip-?seq|cut ?(&|and|n) ?(run|tag)|cutntag|atac|\bwgs\b|\bwgbs\b|\brrbs\b|\bwes\b|exome"
    r"|methyl|bisulfite|proteom|mass spec|crispr|microarray|\bacgh\b|\bsnp\b|hi-?c\b|ribo-?seq"
    r"|\b16s\b|tcr-?seq|bcr-?seq|merip|amplicon|nascent|total ?rna|mirna|small ?rna|\brip\b"
    r"|dart-?seq|ss-?rna|chromatin",
    re.I,
)
# What a tag must also say for an ATAC or immune-repertoire library to count as
# a transcriptome: "[scATAC-seq]" and "[scTCR-seq]" are not, "[Multiome]" and
# "[scRNA-seq and scTCR-seq]" are.
TAG_HAS_RNA = re.compile(r"rna|multiome|gex|gene expression|transcript|cite-?seq", re.I)
TAG_RNASEQ = re.compile(r"rna[- _]?seq|rna sequencing", re.I)


def title_assay(title):
    """'spatial' | 'singlecell' | 'rnaseq' | 'other' | None, from the title.

    Mostly from the trailing [assay] tag; a title that is about chromatin
    accessibility only ("scATAC-seq of first and second-trimester human gonads")
    is 'other' even without a tag.
    """
    m = TITLE_TAG.search(title or "")
    assay = _tag_assay(m.group(1)) if m else None
    if assay is None and _atac_only(title or ""):
        return "other"
    return assay


def _atac_only(text):
    return bool(
        re.search(r"atac", text, re.I)
        and not TAG_HAS_RNA.search(text)
        and not TAG_SPATIAL.search(text)
        and not re.search(r"transcriptom|expression|multiom", text, re.I)
    )


def _tag_assay(tag):
    sp, sc = bool(TAG_SPATIAL.search(tag)), bool(TAG_SC.search(tag))
    if sp and sc:
        return None  # "[scRNA-seq and ST]": the record itself is paired
    if sp:
        return "spatial"
    if re.search(r"atac|tcr|bcr|vdj|repertoire", tag, re.I) and not TAG_HAS_RNA.search(tag):
        return "other"  # chromatin or receptor sequences only, no expression
    if sc:
        return "singlecell"
    if TAG_RNASEQ.search(tag) and not TAG_OTHER.search(tag):
        # "[RNA-seq]", "[RNA-Seq PDX]": GEO's usual name for the bulk companion
        # of a [scRNA-seq] or [ST] subseries, but GSE118389 is a Smart-seq2
        # single-cell series tagged just "[RNA-Seq]". finalize() keeps it when
        # it has no sibling, or when every sibling is a non-expression assay.
        return "rnaseq"
    if TAG_OTHER.search(tag):
        return "other"
    return None


def linkable(row):
    """Records that are about to be dropped must not size or join study groups:
    a 10-series study whose ChIP and ATAC parts are counted looks like
    boilerplate (groups are capped), while without them it is a real study."""
    return title_assay(row["title"]) != "other"


SC_PLATFORM_NAMES = {name for name, _ in SC_PLATFORM} | {"scRNA-seq (unspecified)"}
IMAGING_NAMES = {name for name, _ in IMAGING}
UNSPECIFIED_SPATIAL = "Spatial (unspecified)"


def _platform_in(text):
    img = first_match(IMAGING, text)
    if img:
        return img, "imaging"
    seq = first_match(SEQUENCING, text)
    if seq == UNSPECIFIED_SPATIAL:
        # Spatial, but the text never says how it was made: claiming
        # "sequencing-based" here would put it under the wrong filter.
        return seq, "unspecified"
    return (seq, "sequencing") if seq else (None, None)


def detect_platform(title, text):
    """The spatial platform of a record, preferring what its own title says.

    The abstract of a multi-assay study names every platform used, and because
    imaging platforms are checked first, a "[Visium]" subseries of a study that
    also ran Xenium used to be labelled Xenium. The trailing [assay] tag, then
    the title, then the text decide; a named platform beats "unspecified".
    """
    m = TITLE_TAG.search(title or "")
    fallback = (None, None)
    for src in ([m.group(1)] if m else []) + [title or "", text or ""]:
        found, cls = _platform_in(src)
        if found and found != UNSPECIFIED_SPATIAL:
            return found, cls
        if found and not fallback[0]:
            fallback = (found, cls)
    return fallback


def platform_source(title, platform):
    """Where a spatial platform label came from: 'title', 'abstract' or ''."""
    if not platform or platform == UNSPECIFIED_SPATIAL or platform in SC_PLATFORM_NAMES:
        return ""
    named, _ = detect_platform(title, "")
    return "title" if named == platform else "abstract"


def platform_class_of(platform):
    if platform in IMAGING_NAMES:
        return "imaging"
    return "unspecified" if platform == UNSPECIFIED_SPATIAL else "sequencing"


# Organ-specific subtypes. A breast-cancer study whose abstract mentions lung
# adenocarcinoma in passing must not answer a LUAD search. Subtypes not listed
# here (EGFR-mutant, KRAS-mutant, HER2-positive, IDH, MSI, HPV ...) occur in
# several cancers and are kept whatever the disease.
SUBTYPE_DISEASES = {
    "LUAD": {"NSCLC", "Lung cancer (other)"},
    "LUSC": {"NSCLC", "Lung cancer (other)"},
    "LCNEC / large cell": {"NSCLC", "SCLC", "Lung cancer (other)"},
    "TNBC": {"Breast cancer"},
    "Luminal": {"Breast cancer"},
    "Invasive ductal (IDC)": {"Breast cancer"},
    "Invasive lobular (ILC)": {"Breast cancer"},
    "DCIS": {"Breast cancer"},
    "ccRCC": {"Kidney / RCC"},
    "pRCC": {"Kidney / RCC"},
    "chRCC": {"Kidney / RCC"},
    "CRPC": {"Prostate cancer"},
    "NEPC": {"Prostate cancer"},
    "DIPG / diffuse midline": {"Glioma / GBM"},
    "AML": {"Leukemia"},
    "B-ALL": {"Leukemia"},
    "T-ALL": {"Leukemia"},
    "CLL": {"Leukemia", "Lymphoma"},
    "CML": {"Leukemia"},
    "MDS": {"Leukemia"},
    "DLBCL": {"Lymphoma"},
    "Follicular lymphoma": {"Lymphoma"},
    "Mantle cell lymphoma": {"Lymphoma"},
    "Hodgkin": {"Lymphoma"},
    "ESCC": {"Esophageal cancer", "Head & neck cancer"},
    "Esophageal adenocarcinoma": {"Esophageal cancer", "Gastric cancer"},
    "iCCA": {"Liver / HCC"},
    "eCCA": {"Liver / HCC"},
    "PDAC": {"Pancreatic cancer"},
    "HCC": {"Liver / HCC"},
    "MIBC": {"Bladder cancer"},
    "NMIBC": {"Bladder cancer"},
    "HGSOC": {"Ovarian cancer"},
    "Ovarian clear cell": {"Ovarian cancer"},
    "HNSCC": {"Head & neck cancer"},
    "OSCC": {"Head & neck cancer"},
    "Basal cell carcinoma": {"Skin cancer (non-melanoma)"},
    "cSCC": {"Skin cancer (non-melanoma)"},
}


def compatible_subtypes(subtypes, disease):
    if disease == "Cancer (other / pan-cancer)":
        return list(subtypes)
    return [s for s in subtypes if s not in SUBTYPE_DISEASES or disease in SUBTYPE_DISEASES[s]]

# When the text names no organ, the disease usually does.
DISEASE_TISSUE = {
    "NSCLC": "Lung",
    "SCLC": "Lung",
    "Lung cancer (other)": "Lung",
    "Mesothelioma": "Pleura / peritoneum",
    "Glioma / GBM": "Brain / CNS",
    "Meningioma": "Brain / CNS",
    "Melanoma": "Skin",
    "Leukemia": "Bone marrow / blood",
    "Multiple myeloma": "Bone marrow / blood",
    "Lymphoma": "Lymph node / spleen",
    "Breast cancer": "Breast",
    "Colorectal cancer": "Colon / intestine",
    "Pancreatic cancer": "Pancreas",
    "Prostate cancer": "Prostate",
    "Liver / HCC": "Liver",
    "Gastric cancer": "Stomach",
    "Ovarian cancer": "Ovary / uterus / placenta",
    "Kidney / RCC": "Kidney",
    "Head & neck cancer": "Head & neck",
    "Esophageal cancer": "Esophagus",
    "Bladder cancer": "Bladder / urinary",
    "Thyroid cancer": "Thyroid",
}

GENERIC_CANCER = "Cancer (other / pan-cancer)"


CELLS_RE = re.compile(
    r"([\d]{1,3}(?:[,\.]\d{3})+|\d{4,9})\s*(?:single[- ])?(?:cells|nuclei)\b", re.I
)
CONTROLLED_RE = re.compile(r"\bdbgap\b|\bega\b|controlled access|restricted access", re.I)


def first_match(vocab, text):
    for name, rx in vocab:
        if rx.search(text):
            return name
    return None


def all_matches(vocab, text, limit=6):
    """Every label whose pattern hits — for fields a record can hold several of."""
    out = []
    for name, rx in vocab:
        if rx.search(text):
            out.append(name)
            if len(out) >= limit:
                break
    return out


def best_tissue(text):
    """Pick the tissue with the most keyword hits, not merely the first listed."""
    best, best_score = None, 0
    for name, rx in TISSUE:
        score = len(rx.findall(text))
        if score > best_score:
            best, best_score = name, score
    return best or "Other"


def cell_count(text):
    m = CELLS_RE.search(text)
    if not m:
        return 0
    try:
        v = int(m.group(1).replace(",", "").replace(".", ""))
    except ValueError:
        return 0
    return v if 500 <= v <= 50_000_000 else 0


def organism(taxon_text, fallback_text=""):
    t = taxon_text or fallback_text
    human = bool(re.search(r"Homo sapiens", t, re.I))
    mouse = bool(re.search(r"Mus musculus", t, re.I))
    if human and mouse:
        return "Human & Mouse"
    if human:
        return "Human"
    if mouse:
        return "Mouse"
    return (taxon_text or "Other")[:40] or "Other"


def squash(s, n):
    return re.sub(r"\s+", " ", s or "").strip()[:n]


def classify(title, body, extra=""):
    """Return the derived fields for one record."""
    text = f"{title} || {body} || {extra}"
    platform, platform_class = detect_platform(title, text)
    has_sc = bool(SC_RE.search(text))
    disease_text = DISEASE_NOISE.sub(" ", text)
    cancer = first_match(CANCER, disease_text)
    # The title says what the study is about; the abstract often names other
    # cancers in passing (an EGFR-inhibitor TNBC study that cites lung cancer
    # used to land in NSCLC because NSCLC is first in the list).
    title_cancer = first_match(CANCER, DISEASE_NOISE.sub(" ", title or ""))
    if title_overrides(title_cancer, cancer):
        cancer = title_cancer
    if cancer and PAN_CANCER_TITLE.search(title or ""):
        cancer = "Cancer (other / pan-cancer)"
    other = None if cancer else first_match(OTHER_DISEASE, disease_text)
    if not cancer and not other and HEALTHY.search(text):
        other = "Healthy reference atlas"
    return {
        "platform": platform,
        "platform_class": platform_class,
        "platform_source": platform_source(title, platform),
        "sc_platform": first_match(SC_PLATFORM, text) if has_sc else None,
        "has_sc": has_sc,
        "disease": cancer or other or "Unspecified",
        "disease_group": "Cancer" if cancer else ("Other disease" if other else "Unspecified"),
        # Subtypes are only meaningful once a record reads as cancer; leaving the
        # gate off would tag immunology studies with blood-cancer abbreviations.
        "subtypes": compatible_subtypes(all_matches(SUBTYPES, disease_text), cancer) if cancer else [],
        "tissue": tissue_for(best_tissue(text), cancer or other),
        "cells": cell_count(text),
        "access": "Controlled" if CONTROLLED_RE.search(text) else "Open",
        "study_type": study_type(title, body),
    }


LUNG_SPECIFIC = {"NSCLC", "SCLC", "Mesothelioma"}


def title_overrides(title_cancer, body_cancer):
    """Should the cancer named in the title replace the one from the text?"""
    if not title_cancer or title_cancer == GENERIC_CANCER:
        return False
    # "lung cancer" in the title is less specific than NSCLC in the abstract
    if title_cancer == "Lung cancer (other)" and body_cancer in LUNG_SPECIFIC:
        return False
    return True


def tissue_for(tissue, disease):
    if tissue == "Other" and disease in DISEASE_TISSUE:
        return DISEASE_TISSUE[disease]
    return tissue


def apply_title_assay(row, assay):
    """Make a subseries carry its own assay, not the study's. False = drop it."""
    if assay == "other":
        return False
    if assay == "singlecell":
        row["modality"] = "singlecell"
        m = TITLE_TAG.search(row["title"])
        tag_platform = first_match(SC_PLATFORM, m.group(1)) if m else None
        if tag_platform:
            row["platform"] = row["sc_platform"] = tag_platform
        elif row["platform"] not in SC_PLATFORM_NAMES:
            row["platform"] = row.get("sc_platform") or "scRNA-seq (unspecified)"
        row["platform_class"] = "sequencing"
        row["platform_source"] = ""
    elif assay == "spatial":
        row["modality"] = "spatial"
        if row["platform"] in SC_PLATFORM_NAMES:
            row["platform"] = UNSPECIFIED_SPATIAL
            row["platform_class"] = "unspecified"
            row["platform_source"] = ""
    return True


# --------------------------------------------------------------------------
# Queries
# --------------------------------------------------------------------------
SPATIAL_QUERY = (
    '("spatial transcriptomics"[All Fields] OR "spatially resolved transcriptomics"[All Fields] '
    'OR Visium[All Fields] OR Xenium[All Fields] OR MERFISH[All Fields] OR MERSCOPE[All Fields] '
    'OR CosMx[All Fields] OR GeoMx[All Fields] OR "Slide-seq"[All Fields] OR "Stereo-seq"[All Fields] '
    'OR seqFISH[All Fields] OR STARmap[All Fields] OR "DBiT-seq"[All Fields] OR "Slide-tags"[All Fields] '
    'OR "spatial gene expression"[All Fields]) AND gse[Entry Type] '
    'AND ("Homo sapiens"[Organism] OR "Mus musculus"[Organism])'
)

SC_CANCER_QUERY = (
    '("single cell RNA"[All Fields] OR "scRNA-seq"[All Fields] '
    'OR "single-cell RNA sequencing"[All Fields] OR "snRNA-seq"[All Fields]) '
    'AND (cancer[All Fields] OR tumor[All Fields] OR tumour[All Fields] '
    'OR carcinoma[All Fields] OR malignant[All Fields]) AND gse[Entry Type] '
    'AND "Homo sapiens"[Organism]'
)

# Studies that only ever write "ST" are invisible to SPATIAL_QUERY. A bare "ST"
# term is useless in E-utilities, so search the phrases instead.
ST_QUERY = (
    '("ST data"[All Fields] OR "ST slides"[All Fields] OR "ST spots"[All Fields] '
    'OR "ST sections"[All Fields] OR "ST profiling"[All Fields] '
    'OR "ST analysis"[All Fields] OR "ST experiments"[All Fields]) '
    'AND gse[Entry Type] AND ("Homo sapiens"[Organism] OR "Mus musculus"[Organism])'
)

AE_SPATIAL_QUERY = (
    'spatial transcriptomics OR Visium OR Xenium OR MERFISH OR CosMx OR GeoMx '
    'OR "Slide-seq" OR "Stereo-seq"'
)
AE_SC_QUERY = 'single cell RNA-seq AND (cancer OR tumour OR carcinoma)'


def group_key(text):
    """Key for spotting subseries of one study.

    GEO splits many studies into a SuperSeries plus one subseries per assay —
    GSE189357 holds the scRNA-seq half and GSE189487 the spatial half of the
    same LUAD atlas. esummary does not expose that link (the `relations` field
    comes back empty for both), but the subseries share the study abstract
    verbatim, so the normalised abstract identifies the group.
    """
    t = re.sub(r"\s+", " ", (text or "")).strip().lower()
    return t[:300] if len(t) >= 120 else ""


# ArrayExpress companion studies do not share an abstract: E-MTAB-13526
# ("Single cell RNA-seq atlas of human NSCLC lesions and non-involved tissue")
# and E-MTAB-13530 ("Spatial transcriptomics (10x Visium) of human NSCLC
# lesions and non-involved tissue") describe their own assay. What they share
# is the title once the assay words are removed, and accessions a few numbers
# apart (they are submitted together).
_AE_ASSAY = re.compile(
    r"\([^)]*\)|\[[^\]]*\]|single[- ](cell|nucleus|nuclei)|sc-?rna-?seq|sn-?rna-?seq|rna-?seq"
    r"|spatial(ly)?( resolved)?|transcriptom\w*|visium( hd)?|xenium|cosmx|merfish|geomx|stereo-?seq"
    r"|slide-?seq|10x( genomics)?|chromium|sequencing|profiling|atlas|dataset|data|analysis"
    r"|\b(of|the|and|in|from|by|with|a|an|human|mouse|samples?|tissues?)\b|gene expression",
    re.I,
)
AE_TITLE_SIMILARITY = 0.7
AE_ACCESSION_WINDOW = 200


def _ae_title_core(title):
    words = re.findall(r"[a-z0-9\-]+", _AE_ASSAY.sub(" ", (title or "").lower()))
    return {w for w in words if len(w) > 1}


def _ae_number(acc):
    m = re.match(r"^E-([A-Z]+)-(\d+)$", acc or "")
    return (m.group(1), int(m.group(2))) if m else (None, None)


def ae_title_links(rows):
    """Pairs of ArrayExpress accessions that look like halves of one study."""
    ae = []
    for r in rows:
        if r.get("db") != "ArrayExpress" or not linkable(r):
            continue
        prefix, num = _ae_number(r["accession"])
        core = _ae_title_core(r["title"])
        if prefix and len(core) >= 3:
            ae.append((prefix, num, core, r["accession"]))
    ae.sort()
    links = []
    for i, (p1, n1, c1, a1) in enumerate(ae):
        for p2, n2, c2, a2 in ae[i + 1 :]:
            if p2 != p1 or n2 - n1 > AE_ACCESSION_WINDOW:
                break
            if len(c1 & c2) / len(c1 | c2) >= AE_TITLE_SIMILARITY:
                links.append((a1, a2))
    return links


def geo_title_links(rows):
    """GEO subseries whose titles differ only in the trailing [assay] tag.

    Most subseries share the abstract, but some do not (GSE292469 snRNA-seq and
    GSE342824 Xenium of one osimertinib study), so the title is a second clue.
    """
    by_stem = {}
    for r in rows:
        if r.get("db") != "GEO" or not TITLE_TAG.search(r["title"]) or not linkable(r):
            continue
        stem = TITLE_TAG.sub("", r["title"]).strip().lower()
        if len(stem) >= 30:
            by_stem.setdefault(stem, []).append(r["accession"])
    links = []
    for accs in by_stem.values():
        if 2 <= len(accs) <= 8:
            links += [(accs[0], a) for a in accs[1:]]
    return links


def group_and_promote(rows, edges, max_group=12):
    """Union the edges into studies; set siblings; promote split pairs.

    A subseries that only carries scRNA-seq is still part of a paired study when
    a sibling carries the spatial half. The record keeps its own platform; only
    the modality is promoted, and `pairing` records whether the pairing comes
    from the record itself or from its siblings, so the site can say which.
    """
    present = {r["accession"] for r in rows}
    parent = {}

    def find(a):
        parent.setdefault(a, a)
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in edges:
        if a in present and b in present:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

    for r in rows:
        r["siblings"] = []
        r["pairing"] = "record" if r["modality"] == "paired" else ""

    groups = {}
    for r in rows:
        if r["accession"] in parent:
            groups.setdefault(find(r["accession"]), []).append(r)

    promoted = 0
    for members in groups.values():
        if not 2 <= len(members) <= max_group:
            continue
        accessions = [m["accession"] for m in members]
        has_spatial = any(m["modality"] in ("spatial", "paired") for m in members)
        # sc_platform of a spatial-tagged subseries comes from the abstract the
        # whole study shares, so it says nothing about what that series holds.
        has_sc = any(
            m["modality"] in ("singlecell", "paired")
            or (m["sc_platform"] and title_assay(m["title"]) != "spatial")
            for m in members
        )
        for m in members:
            m["siblings"] = [a for a in accessions if a != m["accession"]]
            if has_spatial and has_sc and m["modality"] != "paired":
                m["modality"] = "paired"
                m["pairing"] = "series"
                promoted += 1
    return promoted


def _components(accessions, edges):
    parent = {a: a for a in accessions}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in edges:
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra
    comps = {}
    for a in accessions:
        comps.setdefault(find(a), []).append(a)
    return list(comps.values())


def finalize(rows, edges):
    """Apply [assay] tags, drop non-spatial/non-single-cell subseries, link.

    Studies are grouped over every record before anything is dropped: the
    abstract edges form a star around the first subseries found, and when that
    one is the bulk companion, dropping it first would orphan its siblings.
    """
    assays = {r["accession"]: title_assay(r["title"]) for r in rows}
    comps = _components([r["accession"] for r in rows], edges)
    comp_of = {a: c for c in comps for a in c}
    kept = []
    for r in rows:
        acc = r["accession"]
        assay = assays[acc]
        if assay == "rnaseq":
            others = [assays[x] for x in comp_of[acc] if x != acc]
            # Inside a multi-series study a plain "[RNA-seq]" is the bulk
            # companion - unless every sibling is a non-expression assay
            # (GSE118389 [RNA-Seq] beside GSE118303 [WES]): then it is the
            # study's only expression data and its own text decides.
            if others and not all(o == "other" for o in others):
                continue
            assay = None
        if apply_title_assay(r, assay):
            kept.append(r)
    kept_acc = {r["accession"] for r in kept}
    kept_edges = []
    for c in comps:
        members = [a for a in c if a in kept_acc]
        kept_edges += [(members[0], a) for a in members[1:]]
    promoted = group_and_promote(kept, kept_edges)
    return kept, len(rows) - len(kept), promoted


def link_siblings(rows):
    """Edges between the parts of one study, from three kinds of evidence: a
    shared abstract (GEO subseries), a title that differs only in the [assay]
    tag (GEO), and an assay-stripped title match (ArrayExpress companions)."""
    by_key = {}
    for r in rows:
        k = r.pop("_k", "")
        if k and linkable(r):
            by_key.setdefault(k, []).append(r["accession"])
    edges = []
    for accs in by_key.values():
        if 2 <= len(accs) <= 8:  # a very large group means the key caught boilerplate
            edges += [(accs[0], a) for a in accs[1:]]
    edges += geo_title_links(rows)
    edges += ae_title_links(rows)
    return edges


# --------------------------------------------------------------------------
# Enrichment: what a series actually deposited
# --------------------------------------------------------------------------
# A study abstract names every platform the study used, and a third of spatial
# records name none at all. The files a series deposited (tissue_positions.csv,
# a Xenium transcripts.parquet, a GeoMx .dcc) and its processing text
# ("Space Ranger", "AtoMx") say what that one series holds. ArrayExpress also
# curates a study type, organism and sample count per study.
#
# Each accession is fetched once and cached in data/enrichment.json; bump
# ENRICH_VERSION when EVIDENCE changes so everything is fetched again. The
# fetch has a time budget: what does not fit is picked up on the next run.
ENRICH_VERSION = 1
ENRICH_PATH = DATA / "enrichment.json"
# 30 minutes: with the ~3-minute harvest this stays inside the workflow's
# 45-minute job limit; whatever does not fit is fetched on the next run.
ENRICH_BUDGET_S = float(os.environ.get("ENRICH_BUDGET_S", "1800"))
ENRICH_WORKERS = 2
ENRICH_PAUSE = 0.6
GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
BIOSTUDIES_API = "https://www.ebi.ac.uk/biostudies/api/v1/"

# Signatures in deposited file names and processing/protocol text. Checked on
# a sample of catalogued records against their known platform before use.
EVIDENCE = _c(
    [
        ("Visium HD", r"square_\d{3}um|binned_outputs|visium[ _-]?hd|tissue_positions\.parquet"),
        ("Xenium", r"xenium|transcripts\.parquet|cells\.parquet|morphology_focus|cell_feature_matrix\.(h5|zarr)"),
        ("CosMx", r"cosmx|atomx|exprmat_file|_tx_file|fov_positions|spatial molecular imager"),
        ("MERFISH / MERSCOPE", r"merscope|merfish|vizgen|cell_by_gene|detected_transcripts"),
        ("GeoMx DSP", r"geomx|\.dcc\b|\.pkc\b|labworksheet|digital spatial profil"),
        ("Stereo-seq", r"stereo-?seq|stomics|\.gef\b|\.gem(\.gz)?\b"),
        ("Slide-seq", r"slide-?seq|bead_?locations|\bpuck"),
        ("Slide-tags", r"slide-?tags"),
        ("DBiT-seq", r"\bdbit"),
        ("Curio Seeker", r"curio"),
        ("Molecular Cartography", r"molecular cartography|resolve bioscience"),
        ("Tomo-seq", r"tomo-?seq"),
        ("ST arrays (legacy)", r"_stdata\.tsv|st_pipeline|stpipeline"),
        ("Visium", r"visium|space ?ranger|tissue_positions|scalefactors_json|tissue_(hires|lowres)_image|[-_]spatial\.tar"),
    ]
)
SOFT_FILE = re.compile(r"^!(?:Series|Sample)_supplementary_file\w*\s*=\s*(\S+)")
SOFT_TEXT = re.compile(
    r"^!(?:Series_overall_design|Sample_(?:title|description|data_processing|extract_protocol\w*))\s*="
)


def evidence_hits(texts):
    return [name for name, rx in EVIDENCE if any(rx.search(t) for t in texts)]


def decide_platform(file_hits, text_hits):
    """One platform from the evidence, or None when it is absent or conflicting."""
    seen = set(file_hits) | set(text_hits)
    if "Visium HD" in seen and seen <= {"Visium", "Visium HD"}:
        return "Visium HD"  # HD deposits carry Visium-style files too
    if len(file_hits) == 1:
        return file_hits[0]
    if not file_hits and len(text_hits) == 1:
        return text_hits[0]
    common = set(file_hits) & set(text_hits)
    return common.pop() if len(common) == 1 else None


def _soft_parts(txt):
    files, text = [], []
    for line in (txt or "").splitlines():
        m = SOFT_FILE.match(line)
        if m:
            name = m.group(1).rsplit("/", 1)[-1].lower()
            if name and name != "none":
                files.append(name)
        elif SOFT_TEXT.match(line):
            text.append(line.split("=", 1)[-1].lower())
    return files, text


def geo_evidence(acc):
    """File names and processing text of a GEO series and its samples."""
    def soft(targ):
        txt = _request_text(f"{GEO_ACC}?acc={acc}&targ={targ}&form=text&view=brief")
        time.sleep(ENRICH_PAUSE)
        return txt

    txt = soft("self")
    if txt is None:
        return None
    files, text = _soft_parts(txt)
    if not decide_platform(evidence_hits(files), evidence_hits(text)):
        samples = soft("gsm")
        if samples is None:
            return None
        more_files, more_text = _soft_parts(samples)
        files, text = files + more_files, text + more_text
    return {"v": ENRICH_VERSION, "f": evidence_hits(files), "t": evidence_hits(text)}


def _walk_sections(node, visit, depth=0):
    if depth > 6 or not node:
        return
    if isinstance(node, list):
        for x in node:
            _walk_sections(x, visit, depth + 1)
        return
    if isinstance(node, dict):
        visit(node)
        _walk_sections(node.get("subsections"), visit, depth + 1)


def ae_evidence(acc):
    """BioStudies study metadata (study type, organism, sample count, protocols) and file names."""
    study = _request(BIOSTUDIES_API + "studies/" + acc)
    time.sleep(ENRICH_PAUSE)
    if not study:
        return None
    out = {"v": ENRICH_VERSION, "st": [], "org": [], "n": 0}
    protocols = []

    def visit(sec):
        kind = sec.get("type", "")
        for a in sec.get("attributes") or []:
            name, value = a.get("name", ""), str(a.get("value", ""))
            if kind == "Study" and name == "Study type":
                out["st"].append(value)
            elif kind == "Study" and name == "Organism":
                out["org"].append(value)
            elif kind == "Samples" and name == "Sample count":
                try:
                    out["n"] = int(value)
                except ValueError:
                    pass
            elif kind == "Protocols" and name in ("Description", "Hardware"):
                protocols.append(value.lower())

    _walk_sections(study.get("section"), visit)
    listing = _request(BIOSTUDIES_API + "files/" + acc + "?start=0&length=500") or {}
    time.sleep(ENRICH_PAUSE)
    names = [str(d.get("path", "")).rsplit("/", 1)[-1].lower() for d in listing.get("data") or []]
    out["f"], out["t"] = evidence_hits(names), evidence_hits(protocols)
    return out


def needs_evidence(row):
    if row["db"] == "ArrayExpress":
        return True
    if title_assay(row["title"]) not in (None, "spatial"):
        return False  # a tagged single-cell/bulk subseries: its tag already decides
    return row["modality"] in ("spatial", "paired") and row.get("platform_source") != "title"


def apply_evidence(row, ev):
    """Use cached evidence on a record. False = the source says it holds no transcriptome."""
    if row["db"] == "ArrayExpress":
        if ev.get("org"):
            row["organism"] = organism("; ".join(ev["org"]))
        if ev.get("n"):
            row["samples"] = ev["n"]
        types = " ; ".join(ev.get("st") or []).lower()
        if types:
            spatial = "spatial" in types
            single = bool(re.search(r"single[- ]?(cell|nucle)", types))
            if re.search(r"atac|chip-?seq|methyl|bisulfite", types) and not spatial \
                    and not re.search(r"rna|transcript", types):
                return False
            # the curated study type says what this one study holds
            if spatial and single:
                row["modality"] = "paired"
            elif spatial:
                row["modality"] = "spatial"
                if row["platform"] in SC_PLATFORM_NAMES:
                    row["platform"], row["platform_class"], row["platform_source"] = UNSPECIFIED_SPATIAL, "unspecified", ""
            elif single:
                row["modality"] = "singlecell"
                if row["platform"] not in SC_PLATFORM_NAMES:
                    row["platform"] = row.get("sc_platform") or "scRNA-seq (unspecified)"
                    row["platform_class"], row["platform_source"] = "sequencing", ""
    if row["modality"] in ("spatial", "paired") and row.get("platform_source") != "title":
        found = decide_platform(ev.get("f") or [], ev.get("t") or [])
        if found:
            row["platform"], row["platform_class"] = found, platform_class_of(found)
            row["platform_source"] = "files" if found in (ev.get("f") or []) else "protocol"
    return True


def load_enrichment():
    try:
        doc = json.loads(ENRICH_PATH.read_text())
        if doc.get("version") == ENRICH_VERSION:
            return doc.get("records", {})
    except Exception:  # noqa: BLE001 - no cache yet, or an unreadable one
        pass
    return {}


def save_enrichment(cache):
    ENRICH_PATH.write_text(
        json.dumps({"version": ENRICH_VERSION, "records": dict(sorted(cache.items()))},
                   separators=(",", ":"))
    )


def enrich(rows, fetch=True):
    """Fetch missing evidence (within the time budget) and apply all cached evidence."""
    cache = load_enrichment()
    todo = [r for r in rows if needs_evidence(r) and r["accession"] not in cache]
    # Most useful first, in case the budget runs out: ArrayExpress (few, and its
    # study type, organism and sample count fix whole records), then spatial
    # records with no platform at all, then those labelled from the abstract.
    todo.sort(key=lambda r: (r["db"] != "ArrayExpress", r["platform"] != UNSPECIFIED_SPATIAL))
    if fetch and todo:
        deadline = time.time() + ENRICH_BUDGET_S
        print(f"enrichment: fetching {len(todo)} records "
              f"({sum(1 for r in todo if r['db'] == 'GEO')} GEO)", flush=True)

        def work(r):
            if time.time() > deadline:
                return r["accession"], None
            fn = ae_evidence if r["db"] == "ArrayExpress" else geo_evidence
            try:
                return r["accession"], fn(r["accession"])
            except Exception as exc:  # noqa: BLE001 - one bad record must not stop the run
                print(f"  ! enrichment failed for {r['accession']}: {exc}", file=sys.stderr)
                return r["accession"], None

        with ThreadPoolExecutor(ENRICH_WORKERS) as pool:
            for i, (acc, ev) in enumerate(pool.map(work, todo), 1):
                if ev:
                    cache[acc] = ev
                if i % 200 == 0:
                    print(f"  enrichment {i}/{len(todo)}", flush=True)
        missing = sum(1 for r in todo if r["accession"] not in cache)
        if missing:
            print(f"  {missing} records left for the next run (budget or network)", flush=True)
        save_enrichment(cache)

    kept, dropped, changed = [], 0, 0
    for r in rows:
        ev = cache.get(r["accession"]) if needs_evidence(r) else None
        if ev:
            before = (r["platform"], r["modality"])
            if apply_evidence(r, ev) is False:
                dropped += 1
                continue
            changed += before != (r["platform"], r["modality"])
        kept.append(r)
    print(f"enrichment: {changed} records relabelled, {dropped} dropped", flush=True)
    return kept


def geo_url(acc):
    return f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"


def ae_url(acc):
    return f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{acc}"


def build():
    rows = []
    seen = set()

    print("GEO: spatial search", flush=True)
    spatial_ids = esearch(SPATIAL_QUERY)
    print(f"  {len(spatial_ids)} ids", flush=True)
    print("GEO: bare-'ST' phrase search", flush=True)
    st_ids = esearch(ST_QUERY)
    print(f"  {len(st_ids)} ids", flush=True)
    spatial_ids = list(dict.fromkeys(spatial_ids + st_ids))
    spatial = esummary(spatial_ids)

    for r in spatial:
        c = classify(r["title"], r["summary"], r["gdstype"])
        if not c["platform"]:
            continue
        if r["acc"] in seen:
            continue
        seen.add(r["acc"])
        rows.append(
            {
                "accession": r["acc"],
                "db": "GEO",
                "title": squash(r["title"], 190),
                "description": squash(r["summary"], 320),
                "_k": group_key(r["summary" if "summary" in r else "content"]),
                "modality": "paired" if c["has_sc"] else "spatial",
                "platform": c["platform"],
                "platform_class": c["platform_class"],
                "platform_source": c["platform_source"],
                "sc_platform": c["sc_platform"],
                "disease": c["disease"],
                "disease_group": c["disease_group"],
                "subtypes": c["subtypes"],
                "study_type": c["study_type"],
                "tissue": c["tissue"],
                "organism": organism(r["taxon"]),
                "samples": r["n"],
                "cells": c["cells"],
                "year": (r["pdat"] or "")[:4],
                "date": (r["pdat"] or "").replace("/", "-"),
                "pmid": str(r["pmid"] or ""),
                "access": c["access"],
                "url": geo_url(r["acc"]),
                "ftp": r.get("ftp", ""),
            }
        )

    print("GEO: single-cell cancer search", flush=True)
    sc_ids = [i for i in esearch(SC_CANCER_QUERY) if i not in set(spatial_ids)]
    print(f"  {len(sc_ids)} new ids", flush=True)
    sc = esummary(sc_ids)

    for r in sc:
        c = classify(r["title"], r["summary"], r["gdstype"])
        if not c["has_sc"] or c["disease_group"] != "Cancer":
            continue
        if r["acc"] in seen:
            continue
        seen.add(r["acc"])
        rows.append(
            {
                "accession": r["acc"],
                "db": "GEO",
                "title": squash(r["title"], 190),
                "description": squash(r["summary"], 320),
                "_k": group_key(r["summary" if "summary" in r else "content"]),
                # A record reached through the cancer/scRNA query can still be
                # a spatial series (it says "ST", not "spatial transcriptomics").
                # Let the classification decide, not the query it arrived on.
                "modality": (
                    ("paired" if c["has_sc"] else "spatial")
                    if c["platform"]
                    else "singlecell"
                ),
                "platform": c["platform"] or c["sc_platform"] or "scRNA-seq (unspecified)",
                "platform_class": c["platform_class"] or "sequencing",
                "platform_source": c["platform_source"],
                "sc_platform": c["sc_platform"],
                "disease": c["disease"],
                "disease_group": c["disease_group"],
                "subtypes": c["subtypes"],
                "study_type": c["study_type"],
                "tissue": c["tissue"],
                "organism": organism(r["taxon"]),
                "samples": r["n"],
                "cells": c["cells"],
                "year": (r["pdat"] or "")[:4],
                "date": (r["pdat"] or "").replace("/", "-"),
                "pmid": str(r["pmid"] or ""),
                "access": c["access"],
                "url": geo_url(r["acc"]),
                "ftp": r.get("ftp", ""),
            }
        )

    print("ArrayExpress", flush=True)
    ae = {r["acc"]: r for r in ae_search(AE_SPATIAL_QUERY)}
    for r in ae_search(AE_SC_QUERY):
        ae.setdefault(r["acc"], r)

    for r in ae.values():
        # E-GEOD-* records are ArrayExpress mirrors of GEO series; skip duplicates.
        if r["acc"].upper().startswith("E-GEOD"):
            continue
        text = f"{r['title']} {r['content']}"
        if not re.search(r"Homo sapiens|Mus musculus", text, re.I):
            continue
        c = classify(r["title"], r["content"])
        if not c["platform"] and not c["has_sc"]:
            continue
        if not c["platform"] and c["disease_group"] != "Cancer":
            continue
        if r["acc"] in seen:
            continue
        seen.add(r["acc"])
        rows.append(
            {
                "accession": r["acc"],
                "db": "ArrayExpress",
                "title": squash(r["title"], 190),
                "description": squash(r["content"], 320),
                "_k": group_key(r["summary" if "summary" in r else "content"]),
                "modality": (
                    ("paired" if c["has_sc"] else "spatial")
                    if c["platform"]
                    else "singlecell"
                ),
                "platform": c["platform"] or c["sc_platform"] or "scRNA-seq (unspecified)",
                "platform_class": c["platform_class"] or "sequencing",
                "platform_source": c["platform_source"],
                "sc_platform": c["sc_platform"],
                "disease": c["disease"],
                "disease_group": c["disease_group"],
                "subtypes": c["subtypes"],
                "study_type": c["study_type"],
                "tissue": c["tissue"],
                "organism": organism("", text),
                "samples": 0,
                "cells": c["cells"],
                "year": (r["rdate"] or "")[:4],
                "date": r["rdate"],
                "pmid": "",
                "access": c["access"],
                "url": ae_url(r["acc"]),
                "ftp": "",
            }
        )

    edges = link_siblings(rows)
    rows = enrich(rows)
    rows, dropped, promoted = finalize(rows, edges)
    print(f"title assay tags: dropped {dropped} non-spatial, non-single-cell subseries", flush=True)
    print(f"\nsubseries pairing: {promoted} records promoted to paired via a sibling series", flush=True)

    rows.sort(key=lambda x: (x["date"] or ""), reverse=True)
    return rows


def summarise(rows):
    def tally(key):
        out = {}
        for r in rows:
            out[r[key] or "Unspecified"] = out.get(r[key] or "Unspecified", 0) + 1
        return dict(sorted(out.items(), key=lambda kv: -kv[1]))

    return {
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "total": len(rows),
        "paired_via_siblings": sum(1 for r in rows if r.get("pairing") == "series"),
        "records_with_siblings": sum(1 for r in rows if r.get("siblings")),
        "by_modality": tally("modality"),
        "by_platform_class": tally("platform_class"),
        "by_platform_source": tally("platform_source"),
        "by_disease_group": tally("disease_group"),
        "by_study_type": tally("study_type"),
        "by_database": tally("db"),
        "by_year": tally("year"),
        "sources": [
            "NCBI GEO (E-utilities esearch/esummary, db=gds)",
            "EBI ArrayExpress (BioStudies search API)",
            "GEO series/sample records and BioStudies study files (platform evidence)",
        ],
    }


def reclassify(rows):
    """Apply the newer rules to an existing catalog without re-harvesting.

    The stored description is truncated, so only rules that are safe on the
    title + stored text are re-applied: study type, title-first disease,
    mesothelioma and pan-cancer labels, tissue fallback, the [assay] tag of
    subseries, and study linking. GEO links found from the full abstract by the
    last harvest are kept as edges. The next full harvest recomputes everything.
    """
    edges = []
    unlinkable = {r["accession"] for r in rows if not linkable(r)}
    for r in rows:
        if r["accession"] not in unlinkable:
            edges += [(r["accession"], s) for s in r.get("siblings", []) if s not in unlinkable]
        # undo the previous sibling promotion; it is recomputed below
        if r.get("pairing") == "series":
            r["modality"] = "singlecell" if r["platform"] in SC_PLATFORM_NAMES else "spatial"

    for r in rows:
        r["study_type"] = study_type(r["title"], r["description"])
        text = f"{r['title']} || {r['description']}"
        if r["disease"] == "Lung cancer (other)" and re.search(r"mesothelioma", text, re.I):
            r["disease"] = "Mesothelioma"
        title_cancer = first_match(CANCER, DISEASE_NOISE.sub(" ", r["title"]))
        if r["disease_group"] == "Cancer" and title_overrides(title_cancer, r["disease"]):
            r["disease"] = title_cancer
        if r["disease_group"] == "Cancer" and PAN_CANCER_TITLE.search(r["title"]):
            r["disease"] = GENERIC_CANCER
        r["tissue"] = tissue_for(r["tissue"], r["disease"])
        r["subtypes"] = compatible_subtypes(r["subtypes"], r["disease"])
        r["date"] = (r["date"] or "").replace("/", "-")
        if "platform_source" not in r:
            r["platform_source"] = platform_source(r["title"], r["platform"])
        # a platform named in the record's own title wins over the abstract's
        if r["platform"] not in SC_PLATFORM_NAMES:
            named, cls = detect_platform(r["title"], "")
            if named and named != UNSPECIFIED_SPATIAL:
                r["platform"], r["platform_class"], r["platform_source"] = named, cls, "title"
        if r["platform"] == UNSPECIFIED_SPATIAL:
            r["platform_class"], r["platform_source"] = "unspecified", ""
    rows = enrich(rows, fetch="--fetch" in sys.argv)
    edges += geo_title_links(rows)
    edges += ae_title_links(rows)
    kept, dropped, promoted = finalize(rows, edges)
    print(f"reclassify: dropped {dropped} non-spatial/non-single-cell subseries, "
          f"{promoted} records paired via siblings", flush=True)
    return kept


def reclassify_main():
    existing = DATA / "datasets.json"
    doc = json.loads(existing.read_text())
    rows = reclassify(doc["datasets"])
    if len(rows) < 200:
        print("reclassify left too few records - refusing to write", file=sys.stderr)
        return 1
    meta = summarise(rows)
    meta["new_since_last_run"] = doc.get("meta", {}).get("new_since_last_run", 0)
    meta["generated"] = doc.get("meta", {}).get("generated", meta["generated"])
    meta["reclassified"] = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    existing.write_text(
        json.dumps({"meta": meta, "datasets": rows}, ensure_ascii=False, separators=(",", ":"))
    )
    (DATA / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2))
    print(f"Reclassified {len(rows)} datasets.")
    return 0


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    if "--reclassify" in sys.argv:
        return reclassify_main()
    rows = build()

    if len(rows) < 200:
        print(
            f"Harvest returned only {len(rows)} records - refusing to overwrite "
            "the existing catalog with a likely-truncated result.",
            file=sys.stderr,
        )
        return 1

    previous = set()
    existing = DATA / "datasets.json"
    if existing.exists():
        try:
            previous = {
                r["accession"] for r in json.loads(existing.read_text())["datasets"]
            }
        except Exception:  # noqa: BLE001
            previous = set()

    new_accessions = [r["accession"] for r in rows if r["accession"] not in previous]

    meta = summarise(rows)
    meta["new_since_last_run"] = len(new_accessions) if previous else 0

    existing.write_text(
        json.dumps({"meta": meta, "datasets": rows}, ensure_ascii=False, separators=(",", ":"))
    )
    (DATA / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2))

    print(f"\nWrote {len(rows)} datasets ({meta['new_since_last_run']} new).")
    print(json.dumps(meta["by_modality"], indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
