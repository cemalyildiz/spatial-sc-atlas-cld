/* Spatial & Single-Cell Atlas — browse, filter, compare, summarise. */
(function () {
  'use strict';

  // ---------------------------------------------------------------- state --
  var DATA = [];
  var META = {};
  var PAGE = 60;
  var shown = PAGE;
  var picked = new Map();      // accession -> record
  var lastFiltered = [];
  var byAccession = Object.create(null);

  var F = {
    q: '',
    sort: 'relevance',
    modality: new Set(),
    platform_class: new Set(),
    platform: new Set(),
    disease_group: new Set(),
    disease: new Set(),
    subtypes: new Set(),
    study_type: new Set(),
    tissue: new Set(),
    organism: new Set(),
    db: new Set(),
    access: new Set(),
    yearMin: null,
    yearMax: null,
    cellsMin: null
  };

  var MODALITY_LABEL = {
    paired: 'Paired (spatial + single-cell)',
    spatial: 'Spatial only',
    singlecell: 'Single-cell only'
  };
  var MODALITY_SHORT = { paired: 'Paired', spatial: 'Spatial', singlecell: 'Single-cell' };
  var CLASS_LABEL = { imaging: 'Imaging-based', sequencing: 'Sequencing-based', unspecified: 'Not stated' };
  // where the platform label of a spatial record came from, strongest first
  var SOURCE_LABEL = {
    title: "Named in the record's own title",
    files: 'Read from the deposited file names',
    protocol: 'Read from the processing / protocol text',
    abstract: 'Named in the study abstract'
  };

  var GROUPS = [
    { key: 'modality',       title: 'Modality',           open: true,  label: function (v) { return MODALITY_LABEL[v] || v; } },
    { key: 'platform_class', title: 'Production method',  open: true,  label: function (v) { return CLASS_LABEL[v] || v; } },
    { key: 'disease_group',  title: 'Disease group',      open: true },
    { key: 'disease',        title: 'Disease / cancer type', open: true },
    { key: 'subtypes',       title: 'Subtype', multi: true },
    { key: 'study_type',     title: 'Study type' },
    { key: 'platform',       title: 'Platform' },
    { key: 'tissue',         title: 'Tissue / organ' },
    { key: 'organism',       title: 'Organism' },
    { key: 'db',             title: 'Source database' },
    { key: 'access',         title: 'Access' }
  ];

  // --------------------------------------------------------------- helpers --
  var $ = function (s, r) { return (r || document).querySelector(s); };
  var $$ = function (s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); };

  function el(tag, cls, text) {
    var n = document.createElement(tag);
    if (cls) n.className = cls;
    if (text != null) n.textContent = text;
    return n;
  }

  function nf(n) {
    if (n == null || isNaN(n)) return '—';
    return Number(n).toLocaleString('en-US');
  }

  function compact(n) {
    if (!n) return '—';
    if (n >= 1e6) return (n / 1e6).toFixed(n >= 1e7 ? 0 : 1).replace(/\.0$/, '') + 'M';
    if (n >= 1e3) return (n / 1e3).toFixed(n >= 1e4 ? 0 : 1).replace(/\.0$/, '') + 'k';
    return String(n);
  }

  var toastTimer;
  function toast(msg) {
    var t = $('#toast');
    t.textContent = msg;
    t.hidden = false;
    clearTimeout(toastTimer);
    toastTimer = setTimeout(function () { t.hidden = true; }, 2200);
  }

  function download(name, text, mime) {
    var blob = new Blob([text], { type: mime || 'text/plain;charset=utf-8' });
    var url = URL.createObjectURL(blob);
    var a = document.createElement('a');
    a.href = url; a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
  }

  function copyText(text, okMsg) {
    function fallback() {
      var ta = document.createElement('textarea');
      ta.value = text;
      ta.style.position = 'fixed';
      ta.style.opacity = '0';
      document.body.appendChild(ta);
      ta.select();
      try { document.execCommand('copy'); toast(okMsg); }
      catch (e) { toast('Could not copy — select the text manually'); }
      ta.remove();
    }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(function () { toast(okMsg); }, fallback);
    } else { fallback(); }
  }

  // ------------------------------------------------------------ filtering --
  function haystack(r) {
    if (r._h) return r._h;
    r._h = (r.accession + ' ' + r.title + ' ' + r.description + ' ' + r.disease + ' ' +
            r.tissue + ' ' + r.platform + ' ' + r.organism + ' ' + (r.sc_platform || '') +
            ' ' + (r.subtypes || []).join(' ') +
            ' ' + (r.siblings || []).join(' ')).toLowerCase();
    return r._h;
  }

  function matchesExcept(r, skipKey) {
    for (var i = 0; i < GROUPS.length; i++) {
      var k = GROUPS[i].key;
      if (k === skipKey) continue;
      var set = F[k];
      if (!set.size) continue;
      var val = r[k];
      if (Array.isArray(val)) {
        // multi-value field: the record matches if it carries any selected value
        var hit = false;
        for (var v = 0; v < val.length; v++) { if (set.has(val[v])) { hit = true; break; } }
        if (!hit) return false;
      } else if (!set.has(val == null ? '' : String(val))) {
        return false;
      }
    }
    if (skipKey !== '_year') {
      var y = parseInt(r.year, 10);
      if (F.yearMin != null && (!y || y < F.yearMin)) return false;
      if (F.yearMax != null && (!y || y > F.yearMax)) return false;
    }
    if (F.cellsMin != null && (r.cells || 0) < F.cellsMin) return false;
    if (F.q) {
      var hs = haystack(r);
      var terms = F.q.toLowerCase().split(/\s+/).filter(Boolean);
      for (var j = 0; j < terms.length; j++) if (hs.indexOf(terms[j]) === -1) return false;
    }
    return true;
  }

  function filtered() {
    return DATA.filter(function (r) { return matchesExcept(r, null); });
  }

  function facetCounts(key) {
    var counts = Object.create(null);
    for (var i = 0; i < DATA.length; i++) {
      var r = DATA[i];
      if (!matchesExcept(r, key)) continue;
      var raw = r[key];
      if (Array.isArray(raw)) {
        for (var a = 0; a < raw.length; a++) {
          if (raw[a]) counts[raw[a]] = (counts[raw[a]] || 0) + 1;
        }
        continue;
      }
      var v = raw == null ? '' : String(raw);
      if (!v) continue;
      counts[v] = (counts[v] || 0) + 1;
    }
    return counts;
  }

  // GEO wrote dates as YYYY/MM/DD and ArrayExpress as YYYY-MM-DD; compare one form
  // so an older catalog still sorts correctly.
  function iso(d) { return (d || '').replace(/\//g, '-'); }

  var SORTERS = {
    relevance: function (a, b) {
      var rank = { paired: 0, spatial: 1, singlecell: 2 };
      var d = rank[a.modality] - rank[b.modality];
      if (d) return d;
      return iso(b.date) < iso(a.date) ? -1 : iso(b.date) > iso(a.date) ? 1 : 0;
    },
    newest: function (a, b) { return iso(b.date) < iso(a.date) ? -1 : iso(b.date) > iso(a.date) ? 1 : 0; },
    oldest: function (a, b) { return iso(a.date) < iso(b.date) ? -1 : iso(a.date) > iso(b.date) ? 1 : 0; },
    cells: function (a, b) { return (b.cells || 0) - (a.cells || 0); },
    samples: function (a, b) { return (b.samples || 0) - (a.samples || 0); },
    // numeric-aware, so GSE99795 comes before GSE100000
    accession: function (a, b) { return a.accession.localeCompare(b.accession, 'en', { numeric: true }); }
  };

  // --------------------------------------------------------------- render --
  function badge(r) {
    var b = el('span', 'badge badge-' + (r.modality === 'singlecell' ? 'sc' : r.modality));
    b.textContent = MODALITY_SHORT[r.modality] || r.modality;
    if (r.modality === 'paired') {
      b.title = r.pairing === 'series'
        ? 'Paired at the study level: a sibling series carries the other half — see "Other series in this study"'
        : 'Both spatial and single-cell data are referenced in this record';
    }
    return b;
  }

  function card(r) {
    var c = el('div', 'card' + (picked.has(r.accession) ? ' is-picked' : ''));
    c.setAttribute('data-acc', r.accession);
    c.setAttribute('role', 'button');
    c.tabIndex = 0;

    var pick = el('div', 'card-pick');
    var cb = document.createElement('input');
    cb.type = 'checkbox';
    cb.checked = picked.has(r.accession);
    cb.setAttribute('aria-label', 'Select ' + r.accession + ' for comparison');
    cb.addEventListener('click', function (ev) {
      ev.stopPropagation();
      togglePick(r, c, cb.checked);
    });
    pick.appendChild(cb);

    var main = el('div', 'card-main');

    var top = el('div', 'card-top');
    top.appendChild(el('span', 'acc', r.accession));
    top.appendChild(badge(r));
    if (r.platform) {
      var t = el('span', 'tech tech-' + (r.platform_class || 'sequencing'), r.platform);
      t.title = (CLASS_LABEL[r.platform_class] || '') + ' platform';
      top.appendChild(t);
    }
    if (r.siblings && r.siblings.length) {
      var sib = el('span', 'sib-chip', '+' + r.siblings.length + ' series');
      sib.title = 'Same study, other series: ' + r.siblings.join(', ');
      top.appendChild(sib);
    }
    if (r.study_type === 'Method / tool development') {
      var mc = el('span', 'method-chip', 'Method');
      mc.title = 'Method / tool development study — the data are often a demonstration set';
      top.appendChild(mc);
    }
    (r.subtypes || []).slice(0, 3).forEach(function (s) {
      var chip = el('span', 'sub-chip', s);
      chip.title = 'Subtype detected in the record text';
      top.appendChild(chip);
    });
    main.appendChild(top);

    main.appendChild(el('h3', null, r.title));

    var meta = el('div', 'meta');
    function bit(label, value) {
      if (!value || value === 'Unspecified' || value === 'Other') return;
      var s = el('span');
      s.appendChild(el('b', null, value));
      s.appendChild(document.createTextNode(' ' + label));
      meta.appendChild(s);
    }
    bit('', r.disease !== 'Unspecified' ? r.disease : '');
    bit('', r.tissue !== 'Other' ? r.tissue : '');
    bit('', r.organism);
    if (r.samples) bit(r.samples === 1 ? 'sample' : 'samples', nf(r.samples));
    if (r.cells) bit('cells', compact(r.cells));
    if (r.year) bit('', r.year);
    meta.appendChild(el('span', null, r.db));
    main.appendChild(meta);

    c.appendChild(pick);
    c.appendChild(main);

    c.addEventListener('click', function () { openDetail(r); });
    c.addEventListener('keydown', function (ev) {
      if (ev.key === 'Enter' || ev.key === ' ') { ev.preventDefault(); openDetail(r); }
    });
    return c;
  }

  function renderResults() {
    var rows = filtered();
    rows.sort(SORTERS[F.sort] || SORTERS.relevance);
    lastFiltered = rows;

    var host = $('#results');
    host.textContent = '';

    var total = rows.length;
    var pairedN = 0;
    for (var i = 0; i < total; i++) if (rows[i].modality === 'paired') pairedN++;

    $('#resultCount').innerHTML = total
      ? '<strong>' + nf(total) + '</strong> dataset' + (total === 1 ? '' : 's') +
        ' · <strong>' + nf(pairedN) + '</strong> paired'
      : 'No datasets match these filters';

    if (!total) {
      var e = el('div', 'empty');
      e.appendChild(el('h3', null, DATA.length ? 'Nothing matches' : 'Catalog not built yet'));
      e.appendChild(el('p', null, DATA.length
        ? 'Try clearing a filter or two — the counts beside each option show what is still reachable.'
        : 'Run the "Refresh catalog" workflow in the repository to harvest datasets from GEO and ArrayExpress.'));
      host.appendChild(e);
      renderChips();
      return;
    }

    var frag = document.createDocumentFragment();
    var limit = Math.min(shown, total);
    for (var j = 0; j < limit; j++) frag.appendChild(card(rows[j]));
    host.appendChild(frag);

    if (limit < total) {
      var more = el('button', 'ghost-btn loadmore', 'Show ' + nf(Math.min(PAGE, total - limit)) + ' more');
      more.addEventListener('click', function () { shown += PAGE; renderResults(); });
      host.appendChild(more);
    }

    renderChips();
  }

  function renderChips() {
    var host = $('#activeChips');
    host.textContent = '';
    var n = 0;

    GROUPS.forEach(function (g) {
      F[g.key].forEach(function (v) {
        n++;
        host.appendChild(chip((g.label ? g.label(v) : v), function () {
          F[g.key].delete(v); onFilterChange();
        }));
      });
    });
    if (F.yearMin != null || F.yearMax != null) {
      n++;
      host.appendChild(chip('Year ' + (F.yearMin || '…') + '–' + (F.yearMax || '…'), function () {
        F.yearMin = F.yearMax = null; onFilterChange();
      }));
    }
    if (F.cellsMin != null) {
      n++;
      host.appendChild(chip('≥ ' + compact(F.cellsMin) + ' cells', function () {
        F.cellsMin = null; onFilterChange();
      }));
    }
    if (F.q) {
      n++;
      host.appendChild(chip('“' + F.q + '”', function () {
        F.q = ''; $('#search').value = ''; onFilterChange();
      }));
    }

    var pill = $('#activeFilterCount');
    pill.textContent = String(n);
    pill.hidden = n === 0;
  }

  function chip(label, onRemove) {
    var c = el('span', 'chip');
    c.appendChild(document.createTextNode(label));
    var b = el('button', null, '×');
    b.setAttribute('aria-label', 'Remove filter ' + label);
    b.addEventListener('click', onRemove);
    c.appendChild(b);
    return c;
  }

  // -------------------------------------------------------------- filters --
  function buildFilters() {
    var host = $('#filterGroups');
    host.textContent = '';

    GROUPS.forEach(function (g) {
      var d = el('details', 'fgroup');
      d.open = !!g.open || F[g.key].size > 0;
      var s = el('summary');
      s.appendChild(document.createTextNode(g.title));
      var n = el('span', 'n');
      s.appendChild(n);
      d.appendChild(s);

      var box = el('div', 'fopts');
      box.setAttribute('data-key', g.key);
      d.appendChild(box);
      host.appendChild(d);
      g._countEl = n;
      g._box = box;
    });

    // year range
    var dy = el('details', 'fgroup');
    dy.open = F.yearMin != null || F.yearMax != null;
    var sy = el('summary');
    sy.appendChild(document.createTextNode('Release year'));
    dy.appendChild(sy);
    var wrap = el('div', 'rangewrap');
    var from = document.createElement('input');
    from.type = 'number'; from.placeholder = 'from'; from.min = '2010'; from.max = '2100';
    from.value = F.yearMin || '';
    from.setAttribute('aria-label', 'Earliest release year');
    var to = document.createElement('input');
    to.type = 'number'; to.placeholder = 'to'; to.min = '2010'; to.max = '2100';
    to.value = F.yearMax || '';
    to.setAttribute('aria-label', 'Latest release year');
    function applyYears() {
      F.yearMin = from.value ? parseInt(from.value, 10) : null;
      F.yearMax = to.value ? parseInt(to.value, 10) : null;
      onFilterChange(true);
    }
    from.addEventListener('change', applyYears);
    to.addEventListener('change', applyYears);
    wrap.appendChild(from);
    wrap.appendChild(el('span', null, '–'));
    wrap.appendChild(to);
    dy.appendChild(wrap);
    host.appendChild(dy);

    // minimum cells
    var dc = el('details', 'fgroup');
    dc.open = F.cellsMin != null;
    var sc = el('summary');
    sc.appendChild(document.createTextNode('Minimum cells'));
    dc.appendChild(sc);
    var cwrap = el('div', 'fopts');
    [0, 10000, 50000, 100000, 500000].forEach(function (v) {
      var lab = el('label', 'fopt');
      var input = document.createElement('input');
      input.type = 'radio'; input.name = 'cellsMin';
      input.checked = (F.cellsMin || 0) === v;
      input.addEventListener('change', function () {
        F.cellsMin = v || null;
        onFilterChange(true);
      });
      lab.appendChild(input);
      lab.appendChild(el('span', 'lbl', v ? '≥ ' + compact(v) : 'Any'));
      cwrap.appendChild(lab);
    });
    dc.appendChild(cwrap);
    // a cell count is only known where the record states one, so say how many do
    var stated = DATA.filter(function (r) { return r.cells > 0; }).length;
    dc.appendChild(el('p', 'fnote', 'Only ' + nf(stated) + ' of ' + nf(DATA.length) +
      ' records state a cell count; any minimum hides the rest.'));
    host.appendChild(dc);

    refreshFacets();
  }

  function refreshFacets() {
    GROUPS.forEach(function (g) {
      var counts = facetCounts(g.key);
      var values = Object.keys(counts);

      // keep active-but-now-zero values visible so they can be switched off
      F[g.key].forEach(function (v) { if (values.indexOf(v) === -1) values.push(v); });

      values.sort(function (a, b) {
        var d = (counts[b] || 0) - (counts[a] || 0);
        return d || a.localeCompare(b);
      });

      g._countEl.textContent = F[g.key].size ? String(F[g.key].size) : '';
      var box = g._box;
      box.textContent = '';

      values.forEach(function (v) {
        var on = F[g.key].has(v);
        var c = counts[v] || 0;
        var lab = el('label', 'fopt' + (c ? '' : ' is-empty') + (on ? ' is-on' : ''));
        var input = document.createElement('input');
        input.type = 'checkbox';
        input.checked = on;
        input.addEventListener('change', function () {
          if (input.checked) F[g.key].add(v); else F[g.key].delete(v);
          onFilterChange();
        });
        lab.appendChild(input);
        var text = g.label ? g.label(v) : v;
        var span = el('span', 'lbl', text);
        span.title = text;
        lab.appendChild(span);
        lab.appendChild(el('span', 'cnt', nf(c)));
        box.appendChild(lab);
      });

      if (!values.length) box.appendChild(el('div', 'fopt is-empty', '—'));
    });
  }

  function onFilterChange(skipFacetRebuild) {
    shown = PAGE;
    if (!skipFacetRebuild) refreshFacets(); else refreshFacets();
    renderResults();
    renderStats();
    writeURL();
  }

  function resetFilters() {
    GROUPS.forEach(function (g) { F[g.key].clear(); });
    F.yearMin = F.yearMax = F.cellsMin = null;
    F.q = '';
    $('#search').value = '';
    buildFilters();
    onFilterChange();
  }

  // ------------------------------------------------------------ URL state --
  var URL_KEYS = { modality: 'm', platform_class: 'pc', platform: 'p', disease_group: 'dg',
                   disease: 'd', subtypes: 'st', study_type: 'sy', tissue: 't', organism: 'o', db: 'db', access: 'ac' };

  function writeURL() {
    var p = new URLSearchParams();
    Object.keys(URL_KEYS).forEach(function (k) {
      if (F[k].size) p.set(URL_KEYS[k], Array.from(F[k]).join('|'));
    });
    if (F.q) p.set('q', F.q);
    if (F.sort !== 'relevance') p.set('sort', F.sort);
    if (F.yearMin != null) p.set('y1', F.yearMin);
    if (F.yearMax != null) p.set('y2', F.yearMax);
    if (F.cellsMin != null) p.set('c', F.cellsMin);
    var view = $('.tab.is-active').getAttribute('data-view');
    if (view !== 'browse') p.set('v', view);
    var qs = p.toString();
    history.replaceState(null, '', qs ? '?' + qs : location.pathname);
  }

  function readURL() {
    var p = new URLSearchParams(location.search);
    Object.keys(URL_KEYS).forEach(function (k) {
      var v = p.get(URL_KEYS[k]);
      if (v) v.split('|').forEach(function (x) { if (x) F[k].add(x); });
    });
    F.q = p.get('q') || '';
    F.sort = p.get('sort') || 'relevance';
    F.yearMin = p.get('y1') ? parseInt(p.get('y1'), 10) : null;
    F.yearMax = p.get('y2') ? parseInt(p.get('y2'), 10) : null;
    F.cellsMin = p.get('c') ? parseInt(p.get('c'), 10) : null;
    $('#search').value = F.q;
    $('#sort').value = F.sort;
    var v = p.get('v');
    if (v) switchView(v);
  }

  // --------------------------------------------------------------- detail --
  function openDetail(r) {
    $('#panelTitle').textContent = r.accession;
    var body = $('#panelBody');
    body.textContent = '';

    var top = el('div', 'card-top');
    top.style.marginBottom = '12px';
    top.appendChild(badge(r));
    if (r.platform) top.appendChild(el('span', 'tech tech-' + (r.platform_class || 'sequencing'), r.platform));
    body.appendChild(top);

    var h = el('h3', null, r.title);
    h.style.margin = '0 0 14px';
    h.style.fontSize = '15px';
    h.style.lineHeight = '1.45';
    body.appendChild(h);

    if (r.description) body.appendChild(el('p', 'detail-desc', r.description));

    var dl = el('dl', 'dl');
    function row(k, v) {
      if (!v && v !== 0) return;
      dl.appendChild(el('dt', null, k));
      dl.appendChild(el('dd', null, v));
    }
    row('Modality', MODALITY_LABEL[r.modality] || r.modality);
    row('Platform', r.platform + (r.platform_class ? ' · ' + CLASS_LABEL[r.platform_class] : ''));
    if (r.platform_source && SOURCE_LABEL[r.platform_source]) row('Platform evidence', SOURCE_LABEL[r.platform_source]);
    if (r.sc_platform && r.modality !== 'singlecell') row('Single-cell assay', r.sc_platform);
    row('Disease', r.disease + (r.disease_group !== 'Unspecified' ? ' (' + r.disease_group + ')' : ''));
    if (r.subtypes && r.subtypes.length) row('Subtype', r.subtypes.join(' · '));
    row('Study type', r.study_type);
    if (r.siblings && r.siblings.length) {
      dl.appendChild(el('dt', null, 'Other series'));
      var dd = el('dd');
      r.siblings.forEach(function (acc, i) {
        if (i) dd.appendChild(document.createTextNode(' · '));
        var other = byAccession[acc];
        if (other) {
          var a = el('a', null, acc);
          a.href = '#';
          a.addEventListener('click', function (ev) { ev.preventDefault(); openDetail(other); });
          a.title = other.title;
          dd.appendChild(a);
        } else {
          dd.appendChild(document.createTextNode(acc));
        }
      });
      dd.appendChild(el('div', 'dd-note',
        r.pairing === 'series'
          ? 'Same study, split across ' + r.db + ' series — this record is counted as paired because of them.'
          : 'Same study, split across ' + r.db + ' series.'));
      dl.appendChild(dd);
    }
    row('Tissue / organ', r.tissue);
    row('Organism', r.organism);
    row('Samples', r.samples ? nf(r.samples) : null);
    row('Cells reported', r.cells ? nf(r.cells) : null);
    row('Released', iso(r.date) || r.year);
    row('Access', r.access);
    row('Source', r.db);
    body.appendChild(dl);

    var links = el('div', 'detail-links');
    var a1 = el('a', 'is-primary', 'Open in ' + r.db);
    a1.href = r.url; a1.target = '_blank'; a1.rel = 'noopener';
    links.appendChild(a1);

    if (r.db === 'GEO') {
      // GEO reports an ftp:// link; browsers dropped that scheme, so use the
      // HTTPS mirror of the same path.
      var base = (r.ftp || ('ftp://ftp.ncbi.nlm.nih.gov/geo/series/' +
        r.accession.slice(0, -3) + 'nnn/' + r.accession + '/'))
        .replace(/^ftp:\/\//, 'https://');
      if (base.slice(-1) !== '/') base += '/';
      // Link the series root, not suppl/ — a series with no supplementary
      // files has no suppl/ directory and that link would 404.
      var a2 = el('a', null, 'Download (FTP)');
      a2.href = base;
      a2.target = '_blank'; a2.rel = 'noopener';
      a2.title = 'Series directory: suppl/ holds processed matrices and images, matrix/ the series matrix';
      links.appendChild(a2);
    }
    if (r.pmid) {
      var a3 = el('a', null, 'PubMed ' + r.pmid);
      a3.href = 'https://pubmed.ncbi.nlm.nih.gov/' + r.pmid + '/';
      a3.target = '_blank'; a3.rel = 'noopener';
      links.appendChild(a3);
    }
    var a4 = el('a', null, 'Copy accession');
    a4.href = '#';
    a4.addEventListener('click', function (ev) {
      ev.preventDefault();
      copyText(r.accession, r.accession + ' copied');
    });
    links.appendChild(a4);
    body.appendChild(links);

    var note = el('p', null,
      'Modality, platform, disease and tissue are auto-classified from the record text — confirm them on the source page before use.');
    note.style.cssText = 'margin:18px 0 0;font-size:11.5px;color:var(--text-muted);';
    body.appendChild(note);

    $('#overlay').hidden = false;
  }

  // -------------------------------------------------------------- compare --
  function togglePick(r, cardEl, on) {
    if (on) {
      if (picked.size >= 6) {
        toast('Six datasets is the comfortable maximum to compare');
        var cb = cardEl && $('input', cardEl);
        if (cb) cb.checked = false;
        return;
      }
      picked.set(r.accession, r);
    } else {
      picked.delete(r.accession);
    }
    if (cardEl) cardEl.classList.toggle('is-picked', on);
    updateTray();
  }

  function updateTray() {
    var tray = $('#tray');
    tray.hidden = picked.size === 0;
    $('#trayCount').textContent = picked.size + ' selected';
  }

  function openCompare() {
    if (!picked.size) return;
    var rows = Array.from(picked.values());
    $('#panelTitle').textContent = 'Comparing ' + rows.length + ' datasets';
    var body = $('#panelBody');
    body.textContent = '';

    var scroll = el('div', 'cmp-scroll');
    var table = el('table', 'cmp');

    var fields = [
      ['Accession', function (r) { return r.accession; }],
      ['Title', function (r) { return r.title; }],
      ['Modality', function (r) { return MODALITY_LABEL[r.modality] || r.modality; }],
      ['Platform', function (r) { return r.platform; }],
      ['Production method', function (r) { return CLASS_LABEL[r.platform_class] || '—'; }],
      ['Platform evidence', function (r) { return SOURCE_LABEL[r.platform_source] || '—'; }],
      ['Disease', function (r) { return r.disease; }],
      ['Subtype', function (r) { return (r.subtypes || []).join(' · ') || '\u2014'; }],
      ['Study type', function (r) { return r.study_type || '\u2014'; }],
      ['Other series', function (r) { return (r.siblings || []).join(' · ') || '\u2014'; }],
      ['Tissue', function (r) { return r.tissue; }],
      ['Organism', function (r) { return r.organism; }],
      ['Samples', function (r) { return r.samples ? nf(r.samples) : '—'; }],
      ['Cells', function (r) { return r.cells ? nf(r.cells) : '—'; }],
      ['Released', function (r) { return iso(r.date) || r.year || '—'; }],
      ['Source', function (r) { return r.db; }]
    ];

    fields.forEach(function (f) {
      var tr = document.createElement('tr');
      tr.appendChild(el('th', null, f[0]));
      rows.forEach(function (r) {
        var td = el('td');
        if (f[0] === 'Accession') {
          var a = el('a', null, r.accession);
          a.href = r.url; a.target = '_blank'; a.rel = 'noopener';
          td.appendChild(a);
        } else {
          td.textContent = f[1](r);
        }
        tr.appendChild(td);
      });
      table.appendChild(tr);
    });

    scroll.appendChild(table);
    body.appendChild(scroll);

    var acts = el('div', 'detail-links');
    acts.style.marginTop = '18px';
    var copy = el('a', null, 'Copy accessions');
    copy.href = '#';
    copy.addEventListener('click', function (ev) {
      ev.preventDefault();
      copyText(rows.map(function (r) { return r.accession; }).join('\n'), 'Accessions copied');
    });
    acts.appendChild(copy);
    body.appendChild(acts);

    $('#overlay').hidden = false;
  }

  // ---------------------------------------------------------------- stats --
  function tally(rows, key) {
    var m = Object.create(null);
    rows.forEach(function (r) {
      var v = r[key];
      if (!v) return;
      m[v] = (m[v] || 0) + 1;
    });
    return Object.keys(m).map(function (k) { return [k, m[k]]; })
      .sort(function (a, b) { return b[1] - a[1]; });
  }

  function tile(k, v, s) {
    var t = el('div', 'tile');
    t.appendChild(el('div', 'k', k));
    t.appendChild(el('div', 'v', v));
    if (s) t.appendChild(el('div', 's', s));
    return t;
  }

  function hbarChart(title, sub, pairs, total, limit) {
    var c = el('div', 'chartcard');
    c.appendChild(el('h3', null, title));
    c.appendChild(el('p', 'sub', sub));
    var wrap = el('div', 'hbars');
    var top = pairs.slice(0, limit || 12);
    var max = top.length ? top[0][1] : 1;
    top.forEach(function (p) {
      var row = el('div', 'hbar');
      var name = el('div', 'name', p[0]);
      name.title = p[0];
      row.appendChild(name);
      var track = el('div', 'track');
      var fill = el('div', 'fill');
      fill.style.width = Math.max(2, (p[1] / max) * 100) + '%';
      track.appendChild(fill);
      track.title = p[0] + ': ' + nf(p[1]) + ' of ' + nf(total) +
        ' (' + ((p[1] / total) * 100).toFixed(1) + '%)';
      row.appendChild(track);
      row.appendChild(el('div', 'val', nf(p[1])));
      wrap.appendChild(row);
    });
    if (!top.length) wrap.appendChild(el('div', 'sub', 'Nothing in the current selection.'));
    c.appendChild(wrap);
    return c;
  }

  function yearChart(rows) {
    var c = el('div', 'chartcard wide');
    c.appendChild(el('h3', null, 'Datasets released per year'));
    c.appendChild(el('p', 'sub', 'By the release date recorded in the source database. The current year is still filling up.'));

    var byYear = {};
    rows.forEach(function (r) {
      var y = parseInt(r.year, 10);
      if (y >= 2010 && y <= 2100) byYear[y] = (byYear[y] || 0) + 1;
    });
    var years = Object.keys(byYear).map(Number).sort(function (a, b) { return a - b; });

    if (!years.length) {
      c.appendChild(el('p', 'sub', 'Nothing in the current selection.'));
      return c;
    }

    var W = 720, H = 220, padL = 34, padR = 8, padT = 18, padB = 26;
    var innerW = W - padL - padR, innerH = H - padT - padB;
    var max = Math.max.apply(null, years.map(function (y) { return byYear[y]; }));
    var step = innerW / years.length;
    var bw = Math.min(48, step * 0.62);

    var ns = 'http://www.w3.org/2000/svg';
    var svg = document.createElementNS(ns, 'svg');
    svg.setAttribute('class', 'vchart');
    svg.setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    svg.setAttribute('role', 'img');
    svg.setAttribute('aria-label', 'Bar chart of dataset counts per release year');

    [0, 0.5, 1].forEach(function (f) {
      var y = padT + innerH - innerH * f;
      var line = document.createElementNS(ns, 'line');
      line.setAttribute('class', 'grid');
      line.setAttribute('x1', padL); line.setAttribute('x2', W - padR);
      line.setAttribute('y1', y); line.setAttribute('y2', y);
      svg.appendChild(line);
      var lab = document.createElementNS(ns, 'text');
      lab.setAttribute('class', 'axis');
      lab.setAttribute('x', padL - 6); lab.setAttribute('y', y + 3);
      lab.setAttribute('text-anchor', 'end');
      var av = Math.round(max * f);
      lab.textContent = av === 0 ? '0' : compact(av);
      svg.appendChild(lab);
    });

    years.forEach(function (y, i) {
      var v = byYear[y];
      var h = Math.max(2, (v / max) * innerH);
      var x = padL + i * step + (step - bw) / 2;
      var yy = padT + innerH - h;

      var rect = document.createElementNS(ns, 'rect');
      rect.setAttribute('class', 'bar');
      rect.setAttribute('x', x); rect.setAttribute('y', yy);
      rect.setAttribute('width', bw); rect.setAttribute('height', h);
      rect.setAttribute('rx', Math.min(4, bw / 2));
      var t = document.createElementNS(ns, 'title');
      t.textContent = y + ': ' + nf(v) + ' datasets';
      rect.appendChild(t);
      svg.appendChild(rect);

      var val = document.createElementNS(ns, 'text');
      val.setAttribute('class', 'vlabel');
      val.setAttribute('x', x + bw / 2); val.setAttribute('y', yy - 5);
      val.setAttribute('text-anchor', 'middle');
      val.textContent = compact(v);
      svg.appendChild(val);

      var lab = document.createElementNS(ns, 'text');
      lab.setAttribute('class', 'axis');
      lab.setAttribute('x', x + bw / 2); lab.setAttribute('y', H - 8);
      lab.setAttribute('text-anchor', 'middle');
      lab.textContent = years.length > 20 ? String(y).slice(2) : String(y);
      svg.appendChild(lab);
    });

    c.appendChild(svg);
    return c;
  }

  function modalityChart(rows) {
    var c = el('div', 'chartcard wide');
    c.appendChild(el('h3', null, 'Modality mix'));
    c.appendChild(el('p', 'sub', 'Paired records reference both spatial and single-cell data in the same deposit.'));

    var order = ['paired', 'spatial', 'singlecell'];
    var colors = ['var(--series-1)', 'var(--series-2)', 'var(--series-3)'];
    var counts = order.map(function (k) {
      return rows.filter(function (r) { return r.modality === k; }).length;
    });
    var total = counts.reduce(function (a, b) { return a + b; }, 0) || 1;

    var bar = el('div', 'stackbar');
    order.forEach(function (k, i) {
      if (!counts[i]) return;
      var seg = el('div');
      seg.style.width = (counts[i] / total * 100) + '%';
      seg.style.background = colors[i];
      seg.title = MODALITY_LABEL[k] + ': ' + nf(counts[i]) +
        ' (' + (counts[i] / total * 100).toFixed(1) + '%)';
      if (counts[i] / total > 0.08) seg.textContent = Math.round(counts[i] / total * 100) + '%';
      bar.appendChild(seg);
    });
    c.appendChild(bar);

    var legend = el('div', 'legend');
    order.forEach(function (k, i) {
      var s = el('span');
      var sw = el('i');
      sw.style.background = colors[i];
      s.appendChild(sw);
      s.appendChild(document.createTextNode(MODALITY_LABEL[k] + ' — ' + nf(counts[i])));
      legend.appendChild(s);
    });
    c.appendChild(legend);
    return c;
  }

  function renderStats() {
    var rows = lastFiltered.length || anyFilterActive() ? lastFiltered : DATA;
    var total = rows.length;

    var tiles = $('#tiles');
    tiles.textContent = '';
    var paired = rows.filter(function (r) { return r.modality === 'paired'; }).length;
    var imaging = rows.filter(function (r) { return r.platform_class === 'imaging'; }).length;
    var cancer = rows.filter(function (r) { return r.disease_group === 'Cancer'; }).length;
    var human = rows.filter(function (r) { return (r.organism || '').indexOf('Human') === 0 || r.organism === 'Human & Mouse'; }).length;
    var withCells = rows.filter(function (r) { return r.cells > 0; });
    var medianCells = 0;
    if (withCells.length) {
      var sorted = withCells.map(function (r) { return r.cells; }).sort(function (a, b) { return a - b; });
      medianCells = sorted[Math.floor(sorted.length / 2)];
    }

    tiles.appendChild(tile('Datasets in view', nf(total), META.total ? 'of ' + nf(META.total) + ' in catalog' : ''));
    tiles.appendChild(tile('Paired spatial + sc', nf(paired), total ? (paired / total * 100).toFixed(0) + '% of view' : ''));
    tiles.appendChild(tile('Imaging-based', nf(imaging), total ? (imaging / total * 100).toFixed(0) + '% of view' : ''));
    tiles.appendChild(tile('Cancer studies', nf(cancer), total ? (cancer / total * 100).toFixed(0) + '% of view' : ''));
    tiles.appendChild(tile('Human', nf(human), total ? (human / total * 100).toFixed(0) + '% of view' : ''));
    tiles.appendChild(tile('Median cells', medianCells ? compact(medianCells) : '—',
      withCells.length ? 'where a count was stated (' + nf(withCells.length) + ')' : 'no counts stated'));

    var grid = $('#chartgrid');
    grid.textContent = '';
    grid.appendChild(yearChart(rows));
    grid.appendChild(modalityChart(rows));
    grid.appendChild(hbarChart('Platforms', 'Detected from the record text.', tally(rows, 'platform'), total, 12));
    grid.appendChild(hbarChart('Disease / cancer type', 'Auto-classified; the pan-cancer bucket holds records with no specific type named.', tally(rows, 'disease'), total, 14));
    grid.appendChild(hbarChart('Tissue / organ', 'Best-scoring tissue keyword per record.', tally(rows, 'tissue'), total, 12));
    grid.appendChild(hbarChart('Source database', 'Where the record was harvested from.', tally(rows, 'db'), total, 6));
  }

  function anyFilterActive() {
    var any = F.q || F.yearMin != null || F.yearMax != null || F.cellsMin != null;
    GROUPS.forEach(function (g) { if (F[g.key].size) any = true; });
    return !!any;
  }

  // ---------------------------------------------------------------- views --
  function switchView(view) {
    $$('.tab').forEach(function (t) {
      var on = t.getAttribute('data-view') === view;
      t.classList.toggle('is-active', on);
      t.setAttribute('aria-selected', on ? 'true' : 'false');
    });
    $$('.view').forEach(function (v) {
      v.classList.toggle('is-active', v.id === 'view-' + view);
    });
    if (view === 'stats') renderStats();
  }

  // --------------------------------------------------------------- export --
  function csvEscape(v) {
    v = v == null ? '' : String(v);
    return /[",\n]/.test(v) ? '"' + v.replace(/"/g, '""') + '"' : v;
  }

  function exportCsv() {
    var cols = ['accession', 'db', 'title', 'modality', 'platform', 'platform_class', 'platform_source',
                'sc_platform', 'pairing', 'siblings', 'disease', 'disease_group', 'subtypes', 'study_type',
                'tissue', 'organism',
                'samples', 'cells', 'year', 'date', 'pmid', 'access', 'url'];
    var lines = [cols.join(',')];
    lastFiltered.forEach(function (r) {
      lines.push(cols.map(function (c) {
        var v = r[c];
        return csvEscape(Array.isArray(v) ? v.join('; ') : v);
      }).join(','));
    });
    download('spatial-sc-atlas-' + lastFiltered.length + '.csv', lines.join('\n'), 'text/csv;charset=utf-8');
    toast(nf(lastFiltered.length) + ' rows exported');
  }

  function exportJson() {
    var clean = lastFiltered.map(function (r) {
      var o = {};
      Object.keys(r).forEach(function (k) { if (k[0] !== '_') o[k] = r[k]; });
      return o;
    });
    download('spatial-sc-atlas-' + clean.length + '.json', JSON.stringify(clean, null, 2), 'application/json');
    toast(nf(clean.length) + ' records exported');
  }

  // ------------------------------------------------------------------ init --
  function wire() {
    $$('.tab').forEach(function (t) {
      t.addEventListener('click', function () {
        switchView(t.getAttribute('data-view'));
        writeURL();
      });
    });

    var searchTimer;
    $('#search').addEventListener('input', function (ev) {
      clearTimeout(searchTimer);
      var v = ev.target.value.trim();
      searchTimer = setTimeout(function () { F.q = v; onFilterChange(); }, 180);
    });

    $('#sort').addEventListener('change', function (ev) {
      F.sort = ev.target.value;
      shown = PAGE;
      renderResults();
      writeURL();
    });

    $('#resetFilters').addEventListener('click', resetFilters);

    $('#filtersToggle').addEventListener('click', function () {
      var open = $('#sidebar').classList.toggle('is-open');
      $('#filtersToggle').setAttribute('aria-expanded', open ? 'true' : 'false');
    });

    $('#exportCsv').addEventListener('click', exportCsv);
    $('#exportJson').addEventListener('click', exportJson);
    $('#copyAcc').addEventListener('click', function () {
      if (!lastFiltered.length) return toast('Nothing to copy');
      copyText(lastFiltered.map(function (r) { return r.accession; }).join('\n'),
               nf(lastFiltered.length) + ' accessions copied');
    });

    $('#openCompare').addEventListener('click', openCompare);
    $('#clearCompare').addEventListener('click', function () {
      picked.clear();
      updateTray();
      renderResults();
    });

    $('#closePanel').addEventListener('click', function () { $('#overlay').hidden = true; });
    $('#overlay').addEventListener('click', function (ev) {
      if (ev.target === $('#overlay')) $('#overlay').hidden = true;
    });
    document.addEventListener('keydown', function (ev) {
      if (ev.key === 'Escape') $('#overlay').hidden = true;
      var tag = (document.activeElement && document.activeElement.tagName) || '';
      if (ev.key === '/' && !/^(INPUT|TEXTAREA|SELECT)$/.test(tag)) {
        ev.preventDefault();
        $('#search').focus();
      }
    });

    var toggle = $('#themeToggle');
    var saved = null;
    try { saved = localStorage.getItem('atlas-theme'); } catch (e) { /* storage may be blocked */ }
    if (saved === 'dark' || saved === 'light') document.documentElement.setAttribute('data-theme', saved);
    toggle.addEventListener('click', function () {
      var cur = document.documentElement.getAttribute('data-theme');
      var isDark = cur === 'dark' ||
        (cur !== 'light' && window.matchMedia('(prefers-color-scheme: dark)').matches);
      var next = isDark ? 'light' : 'dark';
      document.documentElement.setAttribute('data-theme', next);
      try { localStorage.setItem('atlas-theme', next); } catch (e) { /* ignore */ }
    });
  }

  function load() {
    fetch('data/datasets.json', { cache: 'no-cache' })
      .then(function (r) {
        if (!r.ok) throw new Error('HTTP ' + r.status);
        return r.json();
      })
      .then(function (payload) {
        DATA = payload.datasets || [];
        META = payload.meta || {};
        byAccession = Object.create(null);
        DATA.forEach(function (r) { byAccession[r.accession] = r; });
        $('#stamp').textContent = META.generated ? 'Updated ' + META.generated : '';
        readURL();
        buildFilters();
        renderResults();
        renderStats();
      })
      .catch(function () {
        DATA = [];
        META = {};
        $('#stamp').textContent = '';
        buildFilters();
        renderResults();
        renderStats();
      });
  }

  wire();
  load();
})();
