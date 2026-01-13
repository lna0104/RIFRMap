// AI assistance disclosure:
// This script was written/refined with the help of an AI language model (ChatGPT).
//
// Hover overlay on base R plot image (Shiny plotOutput):
// - Horizontal bands: highlight POS / REF / MUT per SPECIES row (ignore clusters)
// - Vertical band: snap to AA columns within clusters (respects dxBlocks gaps)
//   using xcols + dx sent from R.
//
// Expected Shiny message payload:
//   { bands, xcols, dx, plt }
//
// bands columns: species, pos_ymin,pos_ymax, ref_ymin,ref_ymax, mut_ymin,mut_ymax
// xcols columns: cluster, x0, n, xmin, xmax

(function () {
  "use strict";

  // ---------- state ----------
  let SPECIES_BANDS = [];  // [{species, pos:{ymin,ymax}, ref:{}, mut:{}, total:{}}]
  let XCOLS = [];          // [{cluster,x0,n,xmin,xmax}, ...] plot coords (0..1)
  let DX = null;           // AA column step in plot coords (0..1)
  let PLT = null;          // par("plt") fractions: [x1, x2, y1, y2]
  let HOVER_ENABLED = false;

  // cached plot-region px (recomputed only on resize)
  let PR = null;

  // ---------- small utils ----------
  function byId(id) {
    return document.getElementById(id);
  }

  // Robust row decoder:
  // 1) row-wise array: [{...}, ...]
  // 2) column-wise object: {col1:[...], col2:[...]}
  // 3) row-object map: {"0":{...},"1":{...}}
  function toRows(payload) {
    if (!payload) return [];
    if (Array.isArray(payload)) return payload;

    if (payload && typeof payload === "object") {
      const keys = Object.keys(payload);
      if (!keys.length) return [];

      // Case 3: {"0": {...}, "1": {...}}
      if (keys.every(k => String(+k) === k)) {
        return keys
          .map(k => +k)
          .sort((a, b) => a - b)
          .map(k => payload[String(k)]);
      }

      // Case 2: column-wise object: {col: [..], ...}
      const colKeys = keys.filter(k => Array.isArray(payload[k]));
      if (!colKeys.length) return [];

      const n = payload[colKeys[0]].length;
      const out = new Array(n);

      for (let i = 0; i < n; i++) {
        const row = {};
        for (const k of colKeys) row[k] = payload[k][i];
        out[i] = row;
      }
      return out;
    }

    return [];
  }

  function finiteBand(ymin, ymax) {
    ymin = +ymin; ymax = +ymax;
    if (!Number.isFinite(ymin) || !Number.isFinite(ymax)) return null;
    if (ymax < ymin) { const t = ymin; ymin = ymax; ymax = t; }
    return { ymin, ymax };
  }

  function buildSpeciesBandsFromBandsDf(rows) {
    const out = [];
    for (const r of rows) {
      if (!r || r.species == null) continue;

      const pos = finiteBand(r.pos_ymin, r.pos_ymax);
      const ref = finiteBand(r.ref_ymin, r.ref_ymax);
      const mut = finiteBand(r.mut_ymin, r.mut_ymax);

      const ys = [];
      if (pos) ys.push(pos.ymin, pos.ymax);
      if (ref) ys.push(ref.ymin, ref.ymax);
      if (mut) ys.push(mut.ymin, mut.ymax);
      if (!ys.length) continue;

      out.push({
        species: String(r.species),
        pos, ref, mut,
        total: { ymin: Math.min(...ys), ymax: Math.max(...ys) }
      });
    }
    out.sort((a, b) => b.total.ymax - a.total.ymax); // top → bottom
    return out;
  }

  function buildXCols(rows) {
    const out = [];
    for (const r of rows) {
      if (!r) continue;
      const x0 = +r.x0, n = +r.n, xmin = +r.xmin, xmax = +r.xmax;
      if (![x0, n, xmin, xmax].every(Number.isFinite)) continue;
      if (n <= 0) continue;
      out.push({
        cluster: r.cluster != null ? String(r.cluster) : "",
        x0,
        n: Math.trunc(n),
        xmin,
        xmax
      });
    }
    return out;
  }

  function findSpeciesAtY(yPlot) {
    // linear scan; keep it simple (species count usually small-ish)
    for (const s of SPECIES_BANDS) {
      if (yPlot >= s.total.ymin && yPlot <= s.total.ymax) return s;
    }
    return null;
  }

  function bandNameAtY(s, yPlot) {
    if (s.pos && yPlot >= s.pos.ymin && yPlot <= s.pos.ymax) return "pos";
    if (s.ref && yPlot >= s.ref.ymin && yPlot <= s.ref.ymax) return "ref";
    if (s.mut && yPlot >= s.mut.ymin && yPlot <= s.mut.ymax) return "mut";
    return null;
  }

  function findClusterAtX(xPlot) {
    for (const c of XCOLS) {
      if (xPlot >= c.xmin && xPlot <= c.xmax) return c;
    }
    return null;
  }

  // Snap xPlot (0..1) to the nearest AA column center (0..1)
  function snapXToAA(xPlot) {
    if (!Number.isFinite(xPlot) || !Number.isFinite(DX) || DX <= 0 || !XCOLS.length) {
      return null;
    }

    const c = findClusterAtX(xPlot);
    if (!c) return null;

    // index relative to x0
    const t = (xPlot - c.x0) / DX;
    let u = Math.round(t);
    u = Math.max(0, Math.min(c.n - 1, u));
    return c.x0 + u * DX;
  }

  function clearCanvas() {
    const canvas = byId("hoverCanvas");
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.clearRect(0, 0, canvas.width, canvas.height);
  }

  // ---------- Shiny message handler ----------
  if (window.Shiny) {
    Shiny.addCustomMessageHandler("hoverMap", function (msg) {
      // console.log(msg);

      if (msg === "NULL" || msg == null) {
        HOVER_ENABLED = false;
        SPECIES_BANDS = [];
        XCOLS = [];
        DX = null;
        PLT = null;
        clearCanvas();
        return;
      }

      const bandsPayload =
        (msg && msg.bands !== undefined) ? msg.bands :
        (msg && msg.rows  !== undefined) ? msg.rows  :
        msg;

      const xcolsPayload =
        (msg && msg.xcols !== undefined) ? msg.xcols :
        null;

      const bandsRows = toRows(bandsPayload);
      SPECIES_BANDS = buildSpeciesBandsFromBandsDf(bandsRows);

      const xcolsRows = toRows(xcolsPayload);
      XCOLS = buildXCols(xcolsRows);

      DX  = (msg && msg.dx  != null) ? +msg.dx  : null;
      PLT = (msg && msg.plt != null) ? msg.plt  : null;

      HOVER_ENABLED = SPECIES_BANDS.length > 0;

      // console.log("hoverMap:", "bands=", SPECIES_BANDS.length, "xcols=", XCOLS.length, "dx=", DX, "plt=", PLT);
    });
  }

  // ---------- hover drawing ----------
  function setupHoverOverlay() {
    const plotDiv = byId("clusterPlot");
    const canvas = byId("hoverCanvas");
    if (!plotDiv || !canvas) return false;

    const img = plotDiv.querySelector("img");
    if (!img) return false;

    const wrap = canvas.parentElement || plotDiv;
    const ctx = canvas.getContext("2d");
    if (!ctx) return false;

    // ---- sizing / DPR ----
    function plotRegionPx() {
      const W = img.offsetWidth;
      const H = img.offsetHeight;

      if (!PLT || !Array.isArray(PLT) || PLT.length !== 4) {
        return { left: 0, right: W, top: 0, bottom: H };
      }

      const x1 = +PLT[0], x2 = +PLT[1], y1 = +PLT[2], y2 = +PLT[3];
      if (![x1, x2, y1, y2].every(Number.isFinite)) {
        return { left: 0, right: W, top: 0, bottom: H };
      }

      return {
        left: x1 * W,
        right: x2 * W,
        top: (1 - y2) * H,
        bottom: (1 - y1) * H
      };
    }

    function syncCanvas() {
      const w = img.offsetWidth;
      const h = img.offsetHeight;

      // cap DPR for speed (retina can be 2 => 4x pixels)
      const dpr0 = window.devicePixelRatio || 1;
      const dpr = (dpr0 >= 2 ? 2 : 1);

      canvas.style.width = w + "px";
      canvas.style.height = h + "px";
      canvas.width = Math.max(1, Math.round(w * dpr));
      canvas.height = Math.max(1, Math.round(h * dpr));

      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

      PR = plotRegionPx(); // cache plot region
      ctx.clearRect(0, 0, canvas.width, canvas.height);
    }

    function getPR() {
      return PR || plotRegionPx();
    }

    // ---- overlay drawing (only clears plot region) ----
    function drawOverlay(pr, band, xPix) {
      const plotW = pr.right - pr.left;
      const plotH = pr.bottom - pr.top;

      // clear only plot region
      ctx.clearRect(pr.left - 2, pr.top - 2, plotW + 4, plotH + 4);

      // Horizontal band
      if (band) {
        const yTop = pr.top + (1 - band.ymax) * plotH;
        const yBot = pr.top + (1 - band.ymin) * plotH;

        ctx.fillStyle = "rgba(0,0,0,0.06)";
        ctx.fillRect(pr.left, yTop, plotW, yBot - yTop);

        ctx.strokeStyle = "rgba(0,0,0,0.18)";
        ctx.lineWidth = 1;
        ctx.strokeRect(pr.left + 0.5, yTop + 0.5, plotW - 1, (yBot - yTop) - 1);
      }

      // Vertical band width ≈ one AA column (dx) in pixels
      if (Number.isFinite(xPix)) {
        const bandW = (Number.isFinite(DX) && DX > 0)
          ? Math.max(3, DX * plotW * 0.8)   // ~80% of one AA column
          : 10;

        const x = Math.max(pr.left, Math.min(pr.right, xPix));
        const xLeft = x - bandW / 2;

        ctx.fillStyle = "rgba(0,0,0,0.06)";
        ctx.fillRect(xLeft, pr.top, bandW, plotH);

        ctx.strokeStyle = "rgba(0,0,0,0.18)";
        ctx.lineWidth = 1;
        ctx.strokeRect(
          Math.round(xLeft) + 0.5,
          Math.round(pr.top) + 0.5,
          Math.max(1, Math.round(bandW) - 1),
          Math.max(1, Math.round(plotH) - 1)
        );
      }
    }

    // ---- state for throttled drawing ----
    let lastKey = null;     // species::zone
    let lastBand = null;    // {ymin,ymax} in plot coords
    let lastXPlot = null;   // snapped x in plot coords
    let lastXpx = null;     // snapped x in CSS pixels

    // draw scheduling (single RAF for BOTH x and y changes)
    let pending = null; // {band, xpx}
    let rafDraw = false;

    function requestDraw(pr, band, xpx) {
      // avoid redundant draws
      const sameBand =
        (band === lastBand) ||
        (band && lastBand && band.ymin === lastBand.ymin && band.ymax === lastBand.ymax);

      const sameX = (xpx === lastXpx);

      if (sameBand && sameX) return;

      pending = { pr, band, xpx };
      if (rafDraw) return;

      rafDraw = true;
      requestAnimationFrame(() => {
        rafDraw = false;
        if (!pending) return;

        lastBand = pending.band;
        lastXpx  = pending.xpx;

        drawOverlay(pending.pr, lastBand, lastXpx);
        pending = null;
      });
    }

    function updateBandFromY(pr, yPlot) {
      if (!HOVER_ENABLED || !SPECIES_BANDS.length) {
        lastKey = null;
        requestDraw(pr, null, lastXpx);
        return;
      }

      const s = findSpeciesAtY(yPlot);
      if (!s) {
        lastKey = null;
        requestDraw(pr, null, lastXpx);
        return;
      }

      const zone = bandNameAtY(s, yPlot);
      if (!zone || !s[zone]) {
        lastKey = null;
        requestDraw(pr, null, lastXpx);
        return;
      }

      const key = `${s.species}::${zone}`;
      if (key === lastKey) return;

      lastKey = key;
      requestDraw(pr, s[zone], lastXpx);
    }

    // ---- mouse events ----
    wrap.addEventListener("mousemove", function (e) {
      if (!HOVER_ENABLED || !SPECIES_BANDS.length) {
        const pr0 = getPR();
        ctx.clearRect(pr0.left, pr0.top, pr0.right - pr0.left, pr0.bottom - pr0.top);
        return;
      }

      const rect = img.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;

      const pr = getPR();
      if (mx < pr.left || mx > pr.right || my < pr.top || my > pr.bottom) {
        lastKey = null;
        lastBand = null;
        lastXPlot = null;
        lastXpx = null;
        ctx.clearRect(pr.left, pr.top, pr.right - pr.left, pr.bottom - pr.top);
        return;
      }

      const plotW = pr.right - pr.left;
      const plotH = pr.bottom - pr.top;

      // mouse → plot coords (0..1)
      const xPlot = (mx - pr.left) / plotW;
      const yPlot = 1 - ((my - pr.top) / plotH);

      // snap X (cheap) and only update if snapped col changes
      const xSnapPlot = snapXToAA(xPlot);
      if (xSnapPlot !== lastXPlot) {
        lastXPlot = xSnapPlot;
        lastXpx = (xSnapPlot == null) ? null : (pr.left + xSnapPlot * plotW);
        requestDraw(pr, lastBand, lastXpx);
      }

      // update band (also draws via RAF only when it changes)
      updateBandFromY(pr, yPlot);
    });

    wrap.addEventListener("mouseleave", function () {
      const pr = getPR();
      lastKey = null;
      lastBand = null;
      lastXPlot = null;
      lastXpx = null;
      ctx.clearRect(pr.left, pr.top, pr.right - pr.left, pr.bottom - pr.top);
    });

    // Keep canvas aligned with image resizing
    new ResizeObserver(syncCanvas).observe(img);
    window.addEventListener("resize", syncCanvas);
    img.addEventListener("load", syncCanvas);

    syncCanvas();
    return true;
  }

  function start() {
    if (setupHoverOverlay()) return;

    const root = document.body || document.documentElement;
    if (!root) return;

    const obs = new MutationObserver(() => {
      if (setupHoverOverlay()) obs.disconnect();
    });

    obs.observe(root, { childList: true, subtree: true });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", start, { once: true });
  } else {
    start();
  }
})();
