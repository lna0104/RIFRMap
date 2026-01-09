// AI assistance disclosure:
// This script was written/refined with the help of an AI language model (ChatGPT).
//
// Hover overlay on base R plot image (Shiny plotOutput):
// - Horizontal bands: highlight POS / REF / MUT per SPECIES row (ignore clusters)
// - Vertical band: snap to AA columns within clusters (respects dxBlocks gaps)
//   using xcols + dx sent from R.
//
// Expected Shiny message payload:
//   {
//     bands: <data.frame/tibble of per-species y bands>,
//     xcols: <data.frame/tibble of per-cluster x regions>,
//     dx:    <numeric step in plot coords>,
//     plt:   <par("plt") fractions: [x1, x2, y1, y2]>
//   }
//
// bands columns: species, pos_ymin,pos_ymax, ref_ymin,ref_ymax, mut_ymin,mut_ymax
// xcols columns: cluster, x0, n, xmin, xmax

(function () {
  "use strict";

  // ---------- state ----------
  let SPECIES_BANDS = [];  // [{species, pos:{ymin,ymax}, ref:{}, mut:{}, total:{}}]
  let XCOLS = [];          // [{cluster,x0,n,xmin,xmax}, ...] in plot coords (0..1)
  let DX = null;           // AA column step in plot coords (0..1)
  let PLT = null;          // R plot region fractions: [x1, x2, y1, y2]
  let HOVER_ENABLED = false;

  // ---------- small utils ----------
  function byId(id) {
    return document.getElementById(id);
  }

  function isPlainObject(x) {
    return x && typeof x === "object" && !Array.isArray(x);
  }

  function looksLikeRowObjectMap(obj) {
    // Shiny sometimes serializes df as {"0":{...},"1":{...}}
    if (!isPlainObject(obj)) return false;
    const keys = Object.keys(obj);
    if (!keys.length) return false;
    return keys.every(k => String(+k) === k);
  }

  function clearCanvas() {
    const canvas = byId("hoverCanvas");
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.clearRect(0, 0, canvas.width, canvas.height);
  }

  // Robust row decoder:
  // 1) row-wise array: [{...}, ...]
  // 2) column-wise object: {col1:[...], col2:[...]}
  // 3) row-object map: {"0":{...},"1":{...}}
  function toRows(payload) {
    if (!payload) return [];
    if (Array.isArray(payload)) return payload;
  
    // Case 3: {"0": {...}, "1": {...}}
    if (payload && typeof payload === "object" && !Array.isArray(payload)) {
      const keys = Object.keys(payload);
  
      // row-object map (all numeric keys)
      if (keys.length && keys.every(k => String(+k) === k)) {
        return keys
          .map(k => +k)
          .sort((a, b) => a - b)
          .map(k => payload[String(k)]);
      }
  
      // Case 2: column-wise object: {col: [..], ...}
      // Find array-like columns
      const colKeys = keys.filter(k => Array.isArray(payload[k]));
      if (!colKeys.length) return [];
  
      // Determine row count from first array column
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
      if (!r || !r.species) continue;

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
    // sort top → bottom
    out.sort((a, b) => b.total.ymax - a.total.ymax);
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
        x0, n: Math.trunc(n), xmin, xmax
      });
    }
    return out;
  }

  function findSpeciesAtY(yPlot) {
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
  
    // floating index relative to x0
    const t = (xPlot - c.x0) / DX;
  
    // snap to nearest AA column
    let u = Math.round(t);
  
    // hard clamp to [0, n-1]
    u = Math.max(0, Math.min(c.n - 1, u));
  
    // snapped AA column center
    return c.x0 + u * DX;
  }

  // ---------- Shiny message handler ----------
  if (window.Shiny) {
    Shiny.addCustomMessageHandler("hoverMap", function (msg) {
      console.log(msg);
  
      if (msg === "NULL" || msg == null) {
        HOVER_ENABLED = false;
        SPECIES_BANDS = [];
        XCOLS = [];
        DX = null;
        PLT = null;
        clearCanvas();
        return;
      }
  
      // expected {bands, xcols, dx, plt}
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
  
      // enable hover if we at least have horizontal bands
      HOVER_ENABLED = SPECIES_BANDS.length > 0;
  
      console.log(
        "hoverMap:",
        "bands rows =", bandsRows.length,
        "species =", SPECIES_BANDS.length,
        "xcols rows =", xcolsRows.length,
        "xcols =", XCOLS.length,
        "dx =", DX,
        "plt =", PLT
      );
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

    // Clear in device pixels (canvas is DPR-scaled, but ctx is transformed to CSS pixels)
    function clear() {
      ctx.clearRect(0, 0, canvas.width, canvas.height);
    }

    function syncCanvas() {
      const w = img.offsetWidth;
      const h = img.offsetHeight;
      const dpr = window.devicePixelRatio || 1;

      canvas.style.width = w + "px";
      canvas.style.height = h + "px";
      canvas.width = Math.max(1, Math.round(w * dpr));
      canvas.height = Math.max(1, Math.round(h * dpr));

      // Draw in CSS pixels
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      clear();
    }

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

    // ---- overlay drawer ----
    let lastBand = null;     // {ymin,ymax} in plot coords
    let lastXPlot = null;    // snapped x in plot coords
    let lastXpx = null;      // snapped x in CSS pixels
    let lastKey = null;      // species::zone
    let raf = false;
    let pendingY = null;

    function drawOverlay(band, xPix) {
      const pr = plotRegionPx();
      const plotW = pr.right - pr.left;
      const plotH = pr.bottom - pr.top;

      clear();

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

    function handleY(yPlot) {
      if (!HOVER_ENABLED || !SPECIES_BANDS.length) {
        lastKey = null;
        lastBand = null;
        drawOverlay(null, lastXpx);
        return;
      }

      const s = findSpeciesAtY(yPlot);
      if (!s) {
        lastKey = null;
        lastBand = null;
        drawOverlay(null, lastXpx);
        return;
      }

      const zone = bandNameAtY(s, yPlot);
      if (!zone || !s[zone]) {
        lastKey = null;
        lastBand = null;
        drawOverlay(null, lastXpx);
        return;
      }

      const key = `${s.species}::${zone}`;
      if (key === lastKey) return;

      lastKey = key;
      lastBand = s[zone];
      drawOverlay(lastBand, lastXpx);
    }

    function scheduleY(yPlot) {
      pendingY = yPlot;
      if (raf) return;
      raf = true;

      requestAnimationFrame(() => {
        raf = false;
        if (pendingY != null) handleY(pendingY);
      });
    }

    // ---- mouse events (snap X using xcols + dx; throttle Y with rAF) ----
    wrap.addEventListener("mousemove", function (e) {
      if (!HOVER_ENABLED || !SPECIES_BANDS.length) {
        clear();
        return;
      }

      const rect = img.getBoundingClientRect();
      const x = e.clientX - rect.left;
      const y = e.clientY - rect.top;

      const pr = plotRegionPx();
      if (x < pr.left || x > pr.right || y < pr.top || y > pr.bottom) {
        lastKey = null;
        lastBand = null;
        lastXPlot = null;
        lastXpx = null;
        clear();
        return;
      }

      const plotW = pr.right - pr.left;
      const plotH = pr.bottom - pr.top;

      // Convert mouse → plot coords (0..1)
      const xPlot = (x - pr.left) / plotW;
      const yPlot = 1 - ((y - pr.top) / plotH);

      // Snap X to AA column centers using xcols + dx
      const xSnapPlot = snapXToAA(xPlot);

      // Redraw ONLY when snapped column changes
      if (xSnapPlot !== lastXPlot) {
        lastXPlot = xSnapPlot;

        if (xSnapPlot == null) {
          lastXpx = null;
        } else {
          lastXpx = pr.left + xSnapPlot * plotW;
        }

        drawOverlay(lastBand, lastXpx);
      }

      // Band updates throttled (only redraw when band changes)
      scheduleY(yPlot);
    });

    wrap.addEventListener("mouseleave", function () {
      lastKey = null;
      lastBand = null;
      lastXPlot = null;
      lastXpx = null;
      clear();
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
