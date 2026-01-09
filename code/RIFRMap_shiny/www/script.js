// AI assistance disclosure:
// This script were written/refined with the help of an AI language model (ChatGPT).
// Hover overlay on base R plot image (Shiny plotOutput):
// - Horizontal bands: highlight POS / REF / MUT per SPECIES row (ignore clusters)
// - Vertical band: fixed-width filled band + border, same style as horizontal band

(function () {
  "use strict";

  // ---------- state ----------
  let SPECIES_BANDS = []; // [{species, pos:{ymin,ymax}, ref:{}, mut:{}, total:{}}]
  let PLT = null;        // R plot region fractions: [x1, x2, y1, y2]
  let HOVER_ENABLED = false;

  // ---------- utilities ----------
  function byId(id) {
    return document.getElementById(id);
  }

  function toRows(payload) {
    // Accept either:
    // 1) row-wise array: [{species, kind, xmin, xmax, ymin, ymax}, ...]
    // 2) column-wise object: {species:[...], kind:[...], xmin:[...], ...}
    if (Array.isArray(payload)) return payload;

    if (payload && typeof payload === "object" && Array.isArray(payload.species)) {
      const n = payload.species.length;
      const out = new Array(n);
      for (let i = 0; i < n; i++) {
        out[i] = {
          species: payload.species[i],
          kind:    payload.kind?.[i] ?? null,
          xmin:    payload.xmin ? +payload.xmin[i] : NaN,
          xmax:    payload.xmax ? +payload.xmax[i] : NaN,
          ymin:    payload.ymin ? +payload.ymin[i] : NaN,
          ymax:    payload.ymax ? +payload.ymax[i] : NaN
        };
      }
      return out;
    }
    return [];
  }

  function bandInit() {
    return { ymin: +Infinity, ymax: -Infinity };
  }

  function bandAdd(b, r) {
    const ymin = +r.ymin, ymax = +r.ymax;
    if (!Number.isFinite(ymin) || !Number.isFinite(ymax)) return;
    b.ymin = Math.min(b.ymin, ymin);
    b.ymax = Math.max(b.ymax, ymax);
  }

  function buildSpeciesBands(rows) {
    const map = new Map();

    for (const r of rows) {
      if (!r || !r.species || !r.kind) continue;

      let s = map.get(r.species);
      if (!s) {
        s = { species: r.species, pos: bandInit(), ref: bandInit(), mut: bandInit() };
        map.set(r.species, s);
      }

      if (r.kind === "pos") bandAdd(s.pos, r);
      else if (r.kind === "ref") bandAdd(s.ref, r);
      else if (r.kind === "mut") bandAdd(s.mut, r);
    }

    const out = [];
    for (const s of map.values()) {
      if (!Number.isFinite(s.pos.ymin)) s.pos = null;
      if (!Number.isFinite(s.ref.ymin)) s.ref = null;
      if (!Number.isFinite(s.mut.ymin)) s.mut = null;

      const ys = [];
      if (s.pos) ys.push(s.pos.ymin, s.pos.ymax);
      if (s.ref) ys.push(s.ref.ymin, s.ref.ymax);
      if (s.mut) ys.push(s.mut.ymin, s.mut.ymax);
      if (!ys.length) continue;

      s.total = { ymin: Math.min(...ys), ymax: Math.max(...ys) };
      out.push(s);
    }

    // sort top → bottom
    out.sort((a, b) => b.total.ymax - a.total.ymax);
    return out;
  }

  function findSpeciesAtY(bands, y) {
    // y is plot coord (0..1), y increases upward
    for (const s of bands) {
      if (y >= s.total.ymin && y <= s.total.ymax) return s;
    }
    return null;
  }

  function bandNameAtY(s, y) {
    // fixed zones within each species: POS / REF / MUT
    if (s.pos && y >= s.pos.ymin && y <= s.pos.ymax) return "pos";
    if (s.ref && y >= s.ref.ymin && y <= s.ref.ymax) return "ref";
    if (s.mut && y >= s.mut.ymin && y <= s.mut.ymax) return "mut";
    return null;
  }
  
  function clearCanvas() {
    const canvas = document.getElementById("hoverCanvas");
    if (!canvas) return;
  
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
  
    ctx.clearRect(0, 0, canvas.width, canvas.height);
  }

  // ---------- Shiny message handler ----------
  if (window.Shiny) {
    Shiny.addCustomMessageHandler("hoverMap", function (msg) {
      if (msg == "NULL") {
        HOVER_ENABLED = false;
        SPECIES_BANDS = [];
        PLT = null;
  
        clearCanvas();  
  
        console.log("hoverMap: OFF");
        return;
      }
      
      HOVER_ENABLED = true;
      const rowsPayload = (msg && msg.rows !== undefined) ? msg.rows : msg;
      const rows = toRows(rowsPayload);

      SPECIES_BANDS = buildSpeciesBands(rows);
      PLT = (msg && msg.plt) ? msg.plt : null;

      console.log("hoverMap:", "rows =", rows.length, "species =", SPECIES_BANDS.length, "plt =", PLT);
    });
  }

  // ---------- hover drawing ----------
  function setupHoverOverlay() {
    const plotDiv = byId("clusterPlot");
    const canvas  = byId("hoverCanvas");
    if (!plotDiv || !canvas) return false;

    // Shiny plotOutput renders an <img>
    const img = plotDiv.querySelector("img");
    if (!img) return false;

    // overlay should sit in the same wrapper as the img (position:relative on wrapper)
    const wrap = canvas.parentElement || plotDiv;
    const ctx  = canvas.getContext("2d");

    function clear() {
      // canvas is scaled to DPR but ctx is transformed to CSS pixels,
      // so clear in device pixels for correctness:
      ctx.clearRect(0, 0, canvas.width, canvas.height);
    }

    function syncCanvas() {
      const w = img.offsetWidth;
      const h = img.offsetHeight;
      const dpr = window.devicePixelRatio || 1;

      canvas.style.width  = w + "px";
      canvas.style.height = h + "px";
      canvas.width  = Math.max(1, Math.round(w * dpr));
      canvas.height = Math.max(1, Math.round(h * dpr));

      // After this, all drawing coords are in CSS pixels, and lineWidth is fixed in CSS pixels.
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      clear();
    }

    function plotRegionPx() {
      const W = img.offsetWidth;
      const H = img.offsetHeight;

      // If PLT not available, use entire image area
      if (!PLT || !Array.isArray(PLT) || PLT.length !== 4) {
        return { left: 0, right: W, top: 0, bottom: H };
      }

      const x1 = +PLT[0], x2 = +PLT[1], y1 = +PLT[2], y2 = +PLT[3];
      if (![x1, x2, y1, y2].every(Number.isFinite)) {
        return { left: 0, right: W, top: 0, bottom: H };
      }

      return {
        left:   x1 * W,
        right:  x2 * W,
        top:    (1 - y2) * H,
        bottom: (1 - y1) * H
      };
    }

    // ---- overlay drawer: horizontal band + vertical band ----
    let lastBand = null;  // current POS/REF/MUT band object (ymin/ymax)
    let lastXpx  = null;  // current x in CSS pixels, relative to image

    function drawOverlay(band, xPix) {
      const pr = plotRegionPx();
      const plotW = pr.right - pr.left;
      const plotH = pr.bottom - pr.top;

      clear();

      // Horizontal band (same style)
      if (band) {
        const yTop = pr.top + (1 - band.ymax) * plotH;
        const yBot = pr.top + (1 - band.ymin) * plotH;

        ctx.fillStyle = "rgba(0,0,0,0.06)";
        ctx.fillRect(pr.left, yTop, plotW, yBot - yTop);

        ctx.strokeStyle = "rgba(0,0,0,0.18)";
        ctx.lineWidth = 1;
        ctx.strokeRect(pr.left + 0.5, yTop + 0.5, plotW - 1, (yBot - yTop) - 1);
      }

      // Vertical band (fixed width)
      if (Number.isFinite(xPix)) {
        const BAND_W = 25; // ✅ fixed width in CSS pixels (change to taste)
        const x = Math.max(pr.left, Math.min(pr.right, xPix));
        const xLeft = x - BAND_W / 2;

        ctx.fillStyle = "rgba(0,0,0,0.06)";
        ctx.fillRect(xLeft, pr.top, BAND_W, plotH);

        ctx.strokeStyle = "rgba(0,0,0,0.18)";
        ctx.lineWidth = 1;
        ctx.strokeRect(
          Math.round(xLeft) + 0.5,
          Math.round(pr.top) + 0.5,
          BAND_W - 1,
          plotH - 1
        );
      }
    }

    // ---- band selection (throttled) ----
    let lastKey = null;
    let pendingY = null;
    let raf = false;

    function handleMove(yPlot) {
      if (!SPECIES_BANDS.length) {
        lastKey = null;
        lastBand = null;
        drawOverlay(null, lastXpx);
        return;
      }

      const s = findSpeciesAtY(SPECIES_BANDS, yPlot);
      if (!s) {
        lastKey = null;
        lastBand = null;
        drawOverlay(null, lastXpx);
        return;
      }

      const name = bandNameAtY(s, yPlot);
      if (!name || !s[name]) {
        lastKey = null;
        lastBand = null;
        drawOverlay(null, lastXpx);
        return;
      }

      const key = `${s.species}::${name}`;
      if (key === lastKey) return;

      lastKey = key;
      lastBand = s[name];
      drawOverlay(lastBand, lastXpx);
    }

    function schedule(yPlot) {
      pendingY = yPlot;
      if (raf) return;
      raf = true;

      requestAnimationFrame(() => {
        raf = false;
        if (pendingY != null) handleMove(pendingY);
      });
    }

    // ---- mouse events ----
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
        lastXpx = null;
        clear();
        return;
      }

      // store x for vertical band
      lastXpx = x;

      // Convert y within plot region -> plot coords (0..1), y up
      const frac = (y - pr.top) / (pr.bottom - pr.top);
      const yPlot = 1 - frac;

      // schedule band updates (throttled)
      schedule(yPlot);

      // draw immediately so vertical band follows mouse smoothly even if band unchanged
      drawOverlay(lastBand, lastXpx);
    });

    wrap.addEventListener("mouseleave", function () {
      lastKey = null;
      lastBand = null;
      lastXpx = null;
      clear();
    });

    // keep canvas aligned with image resizing
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
