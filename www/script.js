$(document).on("click", "#go_docs", function() {
  Shiny.setInputValue("go_docs", Math.random());
});

$(document).on("click", "#go_explorer", function() {
  Shiny.setInputValue("go_explorer", Math.random());
});

$(document).on("change", ".filter-box input[type='checkbox']", function() {
  const box = $(this).closest(".filter-box");
  box.toggleClass("active", $(this).is(":checked"));
});

$(function () {
  $('[data-toggle="tooltip"]').tooltip();
});

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

  // Currently bound DOM elements. Shiny's clusterPlotUI (server.R) re-renders
  // the whole plotOutput+canvas subtree (fresh <img>/<canvas> nodes) whenever
  // cluster_dims_r() changes, so "bound once at startup" is not safe -- we
  // must keep checking and rebind whenever the live elements differ.
  let CUR = { plotDiv: null, canvas: null, img: null, wrap: null, ctx: null };
  let ROBS = null; // ResizeObserver on CUR.img

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

  function clearCanvas(canvas) {
    canvas = canvas || CUR.canvas;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.clearRect(0, 0, canvas.width, canvas.height);
  }

  function resetInteractionState() {
    lastBand = null;
    lastXPlot = null;
    lastXpx = null;
    pending = null;
  }

  // ---------- plot-region math ----------
  // Always recomputed from the live <img> box + the current PLT -- never
  // cached across events. offsetWidth/offsetHeight reads are cheap, and
  // caching this previously meant a new render's `plt` margins (sent via the
  // hoverMap message) had no effect on hover math until an unrelated resize
  // happened to fire, which is the main source of the reported misalignment.
  function plotRegionPx() {
    const img = CUR.img;
    if (!img) return { left: 0, right: 0, top: 0, bottom: 0 };

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

  // Resize the canvas backing store to match the live <img> box/DPR. Called
  // on bind, on ResizeObserver/window resize, on <img> load, and whenever a
  // hoverMap message arrives (a new render can change `plt` without the
  // <img> box size changing at all).
  function syncCanvasBackingStore() {
    const { img, canvas, ctx } = CUR;
    if (!img || !canvas || !ctx) return;

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
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    resetInteractionState();
  }

  // ---------- overlay drawing ----------
  function drawOverlay(ctx, pr, band, xPix) {
    if (!ctx) return;
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
  let lastBand = null;    // {ymin,ymax} in plot coords
  let lastXPlot = null;   // snapped x in plot coords
  let lastXpx = null;     // snapped x in CSS pixels

  // draw scheduling (single RAF for BOTH x and y changes)
  let pending = null; // {ctx, pr, band, xpx}
  let rafDraw = false;

  function requestDraw(ctx, pr, band, xpx) {
    // avoid redundant draws
    const sameBand =
      (band === lastBand) ||
      (band && lastBand && band.ymin === lastBand.ymin && band.ymax === lastBand.ymax);

    const sameX = (xpx === lastXpx);

    if (sameBand && sameX) return;

    pending = { ctx, pr, band, xpx };
    if (rafDraw) return;

    rafDraw = true;
    requestAnimationFrame(() => {
      rafDraw = false;
      if (!pending) return;

      lastBand = pending.band;
      lastXpx  = pending.xpx;

      drawOverlay(pending.ctx, pending.pr, lastBand, lastXpx);
      pending = null;
    });
  }

  // ---------- mouse events (bound once per live element set, see bindHoverOverlay) ----------
  function onMouseMove(e) {
    const { img, ctx } = CUR;
    if (!img || !ctx) return;

    if (!HOVER_ENABLED || !SPECIES_BANDS.length) {
      const pr0 = plotRegionPx();
      ctx.clearRect(pr0.left, pr0.top, pr0.right - pr0.left, pr0.bottom - pr0.top);
      return;
    }

    const rect = img.getBoundingClientRect();
    const mx = e.clientX - rect.left;
    const my = e.clientY - rect.top;

    const pr = plotRegionPx();
    if (mx < pr.left || mx > pr.right || my < pr.top || my > pr.bottom) {
      resetInteractionState();
      ctx.clearRect(pr.left, pr.top, pr.right - pr.left, pr.bottom - pr.top);
      return;
    }

    const plotW = pr.right - pr.left;
    const plotH = pr.bottom - pr.top;

    // mouse → plot coords (0..1)
    const xPlot = (mx - pr.left) / plotW;
    const yPlot = 1 - ((my - pr.top) / plotH);

    // Resolve the snapped column and the row band from this single mouse
    // position, then issue ONE combined draw. Previously these were two
    // separate requestDraw() calls that each captured the other axis's
    // not-yet-committed value (lastXpx/lastBand are only updated inside the
    // deferred rAF callback) -- when both axes changed on the same event,
    // the second call clobbered the first and silently dropped one axis's
    // highlight for that frame.
    const xSnapPlot = snapXToAA(xPlot);
    const xpx = (xSnapPlot == null) ? null : (pr.left + xSnapPlot * plotW);

    let band = null;
    const s = findSpeciesAtY(yPlot);
    if (s) {
      const zone = bandNameAtY(s, yPlot);
      if (zone && s[zone]) band = s[zone];
    }

    lastXPlot = xSnapPlot;
    requestDraw(ctx, pr, band, xpx);
  }

  function onMouseLeave() {
    const { ctx } = CUR;
    const pr = plotRegionPx();
    resetInteractionState();
    if (ctx) ctx.clearRect(pr.left, pr.top, pr.right - pr.left, pr.bottom - pr.top);
  }

  // ---------- binding ----------
  // Idempotent: safe to call on every DOM mutation. Only does work when the
  // live <img>/<canvas> differ from what we're currently bound to, which
  // happens whenever clusterPlotUI (server.R) re-renders the plot subtree
  // (e.g. a filter change alters the plot's pixel dimensions).
  function bindHoverOverlay() {
    const plotDiv = byId("clusterPlot");
    const canvas = byId("hoverCanvas");
    if (!plotDiv || !canvas) return false;

    const img = plotDiv.querySelector("img");
    if (!img) return false;

    const wrap = canvas.parentElement || plotDiv;

    if (CUR.img === img && CUR.canvas === canvas && CUR.wrap === wrap) {
      return true; // already bound to the current live elements
    }

    const ctx = canvas.getContext("2d");
    if (!ctx) return false;

    if (CUR.wrap) {
      CUR.wrap.removeEventListener("mousemove", onMouseMove);
      CUR.wrap.removeEventListener("mouseleave", onMouseLeave);
    }
    if (ROBS) {
      ROBS.disconnect();
      ROBS = null;
    }

    CUR = { plotDiv, canvas, img, wrap, ctx };

    wrap.addEventListener("mousemove", onMouseMove);
    wrap.addEventListener("mouseleave", onMouseLeave);

    ROBS = new ResizeObserver(syncCanvasBackingStore);
    ROBS.observe(img);
    img.addEventListener("load", syncCanvasBackingStore);
    window.addEventListener("resize", syncCanvasBackingStore);

    syncCanvasBackingStore();
    return true;
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
        resetInteractionState();
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

      // A new render can swap in a fresh <img>/<canvas> (clusterPlotUI) and/or
      // change `plt` margins without changing the <img> box size, so rebind
      // defensively and force a full geometry resync rather than waiting for
      // a resize/load event that may not come.
      bindHoverOverlay();
      syncCanvasBackingStore();

      // console.log("hoverMap:", "bands=", SPECIES_BANDS.length, "xcols=", XCOLS.length, "dx=", DX, "plt=", PLT);
    });
  }

  // ---------- boot ----------
  // Keep watching for the lifetime of the page: clusterPlotUI can replace the
  // plot subtree at any time (not just once at startup), so a one-shot setup
  // is not sufficient.
  function start() {
    bindHoverOverlay();

    const root = document.body || document.documentElement;
    if (!root) return;

    let scheduled = false;
    const obs = new MutationObserver(() => {
      if (scheduled) return;
      scheduled = true;
      requestAnimationFrame(() => {
        scheduled = false;
        bindHoverOverlay();
      });
    });

    obs.observe(root, { childList: true, subtree: true });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", start, { once: true });
  } else {
    start();
  }
})();

