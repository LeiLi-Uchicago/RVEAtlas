/* =============================================================================
 * Pathogen / subtype switch overlay.
 *
 * A self-contained full-screen loading overlay, independent of the `waiter`
 * package (which is used by the per-output withWaiter wrappers — sharing it left
 * a stuck mask on some switches). The server only SHOWS it (custom message
 * `rveSwitchOverlayShow` with {title, subtitle, pathogen}); the client hides it
 * on Shiny's `shiny:idle` (recompute finished), with `shiny:error`,
 * `shiny:disconnected`, and a hard safety timer as fallbacks so it can never
 * stick. Each pathogen gets its own cute cartoon virus.
 * ============================================================================ */
(function () {
  var MIN_MS  = 600;   // keep the indicator on screen at least this long
  var IDLE_MS = 350;   // hide this long after Shiny stops recomputing
  var SAFETY_MS = 15000;
  var HOLD_MS = 1600;  // "hold" mode: wait up to this long for the busy that a
                       // debounced recompute triggers before allowing idle-hide
  var WAIT_MS = 6000;  // if the awaited output never updates, stop waiting after this
  var active = false, shownAt = 0;
  var safetyTimer = null, idleTimer = null, minTimer = null, holdTimer = null, waitTimer = null;
  // When set, idle-hide is suppressed until the named element's DOM content
  // actually changes — used so the mask stays up until, e.g., the Home pathogen
  // text refreshes, which lands a flush (or two) after shiny:idle. We watch the
  // element directly with a MutationObserver (robust to event-timing quirks),
  // and also clear on the matching shiny:value as a backup.
  var waitFor = null, waiting = false, waitObserver = null, waitInitialText = "";
  function stopWaitObserver() {
    if (waitObserver) { waitObserver.disconnect(); waitObserver = null; }
  }
  function contentUpdated() {
    if (!waiting) return;
    waiting = false;
    stopWaitObserver();
    if (waitTimer) { clearTimeout(waitTimer); waitTimer = null; }
    armIdleHide();
  }
  // Only treat the awaited element as "ready" once it holds real, CHANGED text.
  // Shiny blanks the output while recomputing, so we must ignore that transient
  // empty state (and any re-render back to the same text) and wait for the new
  // content to actually land.
  function maybeReady() {
    if (!waiting || !waitFor) return;
    var el = document.getElementById(waitFor);
    var txt = el ? (el.textContent || "").trim() : "";
    if (txt.length > 0 && txt !== waitInitialText) contentUpdated();
  }
  // In "hold" mode the overlay is shown BEFORE the server work starts (e.g. the
  // single-site position change, which is debounced ~800ms). We must not let the
  // idle that follows the mere input flush hide the overlay before that work
  // begins, so we suppress idle-hide until we've seen the recompute's shiny:busy
  // (or HOLD_MS elapses as a fallback for a no-op change).
  var holding = false, sawBusy = false;

  // Per-pathogen palette + spike style for the cartoon. Soft, friendly pastels
  // (COVID is a warm coral rather than an alarming red) so the loading screen
  // never reads as scary.
  var THEME = {
    FLU:   { body: "#F7C066", hi: "#FFE0A6", spike: "#EAA24A", accent: "247,192,102", spikes: 12, tip: "ball" },
    RSV:   { body: "#4CC6B6", hi: "#96E0D6", spike: "#37A99B", accent: "76,198,182",  spikes: 11, tip: "club" },
    COVID: { body: "#F3977E", hi: "#FBC2B1", spike: "#E07C63", accent: "243,151,126", spikes: 12, tip: "club" },
    CHIKV: { body: "#A38FD0", hi: "#CDBFEC", spike: "#8B77BC", accent: "163,143,208", spikes: 18, tip: "thin" }
  };

  // Build a friendly cartoon-virus SVG for a pathogen: a body with a rotating
  // spike corona, glossy highlight, two blinking eyes and a little smile.
  function virusSVG(pid) {
    var t = THEME[pid] || THEME.FLU;
    var cx = 80, cy = 80, inner = 42, outer = (t.tip === "thin") ? 61 : 58;
    var sw = (t.tip === "thin") ? 3 : 4;
    var spikes = "";
    for (var i = 0; i < t.spikes; i++) {
      var a  = (i / t.spikes) * Math.PI * 2;
      var x1 = (cx + Math.cos(a) * inner).toFixed(1), y1 = (cy + Math.sin(a) * inner).toFixed(1);
      var x2 = (cx + Math.cos(a) * outer).toFixed(1), y2 = (cy + Math.sin(a) * outer).toFixed(1);
      spikes += "<line x1=" + x1 + " y1=" + y1 + " x2=" + x2 + " y2=" + y2 +
                " stroke=" + t.spike + " stroke-width=" + sw + " stroke-linecap=round />";
      if (t.tip !== "thin") {
        spikes += "<circle cx=" + x2 + " cy=" + y2 + " r=" + (t.tip === "club" ? 6 : 5) +
                  " fill=" + t.spike + " />";
      }
    }
    return '<svg class="rve-art" viewBox="0 0 160 160" xmlns="http://www.w3.org/2000/svg">' +
             '<g class="rve-spikes">' + spikes + "</g>" +
             "<circle cx=80 cy=80 r=42 fill=" + t.body + " />" +
             "<circle cx=66 cy=63 r=25 fill=" + t.hi + " opacity=0.30 />" +
             '<g class="rve-eyes">' +
               "<circle cx=69 cy=76 r=8 fill=#ffffff />" +
               "<circle cx=91 cy=76 r=8 fill=#ffffff />" +
               "<circle cx=71 cy=78 r=3.8 fill=#24313f />" +
               "<circle cx=93 cy=78 r=3.8 fill=#24313f />" +
               "<circle cx=69.4 cy=75.4 r=1.4 fill=#ffffff />" +
               "<circle cx=91.4 cy=75.4 r=1.4 fill=#ffffff />" +
             "</g>" +
             '<path d="M70 92 Q80 100 90 92" stroke=#24313f stroke-width=3 fill=none stroke-linecap=round />' +
           "</svg>";
  }

  function ensureEl() {
    var el = document.getElementById("rve-switch-overlay");
    if (!el) {
      el = document.createElement("div");
      el.id = "rve-switch-overlay";
      document.body.appendChild(el);
    }
    return el;
  }

  function render(el, msg) {
    var pid = (msg && msg.pathogen) || "FLU";
    var t = THEME[pid] || THEME.FLU;
    el.style.setProperty("--rve-ov-accent", "rgba(" + t.accent + ",0.40)");
    el.innerHTML = virusSVG(pid) +
      '<h3 class="rve-switch-title"></h3>' +
      '<p class="rve-switch-sub"></p>';
    // textContent (not innerHTML) so labels can never inject markup.
    el.querySelector(".rve-switch-title").textContent = (msg && msg.title) || "Loading…";
    el.querySelector(".rve-switch-sub").textContent   = (msg && msg.subtitle) || "";
  }

  function clearTimers() {
    [safetyTimer, idleTimer, minTimer, holdTimer, waitTimer].forEach(function (t) { if (t) clearTimeout(t); });
    safetyTimer = idleTimer = minTimer = holdTimer = waitTimer = null;
  }
  function reallyHide() {
    active = false;
    holding = false;
    sawBusy = false;
    waiting = false;
    waitFor = null;
    stopWaitObserver();
    var el = document.getElementById("rve-switch-overlay");
    if (el) el.classList.remove("rve-visible");
    clearTimers();
  }
  // Arm the "hide after Shiny has been idle for IDLE_MS" timer, unless we're
  // still holding for the recompute's first busy event, or still waiting for a
  // named output to update.
  function armIdleHide() {
    if (!active) return;
    if (holding && !sawBusy) return;
    if (waiting) return;
    if (idleTimer) clearTimeout(idleTimer);
    idleTimer = setTimeout(requestHide, IDLE_MS);
  }
  function requestHide() {
    if (!active) return;
    // Respect a minimum on-screen time so fast switches don't just flicker.
    var elapsed = Date.now() - shownAt;
    if (elapsed >= MIN_MS) { reallyHide(); }
    else { if (minTimer) clearTimeout(minTimer); minTimer = setTimeout(reallyHide, MIN_MS - elapsed); }
  }
  function show(msg) {
    var el = ensureEl();
    render(el, msg);
    clearTimers();
    shownAt = Date.now();
    active = true;
    sawBusy = false;
    holding = !!(msg && msg.hold);
    waitFor = (msg && msg.waitFor) || null;
    stopWaitObserver();
    waiting = !!waitFor;
    el.classList.add("rve-visible");
    safetyTimer = setTimeout(reallyHide, SAFETY_MS);
    if (waiting) {
      var target = document.getElementById(waitFor);
      if (target && window.MutationObserver) {
        waitInitialText = (target.textContent || "").trim();   // the pre-switch text
        waitObserver = new MutationObserver(maybeReady);
        waitObserver.observe(target, { childList: true, subtree: true, characterData: true });
      } else {
        waiting = false;  // can't observe -> don't block the hide
      }
      // Fallback: if the awaited element never changes, stop waiting so the
      // overlay can still idle-hide rather than lingering to SAFETY_MS.
      waitTimer = setTimeout(function () { waiting = false; stopWaitObserver(); armIdleHide(); }, WAIT_MS);
    }
    // Fallback: if the expected recompute never fires a busy (e.g. the user
    // re-picked the same site so nothing changes), stop holding after HOLD_MS so
    // the overlay can still idle-hide instead of lingering until SAFETY_MS.
    if (holding) holdTimer = setTimeout(function () { holding = false; armIdleHide(); }, HOLD_MS);
  }

  // Stay up through the whole busy period (busy cancels a pending hide), then
  // hide once Shiny has been idle for IDLE_MS. No server 'hide' message needed,
  // so it can't get stuck waiting on a flush that never comes.
  $(document).on("shiny:busy", function () {
    sawBusy = true;                                  // the recompute we held for
    if (holdTimer) { clearTimeout(holdTimer); holdTimer = null; }
    if (idleTimer) { clearTimeout(idleTimer); idleTimer = null; }
  });
  $(document).on("shiny:idle", armIdleHide);
  // Backup to the MutationObserver: the awaited output re-rendered.
  $(document).on("shiny:value", function (e) {
    if (!active || !waiting || !waitFor) return;
    var id = (e && (e.name || (e.target && e.target.id)));
    if (id === waitFor) maybeReady();
  });
  // A recompute that errors still reaches shiny:idle, so let the gated idle-hide
  // handle it. Do NOT force-hide here: outputs routinely emit transient
  // shiny:error events (validate() needs, stale-data guards) DURING a switch,
  // and force-hiding on those dropped the mask before the new content rendered.
  $(document).on("shiny:error", function () { if (active) armIdleHide(); });
  $(document).on("shiny:disconnected", reallyHide);

  Shiny.addCustomMessageHandler("rveSwitchOverlayShow", show);
  Shiny.addCustomMessageHandler("rveSwitchOverlayHide", requestHide);
})();
