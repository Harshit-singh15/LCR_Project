// LCRBench shared behavior — sequence strip generator + small UI helpers

const AMINO = "ACDEFGHIKLMNPQRSTVWY";

function randSeq(len) {
  let s = "";
  for (let i = 0; i < len; i++) s += AMINO[Math.floor(Math.random() * AMINO.length)];
  return s;
}

// Build a sequence string with a few "low-complexity" runs (repeated letters)
// spliced in, then wrap runs >=4 identical/near-identical chars in <span class="lcr">
function buildMaskedSequence(totalLen) {
  let seq = "";
  while (seq.length < totalLen) {
    if (Math.random() < 0.22) {
      const letter = AMINO[Math.floor(Math.random() * AMINO.length)];
      const runLen = 4 + Math.floor(Math.random() * 6);
      seq += letter.repeat(runLen);
    } else {
      seq += randSeq(4 + Math.floor(Math.random() * 8));
    }
  }
  return seq.slice(0, totalLen);
}

function maskToHtml(seq) {
  return seq.replace(/([ACDEFGHIKLMNPQRSTVWY])\1{3,}/g, (m) => `<span class="lcr">${m}</span>`);
}

function initSeqStrips() {
  document.querySelectorAll(".seq-strip[data-generate]").forEach((el) => {
    const raw = buildMaskedSequence(240);
    const html = maskToHtml(raw + raw); // doubled for seamless scroll
    const track = document.createElement("div");
    track.className = "seq-track";
    track.innerHTML = html.split("").length ? html : "";
    // Actually render as plain markup (spans already inline), wrap words with spacing
    track.innerHTML = html;
    el.innerHTML = "";
    el.appendChild(track);
  });
}

function initNavToggle() {
  const toggle = document.querySelector(".nav-toggle");
  const links = document.querySelector(".nav-links");
  if (!toggle || !links) return;
  toggle.addEventListener("click", () => {
    links.style.display = links.style.display === "flex" ? "none" : "flex";
  });
}

document.addEventListener("DOMContentLoaded", () => {
  initSeqStrips();
  initNavToggle();
});