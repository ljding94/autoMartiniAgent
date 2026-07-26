#!/usr/bin/env python3
"""Generate ``PROGRESS.html`` (a visual dashboard) from ``PROGRESS.md``.

``PROGRESS.md`` stays the single source of truth (the agent-facing changelog);
this script renders a human-friendly HTML view of it — four-ingredient pipeline
cards, a status timeline, and a remaining-work roadmap — so nothing is maintained
twice. Re-run after editing the markdown::

    python scripts/build_progress_html.py

Parsing contract: the ``## Status log`` markdown table (columns
``date | milestone | status``). Rows with a ``YYYY-MM-DD`` date become timeline
entries; rows with ``—`` become roadmap items, grouped by the ``#N`` ingredient
tag in their milestone text. Everything else in the markdown is ignored — keep
the table well-formed and the dashboard follows.
"""

from __future__ import annotations

import html
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MD_PATH = REPO / "PROGRESS.md"
HTML_PATH = REPO / "PROGRESS.html"

INGREDIENTS = [
    ("1", "Stage", "Atomistic reference data", "🧪"),
    ("2", "Process", "Initial AA→CG mapping", "⚙️"),
    ("3", "Evaluation", "Distribution scorer", "📊"),
    ("4", "Loop", "Iterative repair", "🔁"),
]

DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


# ---------- parsing ----------


def parse_rows(md: str) -> list[tuple[str, str, str]]:
    """Extract (date, milestone, status) rows from the status-log table."""
    rows: list[tuple[str, str, str]] = []
    in_table = False
    for line in md.splitlines():
        s = line.strip()
        low = s.lower()
        if s.startswith("|") and "milestone" in low and "status" in low:
            in_table = True
            continue
        if not in_table:
            continue
        if not s.startswith("|"):
            break  # one contiguous table; first non-row line ends it
        if re.match(r"\|[\s:|-]+\|$", s):
            continue  # |---|---|---| separator
        cells = [c.strip() for c in s.strip().strip("|").split("|")]
        if len(cells) < 3:
            continue
        date, status = cells[0], cells[-1]
        milestone = "|".join(cells[1:-1]).strip()  # tolerate stray pipes
        rows.append((date, milestone, status))
    return rows


def status_kind(status: str) -> str:
    """Map a free-form status string to a small set of badge kinds."""
    s = status.lower()
    if "superseded" in s or "deprecated" in s:
        return "superseded"
    if "resolved" in s or s.startswith("done"):   # "done" wins even if a sub-part is deferred
        return "done"
    if "deferred" in s or "out of scope" in s or "out of headline" in s:
        return "deferred"
    if "not started" in s:
        return "todo"
    if "progress" in s or "partial" in s or "wip" in s:
        return "wip"
    if "contingency" in s:
        return "contingency"
    if "open" in s:
        return "open"
    return "info"


KIND_LABEL = {
    "done": "Done", "wip": "In progress", "todo": "Not started",
    "contingency": "Contingency", "open": "Open", "superseded": "Superseded",
    "deferred": "Deferred", "info": "Info",
}


def inline_md(text: str) -> str:
    """Render the inline markdown used in milestone cells (safe: escape first)."""
    t = html.escape(text)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"\*([^*\n]+)\*", r"<em>\1</em>", t)  # single-asterisk italic (after bold)
    t = re.sub(r"~~([^~]+)~~", r"<del>\1</del>", t)
    t = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', t)
    return t


def extract_goal(md: str) -> str:
    m = re.search(r"\*\*Goal\*\*:\s*(.+)", md)
    return inline_md(m.group(1).strip()) if m else ""


def extract_pitch(md: str) -> str:
    m = re.search(r"\*\*One-line pitch\*\*:\s*(.+)", md)
    return inline_md(m.group(1).strip()) if m else ""


def extract_plan(md: str) -> list[tuple[str, bool]]:
    """Ordered steps from the '## Path to publication' numbered list (text, done)."""
    steps: list[tuple[str, bool]] = []
    in_sec = False
    for line in md.splitlines():
        if line.startswith("## Path to publication"):
            in_sec = True
            continue
        if in_sec:
            if line.startswith("## "):
                break
            m = re.match(r"\s*\d+\.\s+(.*)", line)
            if m:
                text = m.group(1).strip()
                done = bool(re.search(r"—\s*done", text, re.I))
                steps.append((text, done))
    return steps


# ---------- rendering ----------


def badge(kind: str, text: str | None = None) -> str:
    return f'<span class="badge {kind}">{html.escape(text or KIND_LABEL.get(kind, kind))}</span>'


def _card_desc(milestone: str) -> str:
    """Short, clean component label for a card: drop the ``**#N Xxx**:`` prefix,
    parenthetical file-name asides, and code markers (full detail lives in the
    timeline; cards are narrow and long code tokens wrap into ugly boxes)."""
    d = re.sub(r"^\*\*#\d[^*]*\*\*:\s*", "", milestone)
    d = re.sub(r"\s*\([^)]*\)", "", d)   # drop "(classify.py, dispatch.py)" asides
    d = d.replace("`", "")
    return re.sub(r"\s+", " ", d).strip(" —-")


def render_cards(roadmap: list[tuple[str, str, str]]) -> str:
    """One card per ingredient, listing its #N-tagged components + a progress bar."""
    by_ing: dict[str, list[tuple[str, str]]] = {n: [] for n, *_ in INGREDIENTS}
    for _date, milestone, status in roadmap:
        m = re.search(r"#(\d)", milestone)
        if not m or m.group(1) not in by_ing:
            continue
        by_ing[m.group(1)].append((_card_desc(milestone), status_kind(status)))

    cards = []
    for num, short, desc, emoji in INGREDIENTS:
        items = by_ing[num]
        active = [(d, k) for d, k in items if k != "deferred"]   # deferred ≠ pending work
        n_def = len(items) - len(active)
        done = sum(1 for _d, k in active if k == "done")
        total = len(active)
        pct = round(100 * done / total) if total else 0
        li = "\n".join(
            f'<li class="{k}"><span class="dot {k}"></span>{inline_md(d)}</li>'
            for d, k in items
        ) or '<li class="info"><span class="dot info"></span>—</li>'
        label = (f"{done}/{total} built" if total else "—")
        if n_def:
            label += f" · {n_def} deferred"
        cards.append(f"""
      <article class="card">
        <div class="card-head">
          <span class="emoji">{emoji}</span>
          <div><span class="ing-num">#{num}</span><h3>{short}</h3><p class="sub">{html.escape(desc)}</p></div>
        </div>
        <div class="progress"><div class="bar" style="width:{pct}%"></div></div>
        <div class="progress-label">{label}</div>
        <ul class="components">{li}</ul>
      </article>""")
    return "\n".join(cards)


def render_timeline(dated: list[tuple[str, str, str]]) -> str:
    entries = []
    for date, milestone, status in dated:
        kind = status_kind(status)
        entries.append(f"""
        <li class="entry status-{kind}" data-status="{kind}">
          <div class="entry-date">{html.escape(date)}</div>
          <div class="entry-body">
            <div class="entry-top">{badge(kind, _status_short(status))}</div>
            <div class="entry-text">{inline_md(milestone)}</div>
          </div>
        </li>""")
    return "\n".join(entries)


def _status_short(status: str) -> str:
    """Trim verbose 'done (2026-.. , ...)' to a compact badge label."""
    k = status_kind(status)
    return KIND_LABEL.get(k, status)


def render_roadmap(roadmap: list[tuple[str, str, str]]) -> str:
    open_items = [(m, s) for _d, m, s in roadmap if status_kind(s) != "done"]
    if not open_items:
        return "<p>Nothing outstanding 🎉</p>"
    lis = "\n".join(
        f'<li class="status-{status_kind(s)}" data-status="{status_kind(s)}">'
        f'{badge(status_kind(s), _status_short(s))} <span>{inline_md(m)}</span></li>'
        for m, s in open_items
    )
    return f'<ul class="roadmap">{lis}</ul>'


STYLE = """
:root{
  --bg:#0f1420; --panel:#171d2b; --panel2:#1e2637; --ink:#e7ecf5; --muted:#95a1b8;
  --line:#2a3346; --accent:#5b8def;
  --done:#39c07f; --wip:#e8b23a; --todo:#6b7688; --open:#4aa8d8; --sup:#5a6274; --cont:#a07bd6;
}
@media (prefers-color-scheme:light){
  :root{--bg:#f4f6fb;--panel:#fff;--panel2:#f7f9fd;--ink:#1a2233;--muted:#5a6683;--line:#e2e7f0;}
}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
  font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;}
a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}
code{background:rgba(127,140,170,.18);padding:.05em .4em;border-radius:5px;
  font:12.5px/1.4 ui-monospace,SFMono-Regular,Menlo,monospace;overflow-wrap:anywhere;}
del{color:var(--muted)}
.wrap{max-width:1080px;margin:0 auto;padding:32px 22px 80px}
header.top{margin-bottom:26px}
header.top h1{margin:0 0 6px;font-size:26px;letter-spacing:-.02em}
header.top .goal{color:var(--muted);max-width:80ch;margin:0}
.meta{display:flex;gap:18px;flex-wrap:wrap;margin-top:14px;color:var(--muted);font-size:13px}
.meta b{color:var(--ink)}
.pitch{margin:14px 0 0;padding:12px 16px;background:var(--panel);border:1px solid var(--line);
  border-left:3px solid var(--accent);border-radius:10px;font-size:14.5px;max-width:88ch}
ol.plan{list-style:none;margin:0;padding:0;display:grid;gap:10px}
ol.plan li{display:grid;grid-template-columns:30px 1fr;gap:14px;align-items:start;
  background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:12px 15px;font-size:13.5px}
ol.plan li .step-num{width:26px;height:26px;border-radius:99px;background:var(--panel2);
  border:1px solid var(--line);color:var(--muted);font-size:13px;font-weight:700;
  display:flex;align-items:center;justify-content:center}
ol.plan li.done .step-num{background:rgba(57,192,127,.15);color:var(--done);border-color:rgba(57,192,127,.4)}
ol.plan li.done{border-left:3px solid var(--done)}
h2.section{font-size:13px;text-transform:uppercase;letter-spacing:.09em;color:var(--muted);
  margin:38px 0 14px;font-weight:600}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:14px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:16px 16px 8px;min-width:0}
.card-head{display:flex;gap:12px;align-items:flex-start;margin-bottom:12px}
.card-head .emoji{font-size:22px;line-height:1}
.card-head h3{margin:2px 0 0;font-size:17px}
.card-head .ing-num{font-size:11px;font-weight:700;color:var(--accent);letter-spacing:.05em}
.card-head .sub{margin:2px 0 0;color:var(--muted);font-size:12.5px}
.progress{height:7px;background:var(--panel2);border-radius:99px;overflow:hidden;border:1px solid var(--line)}
.progress .bar{height:100%;background:linear-gradient(90deg,var(--done),#59d69a)}
.progress-label{font-size:11.5px;color:var(--muted);margin:5px 0 8px}
ul.components{list-style:none;margin:0;padding:0}
ul.components li{display:flex;gap:8px;align-items:baseline;padding:4px 0;font-size:13px;border-top:1px solid var(--line);min-width:0;overflow-wrap:anywhere}
ul.components li:first-child{border-top:none}
ul.components li.done{color:var(--ink)} ul.components li.todo,ul.components li.contingency{color:var(--muted)}
ul.components li.deferred{color:var(--muted);opacity:.7;font-style:italic}
.dot{width:8px;height:8px;border-radius:99px;flex:none;position:relative;top:1px}
.dot.done{background:var(--done)}.dot.wip{background:var(--wip)}.dot.todo{background:var(--todo)}
.dot.open{background:var(--open)}.dot.contingency{background:var(--cont)}.dot.superseded{background:var(--sup)}.dot.info{background:var(--todo)}
.dot.deferred{background:transparent;box-shadow:inset 0 0 0 1.5px var(--sup)}
.badge{display:inline-block;font-size:11px;font-weight:600;padding:2px 9px;border-radius:99px;
  border:1px solid transparent;white-space:nowrap}
.badge.done{background:rgba(57,192,127,.15);color:var(--done);border-color:rgba(57,192,127,.35)}
.badge.wip{background:rgba(232,178,58,.15);color:var(--wip);border-color:rgba(232,178,58,.35)}
.badge.todo{background:rgba(107,118,136,.15);color:var(--muted);border-color:rgba(107,118,136,.3)}
.badge.open{background:rgba(74,168,216,.15);color:var(--open);border-color:rgba(74,168,216,.35)}
.badge.contingency{background:rgba(160,123,214,.15);color:var(--cont);border-color:rgba(160,123,214,.35)}
.badge.superseded{background:rgba(90,98,116,.12);color:var(--sup);border-color:rgba(90,98,116,.3)}
.badge.info{background:rgba(107,118,136,.15);color:var(--muted)}
.badge.deferred{background:rgba(90,98,116,.1);color:var(--sup);border:1px dashed rgba(90,98,116,.45)}
.filters{display:flex;gap:8px;flex-wrap:wrap;margin-bottom:14px}
.filters button{background:var(--panel);color:var(--muted);border:1px solid var(--line);
  border-radius:99px;padding:5px 13px;font-size:12.5px;cursor:pointer}
.filters button.active{color:var(--ink);border-color:var(--accent);background:var(--panel2)}
ul.timeline{list-style:none;margin:0;padding:0;position:relative}
ul.timeline:before{content:"";position:absolute;left:104px;top:6px;bottom:6px;width:2px;background:var(--line)}
li.entry{display:grid;grid-template-columns:96px 1fr;gap:20px;padding:9px 0;position:relative}
li.entry .entry-date{color:var(--muted);font-size:12.5px;text-align:right;padding-top:3px;font-variant-numeric:tabular-nums}
li.entry .entry-body{background:var(--panel);border:1px solid var(--line);border-radius:12px;padding:11px 14px;position:relative}
li.entry .entry-body:before{content:"";position:absolute;left:-20px;top:14px;width:11px;height:11px;
  border-radius:99px;background:var(--todo);border:2px solid var(--bg)}
li.entry.status-done .entry-body:before{background:var(--done)}
li.entry.status-wip .entry-body:before{background:var(--wip)}
li.entry.status-open .entry-body:before{background:var(--open)}
li.entry.status-superseded .entry-body:before{background:var(--sup)}
li.entry.status-superseded{opacity:.62}
.entry-top{margin-bottom:5px}
.entry-text{font-size:13.5px}
ul.roadmap{list-style:none;margin:0;padding:0;display:grid;gap:8px}
ul.roadmap li{display:flex;gap:11px;align-items:baseline;background:var(--panel);
  border:1px solid var(--line);border-radius:10px;padding:9px 13px;font-size:13.5px}
.hidden{display:none!important}
footer{margin-top:44px;color:var(--muted);font-size:12px;border-top:1px solid var(--line);padding-top:16px}
"""

SCRIPT = """
const btns=document.querySelectorAll('.filters button');
btns.forEach(b=>b.onclick=()=>{
  btns.forEach(x=>x.classList.remove('active'));b.classList.add('active');
  const f=b.dataset.filter;
  document.querySelectorAll('[data-status]').forEach(el=>{
    el.classList.toggle('hidden', f!=='all' && el.dataset.status!==f);
  });
});
"""


def build() -> str:
    md = MD_PATH.read_text()
    rows = parse_rows(md)
    dated = [r for r in rows if DATE_RE.match(r[0])]
    roadmap = [r for r in rows if not DATE_RE.match(r[0])]
    dated.sort(key=lambda r: r[0], reverse=True)

    goal = extract_goal(md)
    pitch = extract_pitch(md)
    plan = extract_plan(md)
    last_updated = dated[0][0] if dated else "—"
    n_done = sum(1 for r in dated if status_kind(r[2]) == "done")
    n_open_roadmap = sum(1 for r in roadmap if status_kind(r[2]) != "done")

    filters = ["all", "done", "wip", "open", "todo", "superseded"]
    filter_btns = "\n".join(
        f'<button data-filter="{f}"{" class=\"active\"" if f=="all" else ""}>'
        f'{"All" if f=="all" else KIND_LABEL.get(f, f)}</button>'
        for f in filters
    )

    pitch_html = f'<p class="pitch">{pitch}</p>' if pitch else ""
    plan_html = "\n".join(
        f'<li class="{"done" if done else ""}">'
        f'<span class="step-num">{"✓" if done else i}</span><div>{inline_md(text)}</div></li>'
        for i, (text, done) in enumerate(plan, start=1)
    )
    plan_section = (
        f'<h2 class="section">Path to publication</h2>\n<ol class="plan">{plan_html}</ol>'
        if plan else ""
    )

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>autoMartiniAgent — Progress</title>
<style>{STYLE}</style>
</head>
<body>
<div class="wrap">
  <header class="top">
    <h1>autoMartiniAgent — Progress</h1>
    <p class="goal">{goal}</p>
    {pitch_html}
    <div class="meta">
      <span>Last update <b>{last_updated}</b></span>
      <span><b>{len(dated)}</b> milestones logged</span>
      <span><b>{n_done}</b> done</span>
      <span><b>{n_open_roadmap}</b> items outstanding</span>
    </div>
  </header>

  <h2 class="section">Pipeline — the four ingredients</h2>
  <div class="cards">{render_cards(roadmap)}</div>

  {plan_section}

  <h2 class="section">Status timeline</h2>
  <div class="filters">{filter_btns}</div>
  <ul class="timeline">{render_timeline(dated)}</ul>

  <h2 class="section">Roadmap — outstanding</h2>
  {render_roadmap(roadmap)}

  <footer>
    Generated from <code>PROGRESS.md</code> by <code>scripts/build_progress_html.py</code>
    — edit the markdown, then re-run the script to refresh this view.
  </footer>
</div>
<script>{SCRIPT}</script>
</body>
</html>
"""


def main() -> None:
    HTML_PATH.write_text(build())
    rows = parse_rows(MD_PATH.read_text())
    dated = sum(1 for r in rows if DATE_RE.match(r[0]))
    print(f"wrote {HTML_PATH.relative_to(REPO)} "
          f"({dated} timeline entries, {len(rows) - dated} roadmap rows)")


if __name__ == "__main__":
    main()
