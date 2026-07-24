#!/usr/bin/env python3
"""Generate ``result_note.html`` (the paper-evidence log) from ``RESULTS.md``.

``RESULTS.md`` is the single source of truth for our findings; this renders a
visual notebook — one card per result, with the figure(s) embedded inline — so at
paper-writing time the claims, numbers, and figures are all in one place. Re-run
after adding a finding::

    python scripts/build_result_note.py

Parsing contract: ``## R<n> · <title>`` starts a result card; the ``- **Field**:
value`` lines under it become the card's fields (``Figure`` may repeat and is
embedded). Any other ``## <title>`` section (e.g. "Methods notes") is rendered as a
plain bulleted card. Text before the first ``##`` is the lead paragraph.
"""

from __future__ import annotations

import html
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MD_PATH = REPO / "RESULTS.md"
HTML_PATH = REPO / "result_note.html"


def inline_md(text: str) -> str:
    t = html.escape(text)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"\*([^*\n]+)\*", r"<em>\1</em>", t)
    t = re.sub(r"~~([^~]+)~~", r"<del>\1</del>", t)
    t = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', t)
    return t


def parse(md: str):
    lines = md.splitlines()
    title = next((l[2:].strip() for l in lines if l.startswith("# ")), "Results")
    intro: list[str] = []
    sections: list[tuple[str, list[str]]] = []
    cur: list[str] | None = None
    seen_title = False
    for line in lines:
        if line.startswith("# ") and not line.startswith("## "):
            seen_title = True
            continue
        if line.startswith("## "):
            cur = []
            sections.append((line[3:].strip(), cur))
            continue
        if cur is not None:
            cur.append(line)
        elif seen_title and line.strip():
            intro.append(line)
    return title, intro, sections


def _fields(body: list[str]) -> tuple[list[tuple[str, str]], list[str]]:
    """(ordered non-figure fields, figure paths) from ``- **Key**: value`` lines."""
    fields: list[tuple[str, str]] = []
    figures: list[str] = []
    for line in body:
        m = re.match(r"^-\s*\*\*(.+?)\*\*:\s*(.*)$", line.strip())
        if not m:
            continue
        key, val = m.group(1).strip(), m.group(2).strip()
        if key.lower() == "figure":
            figures.append(val)
        else:
            fields.append((key, val))
    return fields, figures


def render_figure(path: str) -> str:
    name = Path(path).name
    png = re.sub(r"\.pdf$", ".png", path)
    if (REPO / png).exists():                    # prefer a committed PNG (renders everywhere)
        media = (f'<a href="{html.escape(path)}"><img src="{html.escape(png)}" '
                 f'alt="{html.escape(name)}" loading="lazy"></a>')
    else:                                        # fall back to embedding the PDF
        media = f'<object data="{html.escape(path)}" type="application/pdf"></object>'
    return (f'<figure class="fig">{media}'
            f'<figcaption><a href="{html.escape(path)}">{html.escape(name)} (PDF)</a></figcaption></figure>')


def render_result(heading: str, body: list[str]) -> str:
    fields, figures = _fields(body)
    fmap = {k.lower(): v for k, v in fields}
    parts = [f'<article class="result"><h2>{inline_md(heading)}</h2>']
    if "claim" in fmap:
        parts.append(f'<p class="claim">{inline_md(fmap["claim"])}</p>')
    if figures:
        parts.append('<div class="figs">' + "".join(render_figure(f) for f in figures) + "</div>")
    dl = [f"<dt>{html.escape(k)}</dt><dd>{inline_md(v)}</dd>"
          for k, v in fields if k.lower() not in ("claim", "for the paper")]
    if dl:
        parts.append("<dl>" + "".join(dl) + "</dl>")
    if "for the paper" in fmap:
        parts.append(f'<div class="for-paper"><span class="tag">for the paper</span> '
                     f'{inline_md(fmap["for the paper"])}</div>')
    parts.append("</article>")
    return "\n".join(parts)


def render_plain(heading: str, body: list[str]) -> str:
    items = []
    for line in body:
        s = line.strip()
        if s.startswith("- "):
            items.append(f"<li>{inline_md(s[2:])}</li>")
    inner = f"<ul>{''.join(items)}</ul>" if items else ""
    return f'<article class="result plain"><h2>{inline_md(heading)}</h2>{inner}</article>'


STYLE = """
:root{--bg:#0f1420;--panel:#171d2b;--panel2:#1e2637;--ink:#e7ecf5;--muted:#95a1b8;
  --line:#2a3346;--accent:#5b8def;--good:#39c07f;}
@media (prefers-color-scheme:light){:root{--bg:#f4f6fb;--panel:#fff;--panel2:#f7f9fd;
  --ink:#1a2233;--muted:#5a6683;--line:#e2e7f0;}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
  font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
a{color:var(--accent);text-decoration:none}a:hover{text-decoration:underline}
code{background:rgba(127,140,170,.18);padding:.05em .4em;border-radius:5px;
  font:12.5px/1.4 ui-monospace,SFMono-Regular,Menlo,monospace;overflow-wrap:anywhere}
.wrap{max-width:1000px;margin:0 auto;padding:34px 22px 80px}
header h1{margin:0 0 8px;font-size:26px;letter-spacing:-.02em}
header .intro{color:var(--muted);max-width:82ch}
header .intro p{margin:.4em 0}
.result{background:var(--panel);border:1px solid var(--line);border-radius:14px;
  padding:20px 22px;margin:20px 0}
.result h2{margin:0 0 10px;font-size:19px}
.claim{font-size:16px;font-weight:600;margin:.2em 0 14px;padding-left:12px;
  border-left:3px solid var(--accent)}
.figs{display:grid;gap:14px;margin:12px 0}
.fig{margin:0;border:1px solid var(--line);border-radius:10px;overflow:hidden;background:var(--panel2)}
.fig object{width:100%;height:440px;display:block;border:0}
.fig img{width:100%;display:block}
.fig figcaption{font-size:11.5px;color:var(--muted);padding:6px 10px;border-top:1px solid var(--line)}
dl{display:grid;grid-template-columns:max-content 1fr;gap:6px 16px;margin:12px 0}
dl dt{color:var(--muted);font-weight:600;font-size:13px}
dl dd{margin:0;font-size:14px}
.for-paper{margin-top:12px;padding:11px 14px;background:rgba(91,141,239,.10);
  border:1px solid rgba(91,141,239,.28);border-radius:10px;font-size:13.5px}
.for-paper .tag{display:inline-block;font-size:10.5px;font-weight:700;text-transform:uppercase;
  letter-spacing:.06em;color:var(--accent);margin-right:8px}
.result.plain ul{margin:.4em 0;padding-left:1.1em}
.result.plain li{margin:.4em 0;font-size:14px}
h2.section{font-size:13px;text-transform:uppercase;letter-spacing:.09em;color:var(--muted);
  margin:34px 0 4px;font-weight:600}
footer{margin-top:40px;color:var(--muted);font-size:12px;border-top:1px solid var(--line);padding-top:16px}
"""


def build() -> str:
    md = MD_PATH.read_text()
    title, intro, sections = parse(md)
    intro_html = "\n".join(f"<p>{inline_md(l)}</p>" for l in intro if l.strip())
    n_results = sum(1 for h, _ in sections if re.match(r"^R\d+\b", h))
    cards = []
    for heading, body in sections:
        cards.append(render_result(heading, body) if re.match(r"^R\d+\b", heading)
                     else render_plain(heading, body))
    return f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title><style>{STYLE}</style>
</head><body><div class="wrap">
<header><h1>{html.escape(title)}</h1><div class="intro">{intro_html}</div>
<p class="section">{n_results} results logged · figures embedded from <code>derived/</code></p></header>
{chr(10).join(cards)}
<footer>Generated from <code>RESULTS.md</code> by <code>scripts/build_result_note.py</code>.
Figures are embedded PDFs — open this file in a browser with the repo checked out. </footer>
</div></body></html>
"""


def main() -> None:
    HTML_PATH.write_text(build())
    _, _, sections = parse(MD_PATH.read_text())
    n = sum(1 for h, _ in sections if re.match(r"^R\d+\b", h))
    print(f"wrote {HTML_PATH.relative_to(REPO)} ({n} result cards)")


if __name__ == "__main__":
    main()
