#!/usr/bin/env python3
"""Build the documentation website from docs/*.md and tests/README.md.

    python3 docs/build_site.py

writes docs/site/index.html and one page per Markdown file, with a sidebar
listing the pages and the sections of the current page. The figures and the
manual are copied into docs/site/, so the directory is the whole website (it is
what .github/workflows/docs.yml publishes on GitHub Pages); links to other files
of the repository point at them on GitHub. The Markdown files stay the source:
re-run the script after editing them. Needs the Python package 'markdown'
(pip install --user markdown).
"""

import html
import os
import re
import shutil
import subprocess
from pathlib import Path

import markdown

DOCS = Path(__file__).resolve().parent
SITE = DOCS / "site"
REPO = DOCS.parent
GITHUB = "https://github.com/thierrydaoud/Rocflu/blob"

# Pages in sidebar order: (Markdown file relative to docs/, generated page,
# short title for the sidebar).
PAGES = [
    ("Building.md", "Building.html", "Building"),
    ("Running.md", "Running.html", "Running a case"),
    ("InputFile.md", "InputFile.html", "Input file"),
    ("PseudoTurb_fixes.md", "PseudoTurb_fixes.html", "PseudoTurb fixes"),
    ("../tests/README.md", "tests.html", "Regression tests"),
]

TOP_ITEM = re.compile(r"^([-*+]|\d+\.) ")
FENCE = re.compile(r"^\s*(```|~~~)")


def indent_of(line):
    return len(line) - len(line.lstrip(" "))


def normalise(text):
    """Adapt GitHub-style Markdown to Python-Markdown.

    Python-Markdown keeps a line inside a list item only under a 4-space indent
    and starts a list only after a blank line; the docs indent list content
    (nested items, tables, paragraphs, code) by 2 and start lists straight after
    a paragraph line. Inside a list, indents are doubled; code inside a fence is
    shifted with its fence so its own indentation is kept.
    """
    lines, in_fence, in_list, loose, shift, previous = [], False, False, False, 0, ""
    for line in text.splitlines():
        if in_fence:
            if line.strip():
                line = " " * shift + line
            if FENCE.match(line):
                in_fence = False
        else:
            if not line.strip():
                loose = in_list
            elif not line.startswith(" "):
                if TOP_ITEM.match(line):
                    # A list item after a paragraph line, or after an item that
                    # holds blank lines, needs a blank line before it.
                    if previous.strip() and (loose or not in_list):
                        lines.append("")
                    loose = loose and in_list
                    in_list = True
                else:
                    in_list = loose = False
            shift = indent_of(line) if in_list and line.startswith(" ") else 0
            line = " " * shift + line
            in_fence = bool(FENCE.match(line))
        lines.append(line)
        previous = line
    return "\n".join(lines) + "\n"


def slugify(value, separator):
    """GitHub's heading anchors, so the docs' own #links keep working."""
    value = html.unescape(value).strip().lower()
    value = re.sub(r"[^\w\- ]", "", value)
    return value.replace(" ", separator)


def branch():
    """The branch the links to repository files point at on GitHub."""
    name = os.environ.get("GITHUB_REF_NAME")
    if not name:
        name = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--abbrev-ref", "HEAD"],
                              capture_output=True, text=True).stdout.strip()
    return name if name and name != "HEAD" else "main"


def asset(file):
    """Copy a file of docs/ (a figure, the manual) into the site; its path there."""
    relative = Path(file).relative_to(DOCS)
    (SITE / relative).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(file, SITE / relative)
    return relative.as_posix()


def rewrite_links(body, source):
    """Point links to the Markdown pages at the generated pages, files of docs/
    at their copy in the site, and other files of the repository at GitHub."""
    pages = {os.path.normpath(DOCS / md): page for md, page, _ in PAGES}
    here = (DOCS / source).parent

    def link(match):
        attr, target = match.group(1), match.group(2)
        path, _, anchor = target.partition("#")
        if re.match(r"^[a-z]+:|^#", target):
            return match.group(0)
        file = os.path.normpath(here / path)
        if file in pages:
            target = pages[file] + ("#" + anchor if anchor else "")
        elif Path(file).is_relative_to(DOCS):
            target = asset(file)
        else:
            target = (f"{GITHUB}/{branch()}/{Path(file).relative_to(REPO).as_posix()}"
                      + ("#" + anchor if anchor else ""))
        return f'{attr}="{target}"'

    return re.sub(r'(href|src)="([^"]*)"', link, body)


def drop_contents(body):
    """Remove each page's hand-written 'Contents' list; the sidebar replaces it."""
    return re.sub(r'<h2 id="contents">Contents</h2>\s*<ul>.*?</ul>\s*(?=<h[12])',
                  "", body, count=1, flags=re.S)


def render(text, source):
    md = markdown.Markdown(
        extensions=["tables", "fenced_code", "toc", "sane_lists"],
        extension_configs={"toc": {"slugify": slugify, "toc_depth": "2-3"}},
    )
    body = md.convert(normalise(text))
    return drop_contents(rewrite_links(body, source)), md.toc_tokens


def sidebar(current, toc):
    items = []
    for _, page, short in PAGES:
        here = ' aria-current="page"' if page == current else ""
        items.append(f'<li><a href="{page}"{here}>{html.escape(short)}</a>')
        if page == current and toc:
            items.append(toc_list(toc))
        items.append("</li>")
    home = ' aria-current="page"' if current == "index.html" else ""
    return (f'<a class="home" href="index.html"{home}>Rocflu</a>'
            f'<ul class="pages">{"".join(items)}</ul>')


def toc_list(tokens):
    # The page title (h1) is the top token; list its sections, skipping 'Contents'.
    sections = [t for top in tokens for t in top["children"]] or tokens
    out = ['<ul class="toc">']
    for t in sections:
        if t["id"] == "contents":
            continue
        out.append(f'<li><a href="#{t["id"]}">{t["name"]}</a>')
        if t["children"]:
            out.append("<ul>")
            out += [f'<li><a href="#{c["id"]}">{c["name"]}</a></li>'
                    for c in t["children"]]
            out.append("</ul>")
        out.append("</li>")
    out.append("</ul>")
    return "".join(out)


def write_page(name, title, nav, body):
    (SITE / name).write_text(TEMPLATE.format(
        title=html.escape(title), nav=nav, body=body))


def main():
    SITE.mkdir(exist_ok=True)
    cards = []
    for md_name, page, short in PAGES:
        text = (DOCS / md_name).read_text()
        body, toc = render(text, md_name)
        title = re.search(r"^# (.+)$", text, re.M).group(1).replace("`", "")
        write_page(page, title, sidebar(page, toc), body)
        intro = re.search(r"<p>(.*?)</p>", body, re.S).group(1)
        intro = intro.replace('href="#', f'href="{page}#')
        cards.append(f'<li><a href="{page}"><strong>{html.escape(short)}</strong></a>'
                     f'<p>{intro}</p></li>')
        print(f"wrote {SITE / page}")

    index = ('<h1>Rocflu documentation</h1>'
             f'<ul class="cards">{"".join(cards)}</ul>'
             '<p>The full RocfluMP manual is a PDF: '
             f'<a href="{asset(DOCS / "Rocflu_manual.pdf")}">Rocflu_manual.pdf</a>.</p>')
    write_page("index.html", "Rocflu documentation",
               sidebar("index.html", None), index)
    print(f"wrote {SITE / 'index.html'}")


TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
:root {{
  --bg: #ffffff; --fg: #1f2328; --muted: #59636e; --line: #d1d9e0;
  --side: #f6f8fa; --code: #eff1f3; --link: #0969da;
}}
@media (prefers-color-scheme: dark) {{
  :root {{
    --bg: #0d1117; --fg: #e6edf3; --muted: #9198a1; --line: #3d444d;
    --side: #151b23; --code: #1f2630; --link: #4493f8;
  }}
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0; background: var(--bg); color: var(--fg);
  font: 16px/1.6 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
  display: grid; grid-template-columns: 18rem minmax(0, 1fr);
}}
nav {{
  position: sticky; top: 0; height: 100vh; overflow-y: auto;
  background: var(--side); border-right: 1px solid var(--line);
  padding: 1.25rem 1rem; font-size: 0.9rem;
}}
nav ul {{ list-style: none; margin: 0; padding: 0; }}
nav .pages > li {{ margin: 0.25rem 0; }}
nav .pages > li > a {{ font-weight: 600; }}
nav .toc {{ margin: 0.25rem 0 0.75rem 0.75rem; border-left: 1px solid var(--line); }}
nav .toc li {{ padding-left: 0.6rem; }}
nav .toc ul li {{ padding-left: 1rem; }}
nav .toc a {{ color: var(--muted); }}
nav a[aria-current] {{ color: var(--fg); }}
nav .home {{ display: block; font-size: 1.1rem; font-weight: 700;
  color: var(--fg); margin-bottom: 1rem; }}
main {{ max-width: 56rem; padding: 2rem 2.5rem 4rem; }}
a {{ color: var(--link); text-decoration: none; }}
a:hover {{ text-decoration: underline; }}
h1, h2, h3 {{ line-height: 1.25; scroll-margin-top: 1rem; }}
h1 {{ font-size: 2rem; border-bottom: 1px solid var(--line); padding-bottom: 0.3em; }}
h2 {{ font-size: 1.5rem; border-bottom: 1px solid var(--line); padding-bottom: 0.3em;
  margin-top: 2.5rem; }}
code {{ font: 0.875em/1.45 ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  background: var(--code); padding: 0.15em 0.35em; border-radius: 4px; }}
pre {{ background: var(--code); padding: 0.9rem 1rem; border-radius: 6px;
  overflow-x: auto; }}
pre code {{ padding: 0; background: none; }}
table {{ border-collapse: collapse; display: block; overflow-x: auto;
  margin: 1rem 0; font-size: 0.92rem; }}
th, td {{ border: 1px solid var(--line); padding: 0.35rem 0.7rem; vertical-align: top; }}
th {{ background: var(--side); text-align: left; }}
img {{ max-width: 100%; height: auto; background: #fff; border-radius: 4px; }}
blockquote {{ margin: 1rem 0; padding: 0 1rem; color: var(--muted);
  border-left: 0.25rem solid var(--line); }}
.cards {{ list-style: none; padding: 0; }}
.cards li {{ border: 1px solid var(--line); border-radius: 6px;
  padding: 0.75rem 1rem; margin: 0.75rem 0; }}
.cards p {{ margin: 0.3rem 0 0; color: var(--muted); }}
@media (max-width: 800px) {{
  body {{ display: block; }}
  nav {{ position: static; height: auto; border-right: 0;
    border-bottom: 1px solid var(--line); }}
  main {{ padding: 1.25rem 1rem 3rem; }}
}}
</style>
</head>
<body>
<nav>{nav}</nav>
<main>
{body}
</main>
</body>
</html>
"""

if __name__ == "__main__":
    main()
