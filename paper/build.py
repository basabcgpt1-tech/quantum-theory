"""Build paper.pdf from paper.md: pandoc -> HTML with KaTeX math -> Chromium print.
Needs pandoc, node with playwright, and `npm install katex` run once in this directory."""
import re, subprocess, pathlib

here = pathlib.Path(__file__).parent
REFS = {
 "schwinger": "J. Schwinger, <i>Gauge invariance and mass. II</i>, Phys. Rev. <b>128</b>, 2425–2429 (1962).",
 "lowenstein": "J. H. Lowenstein and J. A. Swieca, <i>Quantum electrodynamics in two dimensions</i>, Ann. Phys. <b>68</b>, 172–195 (1971).",
 "cjs": "S. Coleman, R. Jackiw and L. Susskind, <i>Charge shielding and quark confinement in the massive Schwinger model</i>, Ann. Phys. <b>93</b>, 267 (1975).",
 "coleman76": "S. Coleman, <i>More about the massive Schwinger model</i>, Ann. Phys. <b>101</b>, 239–267 (1976).",
 "hamer82": "C. J. Hamer, J. Kogut, D. P. Crewther and M. M. Mazzolini, <i>The massive Schwinger model on a lattice: background field, chiral symmetry and the string tension</i>, Nucl. Phys. B <b>208</b>, 413–438 (1982).",
 "byrnes": "T. M. R. Byrnes, P. Sriganesh, R. J. Bursill and C. J. Hamer, <i>Density matrix renormalisation group approach to the massive Schwinger model</i>, Phys. Rev. D <b>66</b>, 013002 (2002); arXiv:hep-lat/0202014.",
 "adam": "C. Adam, <i>Massive Schwinger model within mass perturbation theory</i>, Ann. Phys. <b>259</b>, 1–63 (1997); arXiv:hep-th/9704064.",
 "act": "E. Arguello Cruz, G. Tarnopolsky and Y. Xin, <i>Precision study of the massive Schwinger model near quantum criticality</i>, arXiv:2412.01902 (2024); Phys. Rev. D, doi:10.1103/vf7m-kzk2.",
 "coleman75": "S. Coleman, <i>Quantum sine-Gordon equation as the massive Thirring model</i>, Phys. Rev. D <b>11</b>, 2088–2097 (1975).",
 "mandelstam": "S. Mandelstam, <i>Soliton operators for the quantized sine-Gordon equation</i>, Phys. Rev. D <b>11</b>, 3026–3030 (1975).",
 "stevenson": "P. M. Stevenson, <i>Gaussian effective potential: Quantum mechanics</i>, Phys. Rev. D <b>30</b>, 1712 (1984); <i>Gaussian effective potential. II. λφ⁴ field theory</i>, Phys. Rev. D <b>32</b>, 1389 (1985).",
 "cardy86": "J. L. Cardy, <i>Effect of boundary conditions on the operator content of two-dimensional conformally invariant theories</i>, Nucl. Phys. B <b>275</b>, 200–218 (1986).",
 "cardy89": "J. L. Cardy, <i>Boundary conditions, fusion rules and the Verlinde formula</i>, Nucl. Phys. B <b>324</b>, 581–596 (1989).",
 "bcn": "H. W. J. Blöte, J. L. Cardy and M. P. Nightingale, <i>Conformal invariance, the central charge, and universal finite-size amplitudes at criticality</i>, Phys. Rev. Lett. <b>56</b>, 742 (1986).",
 "affleck": "I. Affleck, <i>Universal term in the free energy at a critical point and the conformal anomaly</i>, Phys. Rev. Lett. <b>56</b>, 746 (1986).",
 "cc": "P. Calabrese and J. Cardy, <i>Entanglement entropy and quantum field theory</i>, J. Stat. Mech. (2004) P06002; arXiv:hep-th/0405152.",
 "ks": "J. Kogut and L. Susskind, <i>Hamiltonian formulation of Wilson's lattice gauge theories</i>, Phys. Rev. D <b>11</b>, 395 (1975).",
 "bsk": "T. Banks, L. Susskind and J. Kogut, <i>Strong-coupling calculations of lattice gauge theories: (1+1)-dimensional exercises</i>, Phys. Rev. D <b>13</b>, 1043–1053 (1976).",
 "dkpz": "R. Dempsey, I. R. Klebanov, S. S. Pufu and B. Zan, <i>Discrete chiral symmetry and mass shift in the lattice Hamiltonian approach to the Schwinger model</i>, Phys. Rev. Research <b>4</b>, 043133 (2022); arXiv:2206.05308.",
}
src = (here / "paper.md").read_text()
import os
if os.environ.get("PAPER_AUTHOR"):          # set the author line without editing the source
    src = re.sub(r'^author: ".*"$', f'author: "{os.environ["PAPER_AUTHOR"]}"', src, count=1, flags=re.M)
order = []
def cite(m):
    nums = []
    for k in re.findall(r"@(\w+)", m.group(0)):
        if k not in order:
            order.append(k)
        nums.append(order.index(k) + 1)
    return "[" + ", ".join(map(str, nums)) + "]"
body = re.sub(r"\[@[^\]]+\]", cite, src)
missing = [k for k in order if k not in REFS]
assert not missing, missing
lst = "\n".join(f'<li>{REFS[k]}</li>' for k in order)
body = body.replace("::: {#refs}\n:::", f'<ol class="refs">\n{lst}\n</ol>')
(here / "_build.md").write_text(body)
CSS = """
@page { size: A4; margin: 22mm 20mm 22mm 20mm; }
html { font-family: KaTeX_Main, 'Liberation Serif', 'DejaVu Serif', Georgia, serif; font-size: 10.4pt; line-height: 1.45; color: #111; }
body { max-width: none !important; margin: 0 !important; padding: 0 !important; text-align: justify; hyphens: auto; }
h1.title { font-size: 17pt; line-height: 1.25; text-align: center; margin: 0 0 4pt; hyphens: none; }
p.date { text-align: center; color: #444; margin: 0 0 14pt; }
.abstract { margin: 0 9mm 16pt; font-size: 9.6pt; }
h1 { font-size: 13pt; margin: 18pt 0 6pt; text-align: left; break-after: avoid; }
h2 { font-size: 11pt; margin: 13pt 0 4pt; text-align: left; break-after: avoid; }
p { margin: 0 0 6.5pt; orphans: 3; widows: 3; }
.katex { font-size: 1.04em !important; }
.katex-display { margin: 8pt 0; }
span.math.display { display: block; break-inside: avoid; }
table { border-collapse: collapse; margin: 9pt auto; font-size: 8.6pt; text-align: left; break-inside: avoid;
        border-top: 1.2px solid #111; border-bottom: 1.2px solid #111; }
th { border-bottom: 0.8px solid #111; font-weight: 600; }
th, td { padding: 2.2pt 5.5pt; vertical-align: top; }
caption { caption-side: top; font-size: 8.8pt; text-align: left; margin-bottom: 4pt; color: #222; }
caption::before { content: 'Table ' counter(tab) '. '; font-weight: 600; }
table:has(caption) { counter-increment: tab; }
figure { margin: 11pt 0; break-inside: avoid; counter-increment: fig; }
figure img { width: 100%; }
figcaption { font-size: 8.8pt; text-align: left; color: #222; }
figcaption::before { content: 'Figure ' counter(fig) '. '; font-weight: 600; }
body { counter-reset: fig tab; }
code { font-size: 8.8pt; }
ol.refs { font-size: 9pt; padding-left: 16pt; text-align: left; }
ol.refs li { margin-bottom: 2.5pt; }
ul, ol { margin: 0 0 6.5pt; padding-left: 16pt; }
li { margin-bottom: 2pt; }
.header-section-number { margin-right: 0.6em; }
"""
(here / "_style.html").write_text(f"<style>{CSS}</style>")
katex = here / "node_modules" / "katex" / "dist"
assert katex.exists(), "run `npm install katex` in paper/ first"
subprocess.run(["pandoc", "_build.md", "-s", f"--katex=file://{katex}/", "--number-sections", "-H", "_style.html",
                "--metadata", "pagetitle=Two-dimensional QED: mass gap, confinement, criticality",
                "-o", "paper.html"], cwd=here, check=True)
js = r"""
const { chromium } = require('playwright');
(async () => { const b = await chromium.launch(); const p = await b.newPage({ viewport: { width: 643, height: 900 } });
  await p.goto('file://' + process.argv[2]); await p.waitForTimeout(2500);
  const bad = await p.$$eval('.katex-error', e => e.map(x => x.textContent)); if (bad.length) console.log('KATEX ERRORS', bad);
  const raw = await p.$$eval('span.math', e => e.filter(x => !x.querySelector('.katex')).length); if (raw) console.log('UNRENDERED', raw);
  const wide = await p.evaluate(() => { const w = document.documentElement.clientWidth; return [...document.querySelectorAll('.katex-display, table, figure')].filter(e => e.scrollWidth > w + 2 || e.getBoundingClientRect().right > w + 2).map(e => e.textContent.slice(0, 60)); }); if (wide.length) console.log('TOO WIDE', wide);
  await p.pdf({ path: process.argv[3], format: 'A4', printBackground: true, displayHeaderFooter: true,
    headerTemplate: '<span></span>',
    footerTemplate: '<div style="font-size:8px;width:100%;text-align:center;color:#555"><span class="pageNumber"></span></div>',
    margin: { top: '22mm', bottom: '22mm', left: '20mm', right: '20mm' } });
  await b.close(); })();
"""
(here / "_pdf.js").write_text(js)
subprocess.run(["node", "_pdf.js", str(here / "paper.html"), str(here / "paper.pdf")], cwd=here, check=True)
for f in ("_build.md", "_style.html", "_pdf.js"):
    (here / f).unlink()
print("built", here / "paper.pdf")
