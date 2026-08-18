#!/usr/bin/env python3
"""Generate the lsqvce HTML manuals (user, reference, interpretation) from
authored content, in the toolbox's own CSS. Equations are emitted as LaTeX in
\\(...\\) / \\[...\\] delimiters; STEP 2 (docs_src/render_equations.js) then
pre-renders them to static inline SVG and strips all runtime JS, so the shipped
pages contain no JavaScript and no external references. See docs_src/README.md.
Run from the toolbox root; writes into html/. Reproducible: re-run to rebuild.

Created 2026-07-17 (CEST, UTC+02:00) with assistance from Claude (Anthropic),
model Claude Opus 4.8.
"""
import os

STAMP = ("Generated with assistance from Claude (Anthropic), Claude Opus 4.8, "
         "2026-07-17 (CEST, UTC+02:00).")

# The toolbox page CSS, copied verbatim from the existing html/ pages, plus a
# small equation-numbering helper. MathJax v3 from cdnjs (an allowed CDN in the
# HTML pages; loads when the user has a network, degrades to LaTeX source
# otherwise).
HEAD = '''<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>{title}</title><style>
body {{ font-family: Helvetica, Arial, sans-serif; font-size: 14px;
       color: #000; background: #fff; margin: 0; padding: 0 20px 40px 20px;
       max-width: 900px; line-height: 1.5; }}
h1 {{ font-size: 24px; font-weight: normal; color: #d9541e;
     border-bottom: 1px solid #d9541e; padding-bottom: 6px; margin-top: 24px; }}
h2 {{ font-size: 18px; font-weight: bold; color: #000; margin-top: 28px; }}
h3 {{ font-size: 15px; font-weight: bold; margin-top: 20px; }}
.purpose {{ font-size: 15px; font-style: italic; color: #444;
           margin: 4px 0 20px 0; }}
pre.syntax, pre.example {{
  background: #f3f3f3; border: 1px solid #d3d3d3; border-left: 4px solid #d9541e;
  padding: 10px 12px; font-family: Menlo, Consolas, monospace; font-size: 13px;
  overflow-x: auto; white-space: pre-wrap; }}
pre.example {{ border-left-color: #4a90d9; }}
code {{ font-family: Menlo, Consolas, monospace; background: #f3f3f3;
       padding: 1px 4px; border-radius: 2px; }}
table {{ border-collapse: collapse; margin: 12px 0; width: 100%; }}
th {{ background: #eee; text-align: left; padding: 6px 10px;
     border: 1px solid #ccc; font-size: 13px; }}
td {{ padding: 6px 10px; border: 1px solid #ccc; vertical-align: top;
     font-size: 13px; }}
td.name {{ font-family: Menlo, Consolas, monospace; white-space: nowrap;
          width: 1%; }}
.warn {{ background: #fff4e5; border-left: 4px solid #e8912d; padding: 10px 12px;
        margin: 14px 0; }}
.danger {{ background: #fdecea; border-left: 4px solid #c0392b; padding: 10px 12px;
          margin: 14px 0; }}
.note {{ background: #eaf3fb; border-left: 4px solid #4a90d9; padding: 10px 12px;
        margin: 14px 0; }}
.nav {{ font-size: 13px; color: #666; margin: 14px 0 0 0; }}
.nav a {{ color: #0b5ea8; text-decoration: none; }}
.nav a:hover {{ text-decoration: underline; }}
a {{ color: #0b5ea8; }}
.seealso {{ margin-top: 26px; padding-top: 10px; border-top: 1px solid #ddd;
           font-size: 13px; }}
.footer {{ margin-top: 34px; padding-top: 10px; border-top: 1px solid #ddd;
          font-size: 12px; color: #888; }}
.refs {{ font-size: 12px; color: #333; }}
.refs p {{ margin: 3px 0; text-indent: -22px; padding-left: 22px; }}
.toc {{ font-size: 13px; background: #fafafa; border: 1px solid #e3e3e3;
       padding: 10px 14px; margin: 12px 0; }}
.toc a {{ text-decoration: none; }}
img.eq-inline {{ height: 1.05em; vertical-align: -0.15em; }}
.eq-display {{ text-align: center; margin: 0.9em 0; }}
.eq-display img {{ max-width: 100%; }}
</style>
<!-- Equations are static PNG images (rendered from LaTeX at build time); the
     pages contain no JavaScript and no external references. See docs_src/. -->
</head><body>
<div class="nav"><a href="index.html">lsqvce</a> &rsaquo; {crumb}</div>
'''

FOOT = ('<div class="footer">lsqvce 1.0 &mdash; consolidated diagnostics + '
        'streaming VCE toolbox. {what} ' + STAMP + '</div>\n</body></html>\n')

# Shared reference list (Koch/Teunissen/Baarda/... exactly as in the PDFs).
REFERENCES = r'''<h2 id="refs">References</h2>
<div class="refs">
<p>[1] Baarda, W. (1968). <i>A Testing Procedure for Use in Geodetic Networks.</i> Netherlands Geodetic Commission, Publications on Geodesy, New Series, 2(5). Delft.</p>
<p>[2] Pope, A. J. (1976). <i>The Statistics of Residuals and the Detection of Outliers.</i> NOAA Technical Report NOS 65 NGS 1. Rockville, MD.</p>
<p>[3] Koch, K.-R. (1999). <i>Parameter Estimation and Hypothesis Testing in Linear Models,</i> 2nd ed. Springer, Berlin.</p>
<p>[4] Teunissen, P. J. G. (2006). <i>Testing Theory: an Introduction,</i> 2nd ed. VSSD, Delft.</p>
<p>[5] Teunissen, P. J. G. (2003). <i>Adjustment Theory: an Introduction.</i> VSSD, Delft.</p>
<p>[6] Förstner, W. (1979). Ein Verfahren zur Schätzung von Varianz- und Kovarianzkomponenten. <i>Allgemeine Vermessungs-Nachrichten,</i> 86, 446&ndash;453.</p>
<p>[7] Helmert, F. R. (1907). <i>Die Ausgleichungsrechnung nach der Methode der kleinsten Quadrate,</i> 2nd ed. Teubner, Leipzig.</p>
<p>[8] Grafarend, E. W. &amp; Schaffrin, B. (1993). <i>Ausgleichungsrechnung in linearen Modellen.</i> BI-Wissenschaftsverlag, Mannheim.</p>
<p>[9] Koch, K.-R. &amp; Kusche, J. (2002). Regularization of geopotential determination from satellite data by variance components. <i>Journal of Geodesy,</i> 76(5), 259&ndash;268.</p>
<p>[10] Kusche, J. (2003). A Monte-Carlo technique for weight estimation in satellite geodesy. <i>Journal of Geodesy,</i> 76(11&ndash;12), 641&ndash;652.</p>
<p>[11] Golub, G. H. &amp; Van Loan, C. F. (2013). <i>Matrix Computations,</i> 4th ed. Johns Hopkins University Press, Baltimore.</p>
<p>[12] Björck, Å. (1996). <i>Numerical Methods for Least Squares Problems.</i> SIAM, Philadelphia.</p>
<p>[13] Hutchinson, M. F. (1990). A stochastic estimator of the trace of the influence matrix for Laplacian smoothing splines. <i>Communications in Statistics &ndash; Simulation and Computation,</i> 19(2), 433&ndash;450.</p>
<p>[14] Klees, R., Ditmar, P. &amp; Broersen, P. (2003). How to handle colored observation noise in large least-squares problems. <i>Journal of Geodesy,</i> 76(11&ndash;12), 629&ndash;640.</p>
<p>[15] Ditmar, P., Kuznetsov, V., van Eck van der Sluijs, A. A., Schrama, E. &amp; Klees, R. (2006). A model of the Earth&rsquo;s gravity field. <i>Journal of Geodesy,</i> 79(10&ndash;11), 586&ndash;601.</p>
<p>[16] Durbin, J. (1959). Efficient estimation of parameters in moving-average models. <i>Biometrika,</i> 46(3&ndash;4), 306&ndash;316.</p>
<p>[17] Ljung, G. M. &amp; Box, G. E. P. (1978). On a measure of lack of fit in time series models. <i>Biometrika,</i> 65(2), 297&ndash;303.</p>
<p>[18] Wen, H. Y., Kruizinga, G., Paik, M., et al. (2019). <i>GRACE-FO Level-1 Data Product User Handbook.</i> JPL D-56935, Jet Propulsion Laboratory, Pasadena.</p>
<p>[19] Kornfeld, R. P., Arnold, B. W., Gross, M. A., et al. (2019). GRACE-FO: the Gravity Recovery and Climate Experiment Follow-On mission. <i>Journal of Spacecraft and Rockets,</i> 56(3), 931&ndash;951.</p>
<p>[20] Schaffrin, B. (1997). Reliability measures for correlated observations. <i>Journal of Surveying Engineering,</i> 123(3), 126&ndash;137.</p>
<p>[21] Ellmer, M. (2018). <i>Contributions to GRACE Gravity Field Recovery: Improvements in Dynamic Orbit Integration, Stochastic Modelling of the Antenna Offset Correction, and Co-Estimation of Satellite Orbits.</i> PhD thesis, Graz University of Technology.</p>
</div>'''

import re as _re, os as _os, subprocess as _sp, hashlib as _hl
try:
    from PIL import Image as _Image
except ImportError:
    _Image = None

_EQ_DIR = 'html/eq'
_EQ_PRE = (r'\documentclass[border=1pt,varwidth=6in]{standalone}'
           r'\usepackage{amsmath,amssymb,amsfonts,bm}\begin{document}%s\end{document}')

def _render_eq(tex, kind):
    """Render one LaTeX equation to a cropped transparent PNG; return its href.
    Build-time only: needs pdflatex + pdftoppm + Pillow. Images are cached by a
    content hash so identical equations share one file and rebuilds are cheap."""
    _os.makedirs(_EQ_DIR, exist_ok=True)
    h = _hl.md5((kind + tex).encode()).hexdigest()[:12]
    out = f'{_EQ_DIR}/{h}.png'
    if _os.path.exists(out):
        return out
    body = _re.sub(r'\\tag\{[^}]*\}', '', tex).strip()
    tag = _re.search(r'\\tag\{([^}]*)\}', tex)
    if kind == 'display':
        math = r'$\displaystyle %s$' % body
        if tag:
            math = r'\hbox{%s\qquad(%s)}' % (math, tag.group(1))
        dpi = 260
    else:
        math = r'$%s$' % body
        dpi = 240
    d = _os.path.abspath('html/eq/_tmp')
    _os.makedirs(d, exist_ok=True)
    open(d + '/e.tex', 'w').write(_EQ_PRE % math)
    _sp.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', 'e.tex'],
            cwd=d, capture_output=True, text=True)
    if not _os.path.exists(d + '/e.pdf'):
        raise RuntimeError('pdflatex failed for equation: ' + tex[:60])
    _sp.run(['pdftoppm', '-png', '-r', str(dpi), '-transp', 'e.pdf', 'pg'],
            cwd=d, capture_output=True)
    png = d + '/pg-1.png'
    if _Image is not None and _os.path.exists(png):
        im = _Image.open(png).convert('RGBA')
        bb = im.getbbox()
        if bb:
            im = im.crop(bb)
        im.save(out)
    else:
        _os.replace(png, out)
    return out

def _sub_equations(body):
    """Replace \\[...\\] display and \\(...\\) inline LaTeX with <img> tags."""
    def disp(m):
        href = _render_eq(m.group(1).strip(), 'display')
        return ('<div class="eq-display"><img src="eq/%s" alt="equation"></div>'
                % _os.path.basename(href))
    def inl(m):
        href = _render_eq(m.group(1).strip(), 'inline')
        return '<img class="eq-inline" src="eq/%s" alt="equation">' % _os.path.basename(href)
    body = _re.sub(r'\\\[(.+?)\\\]', disp, body, flags=_re.S)
    body = _re.sub(r'\\\((.+?)\\\)', inl, body, flags=_re.S)
    return body

def write(path, title, crumb, body, what):
    body = _sub_equations(body)
    html = HEAD.format(title=title, crumb=crumb) + body + FOOT.format(what=what)
    open(path, 'w').write(html)
    print(f"wrote {path}  ({len(html)} bytes)")


# ======================================================================
# 1. INTERPRETATION GUIDE  (prose, no equations)
# ======================================================================
INTERP = r'''<h1>Interpretation Guide &mdash; reading the result</h1>
<p class="purpose">What to run, and how to read what comes back. This guide is
about judgement, not syntax: which number to look at first, what a red flag
means, and the two ways least squares will quietly mislead you. The equations
are in the <a href="lsqvce_user_manual.html">User Manual</a>, Part&nbsp;I.</p>

<h2>The order to read things</h2>
<p>Read the report top-down and stop at the first failure &mdash; a break at one
level makes the levels below it meaningless.</p>
<table>
<tr><th>Level</th><th>Reads</th><th>If red</th><th>Why it comes first</th></tr>
<tr><td>0 &mdash; Integrity</td><td>the trace identity</td><td>a data/weighting bug</td>
    <td>if &sum;r&#8341; &ne; n&minus;u the geometry is wrong; nothing else is trustworthy</td></tr>
<tr><td>1 &mdash; Global</td><td>model consistency</td><td>inconsistent</td>
    <td>s&#8320;&#178;/&sigma;&#8320;&#178; too large means blunders or optimistic weights</td></tr>
<tr><td>2 &mdash; Localise</td><td>find the culprit</td><td>snoop</td>
    <td>&mdash; but check contamination first</td></tr>
<tr><td>3 &mdash; Reliability</td><td>what you cannot check</td><td>you are blind here</td>
    <td>low r&#8341; = unverifiable, regardless of residual</td></tr>
<tr><td>4 &mdash; Parameters</td><td>result fragility</td><td>results fragile</td>
    <td>high correlation / conditioning = weakly separated unknowns</td></tr>
</table>

<h2>Two ways least squares lies to you</h2>
<div class="danger"><b>Lie 1 &mdash; &ldquo;the residual is small, so the observation is good.&rdquo;</b>
False. A residual shows only the fraction \(r_i\) of the error present,
\(E\{\Delta v_i\} = -r_i\,b_i\). At \(r_i \approx 0\) even a huge blunder leaves a
near-zero residual. Always pair a residual with its redundancy number. A small
residual at low \(r_i\) is the most dangerous reading in the whole report, not
the most reassuring.</div>
<div class="danger"><b>Lie 2 &mdash; &ldquo;the tau test found no outliers, so the data are clean.&rdquo;</b>
Pope&rsquo;s tau normalises by \(s_0\), which the blunders themselves inflate. Past a
few per-cent contamination the test goes blind &mdash; a 100&sigma; error can pass.
Read the contamination gauge before you trust any tau, and switch to a robust
scale (<code>'test','mad'</code>) when it is high.</div>

<h2>The quadrant you must inspect</h2>
<table>
<tr><th></th><th>low |w| (quiet)</th><th>high |w| (loud)</th></tr>
<tr><td><b>high r (checked)</b></td><td>healthy</td><td>plain outlier &mdash; snoop it</td></tr>
<tr><td><b>low r (unchecked)</b></td><td>blind spot &mdash; do not relax</td>
    <td>leverage point &mdash; small w still hides a big error</td></tr>
</table>
<p>The bottom row is where reputations are lost. If a parameter you care about
depends on bottom-row observations, no amount of clean residuals should reassure
you; strengthen the geometry or add independent data.</p>

<h2>Deciding what to do next</h2>
<h3>Level 1 failed &mdash; is it blunders or the model?</h3>
<p>Run <code>lsqvce.systematic(D)</code> before snooping. If the residuals show
runs, serial correlation, or per-group offsets, the problem is the model (an
unmodelled signal, a wrong weight ratio, coloured noise treated as white) &mdash;
snooping it away just deletes good data around a defect that remains. Fix the
model, then re-solve.</p>
<h3>Diffuse contamination &mdash; reject or reweight?</h3>
<p>Isolated blunders: iterative <code>lsqvce.snoop</code>. Many small
contaminations: robust M-estimation (<code>lsqvce.robust</code>) degrades
gracefully where snooping would chase its own tail. The contamination gauge tells
you which regime you are in.</p>
<h3>Group weights look wrong &mdash; variance components</h3>
<p>If one data type or satellite is systematically over-fit or under-fit, run
<code>lsqvce.varComp</code> / <code>lsqvce.vce</code> and read <code>vc.sigma</code>:
a factor far from 1 means that group was mis-weighted going in. Re-solve at the
converged weights before reading reliability &mdash; MDBs computed at wrong weights
are wrong. Confirm <code>vc.converged</code> is true and that the minimum
redundancy was met for every group. For \(u\) in the thousands, the matrix-free
<a href="lsqvce.vceHutchinson.html">Hutchinson backend</a> gives the same sigmas;
see <a href="lsqvce_choosing.html">Choosing an engine</a>.</p>

<h2>Reading a gravity solution</h2>
<p>For the block-arrow path (<a href="lsqvce.adjustBlock.html">lsqvce.adjustBlock</a> /
<a href="lsqvce.adjustVCE.html">lsqvce.adjustVCE</a>), three habits keep you honest:</p>
<ul>
<li><b>Diagnose through the accumulator bridge</b> (<a href="lsqvce.diagnoseAccum.html">lsqvce.diagnoseAccum</a>),
never the reduced common system alone. The common-parameter redundancies are only
correct once the arc-local parameters have been counted; skip the bridge and every
reliability number is inflated.</li>
<li><b>Map the redundancy of the field, not just its residuals.</b> A region carried
by a single arc is a blind spot; its formal error can look fine while being
uncheckable.</li>
<li><b>Watch for a near-singular arc.</b> An arc with too few observations to fix its
own state poisons the Schur complement; the accumulator flags it &mdash; do not force
it in.</li>
</ul>
<div class="note">A clean bill of health means: Level&nbsp;0 exact, Level&nbsp;1
consistent, no high-w observations left, contamination low, and the parameters you
care about backed by high-redundancy observations with acceptable MDB. Any one of
those missing turns &ldquo;passed&rdquo; into &ldquo;passed, but&hellip;&rdquo;.</div>

<div class="seealso"><b>See also:</b>
<a href="lsqvce_user_manual.html">User Manual</a>,
<a href="lsqvce_reference_manual.html">Reference Manual</a>,
<a href="lsqvce_choosing.html">Choosing an engine</a>,
<a href="lsqvce.diagnose.html">lsqvce.diagnose</a>,
<a href="lsqvce.report.html">lsqvce.report</a>,
<a href="lsqvce.systematic.html">lsqvce.systematic</a></div>
''' + REFERENCES


# ======================================================================
# 2. REFERENCE MANUAL  (symbols, equations, conventions)
# ======================================================================
REFMAN = r'''<h1>Reference Manual &mdash; symbols, equations, conventions</h1>
<p class="purpose">Companion to the <a href="lsqvce_user_manual.html">User Manual</a>.
Equation numbers here are self-contained. Everything is streaming and
weight-based (\(p = 1/\sigma^2\), never variance).</p>

<h2>Symbols</h2>
<table>
<tr><th>Symbol</th><th>Meaning</th></tr>
<tr><td class="name">\(n,\ u,\ r\)</td><td>observations, parameters, redundancy \(r = n-u\)</td></tr>
<tr><td class="name">\(A,\ l,\ e\)</td><td>design (Jacobian), observations, errors</td></tr>
<tr><td class="name">\(C,\ P\)</td><td>observation cofactor and weight, \(P = C^{-1}\)</td></tr>
<tr><td class="name">\(\sigma_0^2,\ s_0^2\)</td><td>a-priori and a-posteriori variance of unit weight</td></tr>
<tr><td class="name">\(N,\ b\)</td><td>normal matrix \(A^{\!\top}PA\) and right side \(A^{\!\top}Pl\)</td></tr>
<tr><td class="name">\(Q_{xx}=N^{-1}\)</td><td>parameter cofactor</td></tr>
<tr><td class="name">\(v,\ Q_{vv}\)</td><td>residuals, residual cofactor</td></tr>
<tr><td class="name">\(r_i,\ h_i\)</td><td>redundancy number and leverage, \(r_i = 1 - h_i\) (whitened)</td></tr>
<tr><td class="name">\(w_i,\ \tau_i\)</td><td>Baarda standardised residual, Pope tau</td></tr>
<tr><td class="name">\(n_c,\ n_l,\ K\)</td><td>common params, local params per arc, number of arcs</td></tr>
</table>

<h2>1. Estimator and dispersion</h2>
<p>The Gauss&ndash;Markov model and its best linear unbiased estimate:</p>
\[ l = Ax + e,\qquad E\{e\}=0,\qquad D\{e\}=\sigma_0^2\,C=\sigma_0^2 P^{-1}, \]
\[ \hat{x} = N^{-1}b,\quad N=A^{\!\top}PA,\quad b=A^{\!\top}Pl,\qquad
   Q_{\hat x\hat x}=N^{-1},\quad D\{\hat x\}=\sigma_0^2 N^{-1}. \tag{1} \]
<p>Solved by whitening \(A_w = RA\) (with \(P = R^{\!\top}R\)) and thin-QR
\(A_w = Q_1 R\), so the condition number is not squared [11,12]. The accumulator
updates \(R\) (or \(N,b\)) batch by batch.</p>

<h2>2. Residuals, redundancy, reliability</h2>
\[ v = A\hat{x} - l,\qquad s_0^2 = \frac{v^{\!\top}Pv}{n-u}. \tag{2} \]
<p>With the hat (influence) matrix \(H = A_w N^{-1} A_w^{\!\top}\), in whitened
coordinates the leverage and partial redundancy of observation \(i\) are</p>
\[ h_i = a_{w,i}^{\!\top} N^{-1} a_{w,i},\qquad r_i = 1 - h_i,\qquad
   \sum_{i=1}^{n} r_i = n - u. \tag{3} \]
<p>The last identity is what Level&nbsp;0 of the report checks exactly. The
response of a residual to a single gross error \(b_i\) in observation \(i\) is
\(E\{\Delta v_i\} = -r_i\,b_i\). Baarda&rsquo;s standardised residual and, when
\(\sigma_0\) is unknown, Pope&rsquo;s tau are</p>
\[ w_i = \frac{-v_i\sqrt{p_i}}{\sigma_0\sqrt{r_i}},\qquad
   \tau_i = \frac{-v_i\sqrt{p_i}}{s_0\sqrt{r_i}}. \tag{4} \]
<p>Internal reliability &mdash; the minimal detectable bias at test size
\(\alpha_0\), power \(\beta_0\) &mdash; and its propagation into the parameters
(external reliability) are</p>
\[ \nabla_{0}l_i = \frac{\delta_0\,\sigma_0}{\sqrt{p_i\,r_i}},\qquad
   \nabla \hat{x} = N^{-1} a_{w,i}\,\nabla_0 l_i. \tag{5} \]
<p>The MDB scales as \(1/\sqrt{r_i}\): at \(r_i\to 0\) it diverges &mdash; nothing
can be checked there. References: redundancy and reliability [1,4,20]; tau [2];
leverage by back-substitution through the stored triangular factor, not the
Hutchinson trace estimator [13] (exact partial redundancies are cheap when \(u\)
is small, and exactness matters for these bounds).</p>

<h2>3. Variance components</h2>
<p>For observation group \(g\), Helmert&rsquo;s variance factor is the group&rsquo;s
weighted residual sum of squares over its partial redundancy [6,7,3]:</p>
\[ \hat{\sigma}_g^2 = \frac{v_g^{\!\top}P_g v_g}{r_g},\qquad
   r_g = n_g - \operatorname{tr}(N^{-1}N_g),\quad N_g = A_g^{\!\top}P_g A_g. \tag{6} \]
<p>Förstner&rsquo;s simplification writes \(r_g = \sum_{i\in g} r_i\). Iterate
\(P_g \leftarrow \hat{\sigma}_g^{-2}P_g\) to unity; only ratios are identifiable,
so one group (or the pooled \(s_0\)) fixes the gauge. VCE for gravity and
regularisation: [9,10]. For \(u\) in the thousands the trace
\(\operatorname{tr}(N^{-1}N_g)\) is estimated matrix-free with \(Z\) Hutchinson
probes (<a href="lsqvce.vceHutchinson.html">lsqvce.vceHutchinson</a>); see
<a href="lsqvce_choosing.html">Choosing an engine</a>.</p>

<h3>3a. Relation to Ellmer (2018)</h3>
<p>The matrix-free trace machinery here is Ellmer&rsquo;s [21, &sect;6.5.6]: the
Hutchinson estimator [13] with Rademacher probes \(z\in\{-1,+1\}^{P}\) and the
&ldquo;stacked&rdquo; form \(\operatorname{tr}(X)\approx
\operatorname{tr}(\bar{Z}^{\!\top}X\bar{Z})\), \(\bar{Z}=Z/\sqrt{Z}\), applied
through the shared Cholesky factor so \(N\) is never inverted. The between-group
update \(\hat{\sigma}_g^2=\hat{\alpha}_g^2\sigma_g^2\) with
\(\hat{\alpha}_g^2 = v_g^{\!\top}P_g v_g / (n_g-\operatorname{tr}(N_gN^{-1}))\) is
his eq. (2.6.12). What this toolbox does <em>differently</em> is worth stating
plainly, because it changes what the estimator is pointed at:</p>
<ul>
<li><b>Different operator.</b> Ellmer applies Hutchinson to the redundancy of the
<em>temporal covariance components</em>, \(s_n=\operatorname{tr}(RV_n)\) over the
lag matrices \(V_n\) of one arc&rsquo;s noise model [21, eqs. (6.5.36)&ndash;(6.5.44)]
&mdash; the SHAPE side. Here the same estimator is applied to the
<em>between-group reduced common normal</em>,
\(\operatorname{tr}(N^{-1}N_g^{cc})\) &mdash; the SCALE side and the contribution
matrix \(C_g=N^{-1}N_g^{cc}\). Same tool, different target.</li>
<li><b>Streaming block-arrow form.</b> Ellmer factors the full normal matrix and
treats \(AN^{-1}A^{\!\top}\) monolithically. Here \(N_g^{cc}z\) is formed by
<em>streaming the arcs of group \(g\)</em> after per-system local Schur
elimination, never materialising any group normal (see &sect;4), so the memory
cost is \(O(\text{probe block})\), not \(O(G\,u_c^2)\).</li>
<li><b>Symmetric group split (new).</b> When one system&rsquo;s arcs span several
VCE groups (e.g. grouping a merged pair-system by observation type), the
system&rsquo;s Schur term cannot be assigned to a single group. We use the
ordering-independent split
\[ N_g^{cc}=N_{cc,g}-\tfrac12\!\left(N_{lc,g}^{\!\top}N_{ll,s}^{-1}N_{lc,s}
   +N_{lc,s}^{\!\top}N_{ll,s}^{-1}N_{lc,g}\right), \]
which partitions the reduced normal exactly (\(\sum_g N_g^{cc}=N\), so the
group traces sum to \(u_c\) and the contribution diagonals to \(\mathbf{1}\)).
Ellmer&rsquo;s groups are uncorrelated (block-diagonal \(\Sigma_{ll}\)), so this
coupling does not arise in the thesis and the split is an addition here.</li>
<li><b>Exact by default; diagonal as well as trace.</b> Ellmer uses Hutchinson as
the production path because forming \(R\) explicitly is prohibitive. Here the
<em>exact</em> per-arc leverage path is the default and the authority, with
Hutchinson opt-in for large \(u_c\); and the same probes also yield the
per-parameter contribution diagonal \(\operatorname{diag}(N^{-1}N_g^{cc})\),
which is beyond the scalar trace the VCE update needs
(<a href="lsqvce.streamContribution.html">lsqvce.streamContribution</a>,
<a href="lsqvce.contribution.html">lsqvce.contribution</a>).</li>
</ul>

<h2>4. Block-arrow system and QR pre-elimination</h2>
<p>Gravity recovery stacks \(K\) arcs (or satellite-days) that share the \(u_c\)
common parameters (the field) but each carry their own \(u_k\) local nuisance
parameters. Ordering locals first, then commons, the normal matrix is bordered
block-diagonal &mdash; a block arrow. Eliminating every local block by its Schur
complement leaves the reduced common system</p>
\[ \Big(N_{cc} - \sum_{k} N_{ck}N_{kk}^{-1}N_{kc}\Big)\hat{x}_c
   = b_c - \sum_{k} N_{ck}N_{kk}^{-1}b_k. \tag{7} \]
<p>Per-arc QR (the numerically preferred path): one Householder QR triangularises
the \(u_k\) local columns and applies the same reflectors to the common columns
and the whitened observations,</p>
\[ Q_k^{\!\top}\big[\,A_k^{w}\ \ A_c^{w}\ \big|\ l_k^{w}\,\big]
   = \begin{bmatrix} R_{kk} & R_{kc} & z_k\\ 0 & S_k & y_k \end{bmatrix}. \tag{8} \]
<p>The trailing rows \((S_k,y_k)\) are the arc&rsquo;s contribution to the common
problem, cleaned of its nuisance parameters &mdash; a QR realisation of the Schur
complement, formed without squaring the design. The accumulator appends each
\((S_k,y_k)\) to a running QR of the reduced common system, so memory is
\(O(u_c^2)\) independent of \(K\). Each local block follows by back-substitution,</p>
\[ \hat{x}_k = R_{kk}^{-1}\big(z_k - R_{kc}\hat{x}_c\big). \tag{9} \]
<div class="warn"><b>Why the locals must be counted, not discarded.</b>
Pre-elimination removes the locals from the solve, but they still spend degrees
of freedom: each arc&rsquo;s locals absorb part of its own residual, lowering the
redundancy numbers of that arc&rsquo;s observations. With <code>'KeepLocal',true</code>
the per-arc factors are retained so
<a href="lsqvce.diagnoseAccum.html">lsqvce.diagnoseAccum</a> rebuilds the full hat
matrix and reports the true \(r_i\) (hence MDB and reliability) rather than the
inflated commons-only values [3,10]. Integrity check:
\(\sum_i r_i = n - (u_c + \sum_k u_k)\).</div>

<h2>5. Whitening coloured noise</h2>
<p>If the noise is a stationary AR(\(p\)) process with polynomial
\(A(q^{-1}) = 1 + a_1 q^{-1} + \dots + a_p q^{-p}\), the FIR recursion</p>
\[ \tilde{e}_t = e_t + a_1 e_{t-1} + \dots + a_p e_{t-p} \tag{10} \]
<p>whitens it; applying the same filter to every column of \(A\) and to \(l\)
realises \(P = C^{-1}\) implicitly &mdash; a banded operator, never the dense
\(n\times n\) matrix. For vector observations a VAR(\(p\)) whitener generalises it
with the exact stationary start-up for the first \(p\) epochs, so the operator
equals the dense block-Cholesky GLS solve row for row. When the covariance
separates in space and time, \(C = C_t \otimes C_s\), the whitener is the
Kronecker product</p>
\[ R = R_t \otimes R_s, \qquad C = R^{-1}R^{-\top}, \tag{11} \]
<p>(<a href="lsqvce.separableWhiten.html">lsqvce.separableWhiten</a>) &mdash; what
makes per-epoch vector data affordable. Orders are chosen by information criteria
and validated by a Ljung&ndash;Box test on the whitened residuals [17].</p>
<p><a href="lsqvce.selectArAcv.html">lsqvce.selectArAcv</a> chooses the AR order
without a candidate grid: one Levinson recursion yields the innovation variance of
every lower order through the reflection coefficients \(\kappa_j\),
\(\sigma_p^2 = r_0\prod_j(1-\kappa_j^2)\), and it returns the smallest order whose
refitted autocovariance is within <code>AcvTol</code> &mdash; the same decorrelation
target the whiteness test uses. It guards two failure modes a single AR cannot fix:
a near-unit root (detrend first) and a plateau above tolerance (sharp 1/rev, 2/rev
lines &mdash; move to <a href="lsqvce.ar2arma.html">lsqvce.ar2arma</a>).</p>
<div class="note"><b>Separable pipeline (two layers).</b> The model is
\(C_g = \sigma_g^2\,D_g\,(R_{t,g}\otimes I_d)\,D_g\). Whitening estimates the
temporal SHAPE \(R_{t,g}\) <i>within</i> each group; variance-component estimation
(&sect;3) estimates the SCALE \(\sigma_g^2\) <i>between</i> groups. The loop
iterates solve &rarr; whiten &rarr; VCE-reweight and converges only when every
stream is white AND every \(\sigma_g^2 \approx 1\). Driven end to end by
<a href="lsqvce.adjustVCE.html">lsqvce.adjustVCE</a>.</div>

<div class="seealso"><b>See also:</b>
<a href="lsqvce_user_manual.html">User Manual</a>,
<a href="lsqvce_choosing.html">Choosing an engine</a>,
<a href="lsqvce.redundancy.html">lsqvce.redundancy</a>,
<a href="lsqvce.vce.html">lsqvce.vce</a>,
<a href="lsqvce.adjustBlock.html">lsqvce.adjustBlock</a>,
<a href="lsqvce.separableWhiten.html">lsqvce.separableWhiten</a></div>
''' + REFERENCES


# ======================================================================
# 3. USER MANUAL  (theory + usage + step-by-step; NEW names)
# ======================================================================
USERMAN = r'''<h1>User Manual &mdash; judging an adjustment, getting the stochastic model right</h1>
<p class="purpose">The full walk-through: the theory each diagnostic implements,
what you supply, how to read the output, the gravity-recovery path, and a
step-by-step recipe. New to the toolbox? The runnable version of the recipe is
<code>examples/demo_workflow.m</code>; to pick an engine see
<a href="lsqvce_choosing.html">Choosing an engine</a>.</p>

<div class="toc"><b>Contents.</b>
<a href="#p1">I. Theory</a> &middot;
<a href="#p2">II. What you supply</a> &middot;
<a href="#p3">III. Reading the diagnosis</a> &middot;
<a href="#p4">IV. Gravity recovery</a> &middot;
<a href="#p5">V. Step-by-step</a> &middot;
<a href="#refs">References</a></div>

<div class="note"><b>The single most important sentence.</b> A residual reveals
only the fraction \(r_i\) of any error present in observation \(i\):
\(E\{\Delta v_i\} = -r_i\,b_i\). At \(r_i \to 0\) a 40&sigma; blunder produces
\(|w|=0.000\). This is not a defect of the estimator; it is the geometry of least
squares, and it is why this toolbox reports reliability, not just residuals.</div>

<h2 id="p1">Part I &mdash; Theory</h2>
<p>Notation follows Koch [3] and Teunissen [5]. Vectors are columns; \(\top\) is
transpose; \(n\) observations, \(u\) parameters, \(r=n-u\) redundancy. The
compact equation set is in the <a href="lsqvce_reference_manual.html">Reference
Manual</a>; the essentials:</p>

<h3>Gauss&ndash;Markov model and BLUE</h3>
\[ l = Ax+e,\quad D\{e\}=\sigma_0^2 P^{-1};\qquad
   \hat{x}=N^{-1}b,\ N=A^{\!\top}PA,\ b=A^{\!\top}Pl,\ D\{\hat x\}=\sigma_0^2 N^{-1}. \]
<p>The streaming accumulator forms \(N,b\) by summing batch contributions
\(A^{\!\top}PA\) and \(A^{\!\top}Pl\), so cost is \(O(u^2)\) memory and
\(O(u^2 n)\) time, independent of \(n\).</p>

<h3>QR versus the normal equations</h3>
<p>Forming \(N\) squares the condition number, \(\operatorname{cond}(N)=
\operatorname{cond}(A_w)^2\). The remedy is to whiten (\(P=R^{\!\top}R\)) and solve
by orthogonal triangularisation \(A_w=Q_1R\). The <code>'qr'</code> backend
accumulates \(R\) batch by batch, keeping QR accuracy at normal-equation memory
[11,12].</p>

<h3>Residuals, redundancy, reliability</h3>
\[ s_0^2=\frac{v^{\!\top}Pv}{n-u},\qquad h_i=a_{w,i}^{\!\top}N^{-1}a_{w,i},\qquad
   r_i=1-h_i,\qquad \sum_i r_i=n-u. \]
<p>The trace identity is what Level&nbsp;0 checks. Baarda&rsquo;s w-test
\(w_i=-v_i\sqrt{p_i}/(\sigma_0\sqrt{r_i})\) is the locally most powerful test for a
single gross error; with \(\sigma_0\) unknown it becomes Pope&rsquo;s tau [1,2].
The minimal detectable bias \(\nabla_0 l_i=\delta_0\sigma_0/\sqrt{p_i r_i}\) scales
as \(1/\sqrt{r_i}\): at \(r_i\to0\) nothing can be checked, and the undetected bias
propagates into the parameters as \(\nabla\hat x=N^{-1}a_{w,i}\nabla_0 l_i\)
[1,4,20]. Leverage is computed by back-substitution through the stored triangular
factor, never the explicit inverse; the Hutchinson trace estimator [13] is used
only for the large-\(u\) variance-component path, where exactness is not needed.</p>

<div class="danger"><b>Pope&rsquo;s tau goes blind under contamination.</b> It
divides by \(s_0\), itself estimated from data that contain the blunders. At ~10%
contamination a 100&sigma; blunder can be undetectable, and iterating does not
rescue it. Use a robust scale (<code>'test','mad'</code>) and read the
contamination gauge before trusting any tau.</div>

<h3>Variance components</h3>
<p>When several groups (GPS positions and KBR range-rates, or two satellites)
enter one adjustment, their relative weighting is rarely known a priori.
Helmert&rsquo;s VCE estimates a factor per group,
\(\hat\sigma_g^2=v_g^{\!\top}P_g v_g/r_g\) with
\(r_g=n_g-\operatorname{tr}(N^{-1}N_g)\), updates \(P_g\leftarrow\hat\sigma_g^{-2}P_g\),
and re-solves until the factors settle at unity. Only ratios are identifiable; one
group fixes the gauge [6,7,3,9,10]. For \(u\) in the thousands the trace is
estimated matrix-free (<a href="lsqvce.vceHutchinson.html">lsqvce.vceHutchinson</a>).</p>

<h3>Block-arrow systems and coloured noise</h3>
<p>Gravity recovery combines \(K\) arcs sharing \(n_c\) common parameters but each
carrying \(n_l\) local ones; ordering common-first gives a bordered block-diagonal
(&ldquo;arrow&rdquo;) normal matrix, and each local block is eliminated by its
Schur complement, leaving only the \(n_c\times n_c\) reduced system in core
(full algebra in the <a href="lsqvce_reference_manual.html#refs">Reference Manual
&sect;4</a>). Coloured noise is whitened by an AR/VAR FIR recursion &mdash; a banded
operator equal to the dense block-Cholesky GLS solve, never the dense \(n\times n\)
matrix &mdash; and separable space&times;time covariance by a Kronecker whitener
(&sect;5). Orders are chosen by <a href="lsqvce.selectArAcv.html">lsqvce.selectArAcv</a>
and validated by Ljung&ndash;Box [14&ndash;17].</p>

<h2 id="p2">Part II &mdash; What you supply</h2>
<p>A batch function <code>batchFcn(k)</code> returns the k-th batch:</p>
<table>
<tr><th>Return</th><th>Size</th><th>Meaning</th></tr>
<tr><td class="name">A</td><td>nk &times; u</td><td>design rows for this batch</td></tr>
<tr><td class="name">l</td><td>nk &times; 1</td><td>observations for this batch</td></tr>
<tr><td class="name">p</td><td>see below</td><td>WEIGHTS (not variances) for this batch</td></tr>
</table>
<table>
<tr><th>Form of p</th><th>Interpretation</th></tr>
<tr><td>scalar</td><td>one weight applied to every row of the batch</td></tr>
<tr><td>vector (nk&times;1)</td><td>per-observation weights \(=1/\sigma_i^2\) (diagonal P)</td></tr>
<tr><td>struct (cov model)</td><td>a coloured-noise model; lsqvce whitens with its factor</td></tr>
</table>
<div class="danger"><b>Weights, not variances.</b> \(p_i=1/\sigma_i^2\). Passing
\(\sigma_i^2\) silently inverts the weighting: your best observations are
down-weighted and every downstream statistic is wrong while still looking
plausible. It is dimensionally valid, so the toolbox cannot detect it.</div>
<p><b>The contract.</b> <i>Deterministic</i> &mdash; <code>batchFcn(k)</code> is
called twice (accumulate, then diagnose) and must return byte-identical data both
times; a mismatch is an error, not a silent misattribution. <i>Self-contained</i>
&mdash; a batch returns all of its own rows. <i>Complete</i> &mdash; every
observation appears in exactly one batch, \(\sum n_k=n\).</p>

<h2 id="p3">Part III &mdash; Reading the diagnosis</h2>
<p>The report is ordered so a failure at one level invalidates the levels below.
Read top to bottom, stop at the first red.</p>
<table>
<tr><th>Level</th><th>Question</th><th>Key output</th></tr>
<tr><td>0 &middot; Integrity</td><td>Does the arithmetic hold?</td><td>\(\sum r_i=n-u\), weights finite</td></tr>
<tr><td>1 &middot; Global</td><td>Is the model consistent?</td><td>\(s_0^2/\sigma_0^2\) vs \(\chi^2\)</td></tr>
<tr><td>2 &middot; Localise</td><td>Where is it broken?</td><td>data snooping, contamination gauge</td></tr>
<tr><td>3 &middot; Reliability</td><td>What can you not check?</td><td>redundancy numbers, MDB, quadrants</td></tr>
<tr><td>4 &middot; Parameters</td><td>What do you actually know?</td><td>correlations, conditioning</td></tr>
</table>
<p>The <a href="lsqvce_interpretation_guide.html">Interpretation Guide</a> covers
the reliability quadrants and the contamination gauge in detail. What each tool
adds: <code>lsqvce.systematic(D)</code> tests residuals for structure a
single-outlier search would miss; <code>lsqvce.snoop</code> is iterative data
snooping with the contamination guard; <code>lsqvce.robust</code> is M-estimation
for diffuse contamination; <code>lsqvce.vce</code>/<code>lsqvce.varComp</code>
rebalance mis-weighted groups.</p>

<h2 id="p4">Part IV &mdash; Gravity recovery: the block-arrow path</h2>
<p>The modern one-call front door is <a href="lsqvce.adjustBlock.html">lsqvce.adjustBlock</a>
(single GLS solve) and <a href="lsqvce.adjustVCE.html">lsqvce.adjustVCE</a> (the
iterative SHAPE+SCALE loop). Build loaders with
<a href="lsqvce.makeLoader.html">lsqvce.makeLoader</a> and hand them over:</p>
<pre class="example">% One-call block-arrow solve + five-level verdict
R = lsqvce.adjustBlock(nc, {Lc, Ld}, days);   % loaders from makeLoader
R.report();

% Iterative: discover the noise model (SHAPE) and group weights (SCALE)
R = lsqvce.adjustVCE(nc, {Lc, Ld}, days) ...
        .vceGroupsBy('obsType') ...
        .shape('maxOrder', 60, 'acvTol', 1e-3) ...
        .scale('maxIter', 8, 'tol', 1e-3) ...
        .run();
R.report();
disp(R.history.sigma);        % per-pass variance components</pre>
<p>Under the hood these drive a <a href="lsqvce.BlockAccumulator.html">lsqvce.BlockAccumulator</a>
(kept reachable on <code>R.acc</code> for experts). The diagnostics bridge
<a href="lsqvce.diagnoseAccum.html">lsqvce.diagnoseAccum</a> rebuilds the full
block-arrow system from the retained per-arc Schur factors to count the redundancy
the local parameters absorb &mdash; without it the common-parameter redundancies
(and every reliability number) are inflated.</p>
<div class="warn"><b>Block-arrow traps.</b> (1) Diagnose through the bridge, never
the reduced common system alone. (2) VCE identifies only ratios &mdash; fix one
group. (3) An arc whose local block is near-singular (too few observations to
determine its own state) poisons the Schur complement; the accumulator flags it
rather than silently absorbing it into the field.</div>

<h2 id="p5">Part V &mdash; Step-by-step</h2>
<h3>A &middot; Flat adjustment (single system)</h3>
<ol>
<li><b>Wrap your data as a batch function.</b> Decide the batching (by arc, day,
file); make <code>batchFcn(k)</code> return <code>[A,l,p]</code>. Keep it deterministic.</li>
<li><b>Get the weights right.</b> \(p=1/\sigma^2\). For coloured noise pass a
covariance-model struct instead of a vector and let lsqvce whiten.</li>
<li><b>Solve and judge in one call.</b>
<code>[sol, D, V] = lsqvce.adjust(u, batchFcn, nBatch);</code></li>
<li><b>Read the report top-down.</b> Stop at the first red. A Level&nbsp;0 red is a
data/weighting bug, not a blunder &mdash; fix it before reading on.</li>
<li><b>If Level&nbsp;1 fails, localise before rejecting.</b> Run
<code>lsqvce.systematic(D)</code> first: a structured error is a model defect, and
snooping it away mutilates good data.</li>
<li><b>Snoop, gauge open.</b>
<code>[sol,D,REJ] = lsqvce.snoop(sol, D, batchFcn, 'test','mad');</code></li>
<li><b>Before believing the parameters, read reliability.</b> Check MDB and the
quadrants for every parameter you care about. A clean residual at \(r_i\approx0\)
proves nothing.</li>
</ol>
<h3>B &middot; Gravity recovery (many arcs, VCE)</h3>
<ol>
<li>Write one <a href="lsqvce.makeLoader.html">makeLoader</a> loader per
satellite/arc returning whitened batches; give the common count <code>nc</code>
and per-system <code>numLocal</code>.</li>
<li>For a single solve call <a href="lsqvce.adjustBlock.html">lsqvce.adjustBlock</a>;
to discover the noise model and group weights, use
<a href="lsqvce.adjustVCE.html">lsqvce.adjustVCE</a> with
<code>.vceGroupsBy(...).shape(...).scale(...).run()</code>.</li>
<li>Inspect <code>R.history.sigma</code> for a group being systematically up- or
down-weighted, and confirm convergence.</li>
<li><b>Trust the field only where it is checked.</b> Map the common-parameter
redundancies; regions carried by a single arc are blind spots regardless of how
small their residuals are.</li>
</ol>
<div class="note"><b>Validate your build.</b>
<code>results = run_all_lsqvce();</code> runs the whole suite. The Python mirrors
under <code>validation/</code> reproduce the numerical properties the MATLAB tests
assert (redundancy identity, whitening equals dense GLS, ACV round-trip).</div>

<div class="seealso"><b>See also:</b>
<a href="lsqvce_reference_manual.html">Reference Manual</a>,
<a href="lsqvce_interpretation_guide.html">Interpretation Guide</a>,
<a href="lsqvce_choosing.html">Choosing an engine</a>,
<a href="lsqvce.adjust.html">lsqvce.adjust</a>,
<a href="lsqvce.adjustVCE.html">lsqvce.adjustVCE</a></div>
''' + REFERENCES


# ---- emit --------------------------------------------------------------
# STOP. This generator is currently BEHIND its own output (HANDOVER 1.3b).
#   * user-manual Parts IV-B / IV-C (the interface guide, task 2) exist ONLY in
#     html/lsqvce_user_manual.html -- running this would DELETE them;
#   * the bodies below still contain the pre-rename fluent methods
#     (vceGroupsBy / .shape( / .scale() and two stale `nc` examples, which were
#     fixed in the OUTPUT -- running this would REINTRODUCE them, and
#     audit_tree.py check 13 would then fire on the regenerated HTML.
# THIS FILE WAS NOT RETIRED WITH THE OLD PDFs, and that is deliberate. The master
# PDF (docs_src/build_master_pdf.py) replaced the three PDF manuals, which were
# DELETED. The four HTML manuals this generates are still SHIPPED AND LIVE --
# lsqvce_{user_manual,reference_manual,interpretation_guide,choosing}.html, all
# linked from helptoc.xml and html/index.html. They are what `doc lsqvce` opens.
#
# So this is not obsolete; it is a trap with a guard on it. Running it makes the
# doc set WORSE, so it refuses. Clear _STALE_GUARD only after back-porting
# Parts IV-B/IV-C and the fluent renames into the USERMAN body below, and then
# re-run audit_tree.py checks 13 and 15 on the regenerated HTML.
_STALE_GUARD = True

if __name__ == '__main__' and _STALE_GUARD:
    raise SystemExit(
        'build_manuals.py is STALE and would revert shipped fixes.\n'
        'It is NOT obsolete: the four HTML manuals it generates are still live '
        'and linked from doc lsqvce -- only the three old PDFs were retired.\n'
        'See HANDOVER 1.4C. Set _STALE_GUARD = False only after back-porting '
        'user-manual Parts IV-B/IV-C and the fluent-method renames.')

if __name__ == '__main__':
    write('html/lsqvce_interpretation_guide.html',
          'lsqvce Interpretation Guide', 'Interpretation Guide',
          INTERP, 'Interpretation Guide moved from PDF to HTML.')
    write('html/lsqvce_reference_manual.html',
          'lsqvce Reference Manual', 'Reference Manual',
          REFMAN, 'Reference Manual moved from PDF to HTML (equations rendered from LaTeX to static PNG images).')
    write('html/lsqvce_user_manual.html',
          'lsqvce User Manual', 'User Manual',
          USERMAN, 'User Manual moved from PDF to HTML (equations rendered from LaTeX to static PNG images; names updated).')
