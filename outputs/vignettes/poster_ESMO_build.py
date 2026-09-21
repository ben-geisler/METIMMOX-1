# -*- coding: utf-8 -*-
"""Build outputs/vignettes/poster_ESMO.pptx (ESMO 2026 poster 4355P).

Unpacks outputs/vignettes/poster_ESMO_example.pptx, reuses its package (theme,
master, layouts, UiO and Ahus logos) and writes a new ppt/slides/slide1.xml with
the same box styles (rounded rectangles with the 0093C8 border, solid header
bars, Calibri body, Arial title) but the 4355P content. The three embedded
figures come from data/output/poster_ESMO/ (run poster_ESMO_figures.R first).
Requires Python 3 with only the standard library.

    Rscript outputs/vignettes/poster_ESMO_figures.R
    python outputs/vignettes/poster_ESMO_build.py
"""
import os
import re
import shutil
import zipfile
from xml.sax.saxutils import escape

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
EXAMPLE_PPTX = os.path.join(ROOT, "outputs", "vignettes", "poster_ESMO_example.pptx")
WORK = os.environ.get("POSTER_OUT", os.path.join(ROOT, "data", "output", "poster_ESMO"))
EX = os.path.join(WORK, "example")
BUILD = os.path.join(WORK, "build")
FIGS = WORK
OUT = os.path.join(ROOT, "outputs", "vignettes", "poster_ESMO.pptx")

E = 914400


def emu(inches):
    return int(round(inches * E))


TITLE_SZ = int(os.environ.get("TITLE_SZ", "60"))
BODY_SZ = int(os.environ.get("BODY_SZ", "32"))
AUTHOR_SZ = 44
AFFIL_SZ = 32
TABLE_SZ = 28
FOOT_SZ = 24
CAPTION_SZ = 24

ARIAL = ('<a:latin typeface="Arial" panose="020B0604020202020204" pitchFamily="34" charset="0"/>'
         '<a:cs typeface="Arial" panose="020B0604020202020204" pitchFamily="34" charset="0"/>')

_next_id = [100]


def nid():
    _next_id[0] += 1
    return _next_id[0]


def t(text):
    s = escape(text)
    if text != text.strip():
        return '<a:t xml:space="preserve">%s</a:t>' % s
    return "<a:t>%s</a:t>" % s


def run(text, sz, b=False, i=False, sup=False, color=None, font=None, lang="en-US"):
    attrs = 'lang="%s" sz="%d"' % (lang, sz * 100)
    if b:
        attrs += ' b="1"'
    if i:
        attrs += ' i="1"'
    if sup:
        attrs += ' baseline="30000"'
    inner = ""
    if color:
        inner += '<a:solidFill><a:srgbClr val="%s"/></a:solidFill>' % color
    if font == "Arial":
        inner += ARIAL
    return "<a:r><a:rPr %s dirty=\"0\">%s</a:rPr>%s</a:r>" % (attrs, inner, t(text))


def para(runs, algn=None, bullet=False, spc_after=None, marL=None, marR=None, indent=None,
         end_sz=None, lnSpc=None):
    ppr_attrs = ""
    if bullet:
        ppr_attrs += ' marL="571500" indent="-571500"'
    else:
        if marL is not None:
            ppr_attrs += ' marL="%d"' % marL
        if indent is not None:
            ppr_attrs += ' indent="%d"' % indent
    if marR is not None:
        ppr_attrs += ' marR="%d"' % marR
    if algn:
        ppr_attrs += ' algn="%s"' % algn
    inner = ""
    if lnSpc is not None:
        inner += '<a:lnSpc><a:spcPct val="%d"/></a:lnSpc>' % lnSpc
    if spc_after is not None:
        inner += '<a:spcAft><a:spcPts val="%d"/></a:spcAft>' % (spc_after * 100)
    if bullet:
        inner += ('<a:buFont typeface="Courier New" panose="02070309020205020404" '
                  'pitchFamily="49" charset="0"/><a:buChar char="o"/>')
    ppr = "<a:pPr%s>%s</a:pPr>" % (ppr_attrs, inner) if (ppr_attrs or inner) else ""
    end = '<a:endParaRPr lang="en-US" sz="%d" dirty="0"/>' % (end_sz * 100) if end_sz else ""
    return "<a:p>%s%s%s</a:p>" % (ppr, "".join(runs), end)


def xfrm(x, y, w, h):
    return ('<a:xfrm><a:off x="%d" y="%d"/><a:ext cx="%d" cy="%d"/></a:xfrm>'
            % (emu(x), emu(y), emu(w), emu(h)))


BORDER_STYLE = ('<p:style><a:lnRef idx="2"><a:schemeClr val="accent1"/></a:lnRef>'
                '<a:fillRef idx="1"><a:schemeClr val="lt1"/></a:fillRef>'
                '<a:effectRef idx="0"><a:schemeClr val="accent1"/></a:effectRef>'
                '<a:fontRef idx="minor"><a:schemeClr val="dk1"/></a:fontRef></p:style>')


def round_box(name, x, y, w, h, paras, anchor="ctr", lIns=None, rIns=None, tIns=None, bIns=None):
    """Rounded rectangle with the 0093C8 border and transparent fill (example shape 13)."""
    ins = ""
    if lIns is not None:
        ins += ' lIns="%d"' % emu(lIns)
    if rIns is not None:
        ins += ' rIns="%d"' % emu(rIns)
    if tIns is not None:
        ins += ' tIns="%d"' % emu(tIns)
    if bIns is not None:
        ins += ' bIns="%d"' % emu(bIns)
    return (
        '<p:sp><p:nvSpPr><p:cNvPr id="%d" name="%s"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>'
        '<p:spPr>%s<a:prstGeom prst="roundRect"><a:avLst/></a:prstGeom>'
        '<a:solidFill><a:schemeClr val="lt1"><a:alpha val="0"/></a:schemeClr></a:solidFill>'
        '<a:ln w="28575"><a:solidFill><a:srgbClr val="0093C8"/></a:solidFill></a:ln></p:spPr>%s'
        '<p:txBody><a:bodyPr%s rtlCol="0" anchor="%s"><a:noAutofit/></a:bodyPr><a:lstStyle/>%s</p:txBody></p:sp>'
        % (nid(), escape(name), xfrm(x, y, w, h), BORDER_STYLE, ins, anchor, "".join(paras))
    )


def header_box(name, x, y, w, h, text, sz=66):
    """Solid 0093C8 rounded header bar with white bold Calibri text (example shape 21)."""
    return (
        '<p:sp><p:nvSpPr><p:cNvPr id="%d" name="%s"/><p:cNvSpPr/><p:nvPr/></p:nvSpPr>'
        '<p:spPr>%s<a:prstGeom prst="roundRect"><a:avLst/></a:prstGeom>'
        '<a:solidFill><a:srgbClr val="0093C8"/></a:solidFill>'
        '<a:ln w="12700" cap="flat" cmpd="sng" algn="ctr"><a:solidFill><a:srgbClr val="0093C8"/></a:solidFill>'
        '<a:prstDash val="solid"/><a:miter lim="800000"/></a:ln><a:effectLst/></p:spPr>'
        '<p:txBody><a:bodyPr rtlCol="0" anchor="ctr"/><a:lstStyle/>'
        '<a:p><a:pPr algn="ctr"/><a:r><a:rPr lang="en-US" sz="%d" b="1" dirty="0">'
        '<a:solidFill><a:prstClr val="white"/></a:solidFill>'
        '<a:latin typeface="Calibri" panose="020F0502020204030204"/><a:ea typeface="+mn-ea"/><a:cs typeface="+mn-cs"/>'
        '</a:rPr>%s</a:r></a:p></p:txBody></p:sp>'
        % (nid(), escape(name), xfrm(x, y, w, h), sz * 100, t(text))
    )


def text_box(name, x, y, w, h, paras, anchor="t", autofit=False, lIns=None, rIns=None):
    ins = ""
    if lIns is not None:
        ins += ' lIns="%d"' % emu(lIns)
    if rIns is not None:
        ins += ' rIns="%d"' % emu(rIns)
    fit = "<a:spAutoFit/>" if autofit else "<a:noAutofit/>"
    return (
        '<p:sp><p:nvSpPr><p:cNvPr id="%d" name="%s"/><p:cNvSpPr txBox="1"/><p:nvPr/></p:nvSpPr>'
        '<p:spPr>%s<a:prstGeom prst="rect"><a:avLst/></a:prstGeom><a:noFill/></p:spPr>'
        '<p:txBody><a:bodyPr wrap="square"%s rtlCol="0" anchor="%s">%s</a:bodyPr><a:lstStyle/>%s</p:txBody></p:sp>'
        % (nid(), escape(name), xfrm(x, y, w, h), ins, anchor, fit, "".join(paras))
    )


def pic(name, rid, x, y, w, h, src_rect=None, descr=None):
    sr = ""
    if src_rect:
        sr = "<a:srcRect %s/>" % " ".join('%s="%d"' % (k, v) for k, v in src_rect.items())
    d = ' descr="%s"' % escape(descr) if descr else ""
    return (
        '<p:pic><p:nvPicPr><p:cNvPr id="%d" name="%s"%s/><p:cNvPicPr><a:picLocks noChangeAspect="1"/></p:cNvPicPr>'
        '<p:nvPr/></p:nvPicPr><p:blipFill><a:blip r:embed="%s"/>%s<a:stretch><a:fillRect/></a:stretch></p:blipFill>'
        '<p:spPr>%s<a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr></p:pic>'
        % (nid(), escape(name), d, rid, sr, xfrm(x, y, w, h))
    )


# ----------------------------------------------------------------------------
# Table
# ----------------------------------------------------------------------------
def tc(cell, ncols_span=1, hmerge=False, rowspan=1, vmerge=False):
    """cell: dict(text|runs, b, algn, fill, top, bottom, i)"""
    attrs = ""
    if ncols_span > 1:
        attrs += ' gridSpan="%d"' % ncols_span
    if hmerge:
        attrs += ' hMerge="1"'
    if rowspan > 1:
        attrs += ' rowSpan="%d"' % rowspan
    if vmerge:
        attrs += ' vMerge="1"'
    runs = cell.get("runs")
    if runs is None:
        txt = cell.get("text", "")
        runs = [run(txt, TABLE_SZ, b=cell.get("b", False), i=cell.get("i", False))] if txt else []
    p = para(runs, algn=cell.get("algn"), end_sz=TABLE_SZ)
    lines = ""
    for side in ("L", "R"):
        lines += "<a:ln%s><a:noFill/></a:ln%s>" % (side, side)
    for side, key in (("T", "top"), ("B", "bottom")):
        wdt = cell.get(key)
        if wdt:
            lines += ('<a:ln%s w="%d" cap="flat" cmpd="sng" algn="ctr"><a:solidFill><a:srgbClr val="000000"/>'
                      '</a:solidFill><a:prstDash val="solid"/><a:round/></a:ln%s>' % (side, wdt, side))
        else:
            lines += "<a:ln%s><a:noFill/></a:ln%s>" % (side, side)
    fill = ('<a:solidFill><a:srgbClr val="%s"/></a:solidFill>' % cell["fill"]) if cell.get("fill") else "<a:noFill/>"
    return ('<a:tc%s><a:txBody><a:bodyPr/><a:lstStyle/>%s</a:txBody>'
            '<a:tcPr marL="60000" marR="60000" marT="18000" marB="18000" anchor="ctr">%s%s</a:tcPr></a:tc>'
            % (attrs, p, lines, fill))


def table(name, x, y, col_widths, rows, row_h):
    """rows: list of lists of tc() strings (already merged/spanned)."""
    grid = "".join('<a:gridCol w="%d"/>' % emu(w) for w in col_widths)
    trs = "".join('<a:tr h="%d">%s</a:tr>' % (emu(row_h), "".join(r)) for r in rows)
    h = row_h * len(rows)
    return (
        '<p:graphicFrame><p:nvGraphicFramePr><p:cNvPr id="%d" name="%s"/>'
        '<p:cNvGraphicFramePr><a:graphicFrameLocks noGrp="1"/></p:cNvGraphicFramePr><p:nvPr/></p:nvGraphicFramePr>'
        '<p:xfrm><a:off x="%d" y="%d"/><a:ext cx="%d" cy="%d"/></p:xfrm>'
        '<a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/table">'
        '<a:tbl><a:tblPr firstRow="1" bandRow="1"/><a:tblGrid>%s</a:tblGrid>%s</a:tbl>'
        '</a:graphicData></a:graphic></p:graphicFrame>'
        % (nid(), escape(name), emu(x), emu(y), emu(sum(col_widths)), emu(h), grid, trs)
    )


# ----------------------------------------------------------------------------
# Content
# ----------------------------------------------------------------------------
X = "\u00d7"      # multiplication sign
GE = "\u2265"     # >=
ARR = "\u2192"    # ->
EN = "\u2013"     # en dash

shapes = []

# --- Title box -----------------------------------------------------------------
title_p1 = para([run("Biomarker signals in MSS/pMMR metastatic colorectal cancer (mCRC) are predominantly "
                     "prognostic, not predictive", TITLE_SZ, b=True, font="Arial")], algn="ctr")
title_p2 = para([run(EN + " lessons from the randomized phase-2 METIMMOX trial", TITLE_SZ, b=True,
                     font="Arial")], algn="ctr", spc_after=6)

authors = [("Benjamin P. Geisler", "1,2@"), ("Megan Othus", "3"), ("Sebastian Meltzer", "2"),
           ("Paula A. Bousquet", "2"), ("Emily A. Burger", "1,4"), ("Eline Aas", "1,5"),
           ("Anne H. Ree", "2,6")]
author_runs = []
for k, (nm, aff) in enumerate(authors):
    if k > 0:
        author_runs.append(run(", ", AUTHOR_SZ, font="Arial"))
    author_runs.append(run(nm, AUTHOR_SZ, font="Arial"))
    author_runs.append(run(aff, AUTHOR_SZ, sup=True, font="Arial"))
author_p = para(author_runs, algn="ctr", spc_after=4)

affils = [("1", "Department of Health Management and Health Economics, University of Oslo, Oslo, Norway; "),
          ("2", "Department of Oncology, Akershus University Hospital, L\u00f8renskog, Norway; "),
          ("3", "Fred Hutchinson Cancer Center, Seattle WA, United States; "),
          ("4", "Center for Health Decision Science, Harvard T. H. Chan School of Public Health, Boston MA, United States; "),
          ("5", "Division for Health Services, Norwegian Institute of Public Health, Oslo, Norway; "),
          ("6", "Institute of Clinical Medicine, University of Oslo, Oslo, Norway.")]
affil_runs = []
for num, txt in affils:
    affil_runs.append(run(num, AFFIL_SZ, sup=True, font="Arial"))
    affil_runs.append(run(txt, AFFIL_SZ, font="Arial"))
affil_p = para(affil_runs, algn="ctr", marL=emu(3.0), marR=emu(3.0))

shapes.append(round_box("Title box", 0.35, 0.31, 55.29, 4.22,
                        [title_p1, title_p2, author_p, affil_p], anchor="ctr"))

# --- Poster number (example: 72pt bold B41672) ----------------------------------
shapes.append(text_box("Poster number", 0.86, 3.22, 2.9, 1.31,
                       [para([run("4355P", 72, b=True, color="B41672")])], autofit=True))

# --- Section headers ------------------------------------------------------------
COL_A_X, COL_A_W = 0.31, 15.30      # left column (Background, Methods)
COL_B_X, COL_B_W = 16.11, 16.26     # middle column (Results text, forest plot)
COL_C_X, COL_C_W = 32.87, 22.72     # right column (Table 1, Kaplan-Meier)
shapes.append(header_box("Background header", COL_A_X, 4.89, COL_A_W, 1.08, "Background"))
shapes.append(header_box("Results header", COL_B_X, 4.88, 55.59 - COL_B_X, 1.10, "Results"))

# --- Background -----------------------------------------------------------------
BS = BODY_SZ
SA = 12  # paragraph spacing after (pt)


def bullet(runs):
    return para(runs, bullet=True, spc_after=SA, end_sz=BS)


bg = [
    bullet([run("Immune checkpoint blockade is ineffective in unselected MSS/pMMR mCRC.", BS)]),
    bullet([run("The randomized phase-2 METIMMOX trial (NCT03388190; EudraCT 2017-001845-29) compared alternating "
                "oxaliplatin-based chemotherapy (FLOX) and nivolumab with FLOX alone; progression-free survival "
                "(PFS) was identical.", BS)]),
    bullet([run("Post hoc analyses suggested benefit in three biomarker subgroups but conflated prognostic and "
                "predictive effects.", BS)]),
    bullet([run("We evaluated treatment-effect heterogeneity in a causally motivated regression framework.", BS)]),
]
shapes.append(round_box("Background box", 0.40, 6.29, COL_A_W - 0.09, 5.50, bg, anchor="ctr"))

# --- Methods --------------------------------------------------------------------
shapes.append(header_box("Methods header", COL_A_X, 12.05, COL_A_W, 1.07, "Methods"))
me = [
    bullet([run("Complete-case cohort: 68 patients with both biomarkers, covariates and both endpoints "
                "(59 deaths; 63 progression-or-death events).", BS)]),
    bullet([run("Pre-immunotherapy biomarkers: TMB " + GE + "9 mut/Mb or ", BS), run("BRAF", BS, i=True),
            run(" V600E (baseline NGS), and suppressed systemic inflammation (CRP <5 mg/L) after the initial "
                "2 FLOX cycles common to both arms (week 4, before the first nivolumab dose).", BS)]),
    bullet([run("Tumor lesion reduction (TLR) " + GE + "10% at the first on-treatment CT was adjudicated an "
                "intermediate endpoint (T " + ARR + " TLR " + ARR + " PFS) and analyzed separately with a "
                "week-9 landmark.", BS)]),
    bullet([run("One unified Firth-corrected Cox model per endpoint (PFS and overall survival [OS]; ", BS),
            run("Figure 1", BS, b=True),
            run("): age, sex, treatment, biomarker main effects and treatment " + X + " biomarker interactions; "
                "profile-likelihood 95% CIs and penalized likelihood-ratio tests (PLRT).", BS)]),
    bullet([run("Ridge-penalized Cox regression (cross-validated penalty) as co-primary stability estimator.", BS)]),
]
METH_W = COL_A_W - 0.05
shapes.append(round_box("Methods box", 0.36, 13.40, METH_W, 15.45, me, anchor="t", tIns=0.12))
DAG_W, DAG_H = 6.8, 6.8 * (4.8 / 7.0)
dag_x = 0.36 + (METH_W - DAG_W) / 2
DAG_CAP_H = 0.85
dag_y = 28.85 - 0.05 - DAG_CAP_H - DAG_H
shapes.append(pic("Simplified DAG", "rId7", dag_x, dag_y, DAG_W, DAG_H,
                  descr="Simplified causal diagram: treatment, CRP, TMB/BRAF, age/sex and survival"))
shapes.append(text_box("DAG caption", 1.24, dag_y + DAG_H, 13.5, DAG_CAP_H, [
    para([run("Figure 1. Simplified causal diagram. ", CAPTION_SZ, b=True),
          run("Solid arrows: assumed causal effects; dashed arrow: treatment effect hypothesized to be "
              "modified by CRP and TMB/", CAPTION_SZ), run("BRAF", CAPTION_SZ, i=True),
          run(" (interaction terms).", CAPTION_SZ)], algn="ctr")], anchor="t"))

# --- Results text ---------------------------------------------------------------
re_ = [
    bullet([run("68 patients had complete biomarker and covariate data (mean age 64 years; 47% female; "
                "59 deaths; 63 progression-or-death events).", BS)]),
    bullet([run("No detectable marginal treatment effect: OS HR 0.96 (95% CI 0.58" + EN + "1.61); "
                "PFS HR 0.80 (0.49" + EN + "1.33).", BS)]),
    bullet([run("CRP <5 mg/L and TMB/", BS), run("BRAF", BS, i=True),
            run(" positivity were prognostic: both were associated with longer PFS (HR 0.41 and 0.49), "
                "and low CRP also with longer OS (HR 0.49) (", BS), run("Table 1", BS, b=True), run(", ", BS),
            run("Figure 2", BS, b=True), run(").", BS)]),
    bullet([run("No treatment " + X + " biomarker interaction reached statistical significance (", BS),
            run("Figure 3", BS, b=True),
            run("). The CRP " + X + " treatment interaction for PFS was the most directionally consistent "
                "finding (Firth HR 0.42, 95% CI 0.11" + EN + "1.68; PLRT p = 0.21; ridge HR 0.33) but was "
                "attenuated for OS (Firth HR 0.65; ridge HR 0.40).", BS)]),
    bullet([run("Week-4 CRP positivity was imbalanced by chance (47% with FLOX/nivolumab vs 19% with FLOX "
                "alone; p = 0.02); the CRP " + X + " treatment interaction is estimated against this contrast.", BS)]),
    bullet([run("TLR " + GE + "10% was strongly associated with both endpoints (landmark HR 0.35 for OS and "
                "PFS) but was more frequent with FLOX alone (75.9% vs 52.8% with FLOX/nivolumab; p = 0.07), "
                "and the TLR " + X + " treatment interaction pointed towards less, not more, benefit from "
                "nivolumab (HR >1): TLR marks chemotherapy-responsive disease rather than immunotherapy benefit.", BS)]),
    bullet([run("Proportional hazards held globally (Schoenfeld p = 0.25 for OS, 0.12 for PFS).", BS)]),
]
shapes.append(round_box("Results text box", COL_B_X, 6.29, COL_B_W, 13.62, re_, anchor="ctr"))

# --- Figure 2 (forest) ----------------------------------------------------------
FOR_W = COL_B_W
FOR_H = FOR_W * (8.5 / 16.0)
shapes.append(pic("Figure 3 forest plot", "rId6", COL_B_X, 28.85 - FOR_H, FOR_W, FOR_H,
                  descr="Forest plot of treatment by biomarker interaction hazard ratios"))

# --- Table 1 --------------------------------------------------------------------
TX, TY, TW = COL_C_X, 6.29, COL_C_W
shapes.append(text_box("Table 1 caption", TX, TY, TW, 0.66, [
    para([run("Table 1. ", TABLE_SZ, b=True),
          run("Hazard ratios (95% profile-likelihood CI) from Firth-corrected Cox regression", TABLE_SZ)])],
    anchor="t"))

col_w = [7.42, 4.00, 1.60, 2.00, 4.00, 1.60, 2.10]
assert abs(sum(col_w) - TW) < 1e-6
HDR = "D9D9D9"
GRP = "EFEFEF"
THICK, THIN = 19050, 12700


def hdr(text, span=1, b=True, algn="ctr", fill=HDR, top=None, bottom=None, hmerge=False, rowspan=1, vmerge=False):
    return tc({"text": text, "b": b, "algn": algn, "fill": fill, "top": top, "bottom": bottom},
              ncols_span=span, hmerge=hmerge, rowspan=rowspan, vmerge=vmerge)


def group_row(text):
    return [tc({"text": text, "b": True, "algn": "l", "fill": GRP}, ncols_span=7)] + \
           [tc({"fill": GRP}, hmerge=True) for _ in range(6)]


def data_row(term, os_hr, os_p, os_r, pfs_hr, pfs_p, pfs_r, bottom=None, term_runs=None):
    cells = [tc({"text": term, "runs": term_runs, "algn": "l", "bottom": bottom})]
    for v in (os_hr, os_p, os_r, pfs_hr, pfs_p, pfs_r):
        cells.append(tc({"text": v, "algn": "ctr", "bottom": bottom}))
    return cells


rows = [
    [hdr("Term", rowspan=2, algn="l", top=THICK), hdr("Overall survival", span=3, top=THICK),
     hdr("", hmerge=True, top=THICK), hdr("", hmerge=True, top=THICK),
     hdr("Progression-free survival", span=3, top=THICK),
     hdr("", hmerge=True, top=THICK), hdr("", hmerge=True, top=THICK)],
    [hdr("", vmerge=True, bottom=THIN), hdr("HR (95% CI)", bottom=THIN), hdr("p", bottom=THIN),
     hdr("Ridge HR", bottom=THIN), hdr("HR (95% CI)", bottom=THIN), hdr("p", bottom=THIN),
     hdr("Ridge HR", bottom=THIN)],
    group_row("Marginal associations (one term per model, n = 68)"),
    data_row("Treatment (FLOX/nivolumab vs FLOX)", "0.96 (0.58" + EN + "1.61)", "0.87", EN,
             "0.80 (0.49" + EN + "1.33)", "0.39", EN),
    data_row("CRP <5 mg/L at week 4", "0.49 (0.27" + EN + "0.84)", "0.009", EN,
             "0.41 (0.23" + EN + "0.71)", "0.001", EN),
    data_row("TMB " + GE + "9 mut/Mb or BRAF V600E", "0.68 (0.40" + EN + "1.13)", "0.14", EN,
             "0.49 (0.28" + EN + "0.84)", "0.009", EN,
             term_runs=[run("TMB " + GE + "9 mut/Mb or ", TABLE_SZ), run("BRAF", TABLE_SZ, i=True),
                        run(" V600E", TABLE_SZ)]),
    data_row("TLR " + GE + "10% (week-9 landmark)\u2020", "0.35 (0.21" + EN + "0.61)", "<0.001", EN,
             "0.35 (0.20" + EN + "0.63)", "<0.001", EN),
    group_row("Unified model: treatment " + X + " biomarker interactions (PLRT p)"),
    data_row("CRP " + X + " treatment", "0.65 (0.19" + EN + "2.60)", "0.52", "0.40",
             "0.42 (0.11" + EN + "1.68)", "0.21", "0.33"),
    data_row("TMB/BRAF " + X + " treatment", "0.95 (0.32" + EN + "2.82)", "0.92", "0.86",
             "0.66 (0.21" + EN + "2.05)", "0.47", "0.51",
             term_runs=[run("TMB/", TABLE_SZ), run("BRAF", TABLE_SZ, i=True), run(" " + X + " treatment", TABLE_SZ)]),
    group_row("Exploratory TLR responder model (week-9 landmark, adjusted)\u2020"),
    data_row("TLR " + X + " treatment", "2.47 (0.75" + EN + "7.75)", "0.14", EN,
             "1.88 (0.52" + EN + "6.38)", "0.33", EN, bottom=THICK),
]
ROW_H = 0.60
shapes.append(table("Table 1", TX, TY + 0.72, col_w, rows, ROW_H))
table_bottom = TY + 0.72 + ROW_H * len(rows)

foot_runs = [
    run("Unified model: ", FOOT_SZ, b=True),
    run("age + sex + treatment + CRP + TMB/", FOOT_SZ), run("BRAF", FOOT_SZ, i=True),
    run(" + CRP " + X + " treatment + TMB/", FOOT_SZ), run("BRAF", FOOT_SZ, i=True),
    run(" " + X + " treatment; p = penalized likelihood-ratio test. ", FOOT_SZ),
    run("Ridge: ", FOOT_SZ, b=True),
    run("L2-penalized Cox regression (biomarker main effects penalized, 10-fold cross-validated lambda), point "
        "estimates only. ", FOOT_SZ),
    run("\u2020", FOOT_SZ, b=True),
    run("TLR-classified patients (n = 65); time measured from the week-9 landmark (OS cohort n = 65, PFS cohort "
        "n = 60). The responder model is adjusted for age, sex, CRP and TMB/", FOOT_SZ), run("BRAF", FOOT_SZ, i=True),
    run("; because TLR is measured after randomization, its interaction is a responder-stratified contrast, not a "
        "treatment-selection estimate. ", FOOT_SZ),
    run("CI", FOOT_SZ, b=True), run(", confidence interval; ", FOOT_SZ),
    run("CRP", FOOT_SZ, b=True), run(", C-reactive protein; ", FOOT_SZ),
    run("HR", FOOT_SZ, b=True), run(", hazard ratio; ", FOOT_SZ),
    run("TLR", FOOT_SZ, b=True), run(", tumor lesion reduction; ", FOOT_SZ),
    run("TMB", FOOT_SZ, b=True), run(", tumor mutational burden.", FOOT_SZ),
]
FOOT_H = 2.05
shapes.append(round_box("Table 1 footnote", TX, table_bottom + 0.12, TW, FOOT_H,
                        [para(foot_runs, end_sz=FOOT_SZ)], anchor="ctr", lIns=0.12, rIns=0.12))
foot_bottom = table_bottom + 0.12 + FOOT_H

# --- Figure 1 (KM) --------------------------------------------------------------
KM_TOP = foot_bottom + 0.25
KM_BOTTOM = 28.85
KM_H = KM_BOTTOM - KM_TOP
KM_W = KM_H * (24.0 / 13.5)
if KM_W > TW:
    KM_W = TW
    KM_H = KM_W * (13.5 / 24.0)
km_x = TX + (TW - KM_W) / 2
shapes.append(pic("Figure 2 Kaplan-Meier", "rId5", km_x, KM_TOP, KM_W, KM_H,
                  descr="Kaplan-Meier curves by biomarker status and treatment arm"))

# --- Conclusion -----------------------------------------------------------------
shapes.append(header_box("Conclusion header", 0.40, 29.15, 55.23, 0.95, "Conclusion", sz=54))
concl = [
    para([run("In this causally structured exploratory analysis, biomarker signals in MSS/pMMR mCRC were "
              "predominantly prognostic rather than predictive. The CRP " + X + " treatment interaction for PFS "
              "warrants evaluation in adequately powered prospective studies. The causal framework offers a "
              "template for transparent post hoc biomarker analysis of phase-2 trials.", 44, b=True)],
         algn="ctr"),
]
shapes.append(round_box("Conclusion box", 0.40, 30.30, 55.23, 2.20, concl, anchor="ctr"))

# --- Footer ---------------------------------------------------------------------
shapes.append(pic("Ahus logo", "rId4", 0.37, 32.79, 10.14, 1.31, src_rect={"b": 50329},
                  descr="Akershus University Hospital Cancer Center logo"))
foot = [
    para([run("Disclosure: ", 24, b=True),
          run("B.P. Geisler has no potential conflicts of interest to declare. ", 24),
          run("Funding: ", 24, b=True),
          run("Study sponsored by European Union, Research Council of Norway, Norwegian Cancer Society, "
              "South-Eastern Regional Health Authority, and Bristol-Myers Squibb and Sigma2 AS (in-kind donations). "
              "Views and opinions expressed are, however, those of the authors only and do not necessarily reflect "
              "those of the European Union or the European Health and Digital Executive Agency (HaDEA).", 24)]),
]
shapes.append(text_box("Disclosure and funding", 10.9, 32.60, 27.3, 1.60, foot, anchor="ctr"))
shapes.append(text_box("Email", 38.3, 32.72, 10.4, 0.90, [
    para([run("@ ", 44, b=True, sup=True), run("b.p.geisler@medisin.uio.no", 44, b=True)])], anchor="ctr"))
shapes.append(pic("UiO logo", "rId3", 48.87, 32.71, 6.74, 1.13, descr="University of Oslo logo"))

# ----------------------------------------------------------------------------
# Assemble slide XML
# ----------------------------------------------------------------------------
slide_xml = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<p:sld xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" '
    'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" '
    'xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main">'
    '<p:cSld><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr>'
    '<p:grpSpPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="0" cy="0"/><a:chOff x="0" y="0"/><a:chExt cx="0" cy="0"/></a:xfrm></p:grpSpPr>'
    + "".join(shapes) +
    '</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>'
)

# ----------------------------------------------------------------------------
# Package
# ----------------------------------------------------------------------------
for f in ("poster_esmo_km.png", "poster_esmo_forest.png", "poster_esmo_dag.png"):
    if not os.path.exists(os.path.join(FIGS, f)):
        raise SystemExit("Missing %s in %s: run outputs/vignettes/poster_ESMO_figures.R first." % (f, FIGS))
if os.path.exists(EX):
    shutil.rmtree(EX)
with zipfile.ZipFile(EXAMPLE_PPTX) as zf:
    zf.extractall(EX)
if os.path.exists(BUILD):
    shutil.rmtree(BUILD)
shutil.copytree(EX, BUILD)

media = os.path.join(BUILD, "ppt", "media")
for old in ("image3.emf", "image4.png", "image5.png"):
    os.remove(os.path.join(media, old))
shutil.copy(os.path.join(FIGS, "poster_esmo_km.png"), os.path.join(media, "image3.png"))
shutil.copy(os.path.join(FIGS, "poster_esmo_forest.png"), os.path.join(media, "image4.png"))
shutil.copy(os.path.join(FIGS, "poster_esmo_dag.png"), os.path.join(media, "image5.png"))

with open(os.path.join(BUILD, "ppt", "slides", "slide1.xml"), "w", encoding="utf-8") as fh:
    fh.write(slide_xml)

rels = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/slideLayout" Target="../slideLayouts/slideLayout1.xml"/>'
    '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/notesSlide" Target="../notesSlides/notesSlide1.xml"/>'
    '<Relationship Id="rId3" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="../media/image1.png"/>'
    '<Relationship Id="rId4" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="../media/image2.png"/>'
    '<Relationship Id="rId5" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="../media/image3.png"/>'
    '<Relationship Id="rId6" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="../media/image4.png"/>'
    '<Relationship Id="rId7" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" Target="../media/image5.png"/>'
    '</Relationships>'
)
with open(os.path.join(BUILD, "ppt", "slides", "_rels", "slide1.xml.rels"), "w", encoding="utf-8") as fh:
    fh.write(rels)

# content types: drop the emf default (no emf left in the package)
ct_path = os.path.join(BUILD, "[Content_Types].xml")
with open(ct_path, encoding="utf-8") as fh:
    ct = fh.read()
ct = ct.replace('<Default Extension="emf" ContentType="image/x-emf"/>', "")
with open(ct_path, "w", encoding="utf-8") as fh:
    fh.write(ct)

# stale thumbnail of the example poster: remove part and relationship
thumb = os.path.join(BUILD, "docProps", "thumbnail.jpeg")
if os.path.exists(thumb):
    os.remove(thumb)
rels_path = os.path.join(BUILD, "_rels", ".rels")
with open(rels_path, encoding="utf-8") as fh:
    top_rels = fh.read()
top_rels = re.sub(r'<Relationship [^>]*Target="docProps/thumbnail\.jpeg"[^>]*/>', "", top_rels)
with open(rels_path, "w", encoding="utf-8") as fh:
    fh.write(top_rels)

# core properties: title
core_path = os.path.join(BUILD, "docProps", "core.xml")
with open(core_path, encoding="utf-8") as fh:
    core = fh.read()
core = re.sub(r"<dc:title>.*?</dc:title>",
              "<dc:title>%s</dc:title>" % escape("ESMO 2026 poster 4355P: Biomarker signals in MSS/pMMR mCRC are "
                                                  "predominantly prognostic, not predictive (METIMMOX)"), core)
with open(core_path, "w", encoding="utf-8") as fh:
    fh.write(core)

# zip
if os.path.exists(OUT):
    os.remove(OUT)
with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as zf:
    zf.write(ct_path, "[Content_Types].xml")
    for root, _dirs, files in os.walk(BUILD):
        for f in files:
            full = os.path.join(root, f)
            arc = os.path.relpath(full, BUILD).replace(os.sep, "/")
            if arc == "[Content_Types].xml":
                continue
            zf.write(full, arc)
print("wrote", OUT, os.path.getsize(OUT), "bytes;", len(shapes), "shapes")
print("table bottom %.2f, footnote bottom %.2f, KM %.2f-%.2f (w %.2f)" % (table_bottom, foot_bottom, KM_TOP, KM_BOTTOM, KM_W))
