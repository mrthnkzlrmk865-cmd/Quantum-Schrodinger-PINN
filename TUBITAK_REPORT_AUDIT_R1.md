# TÜBİTAK Report Audit — R1

| Item | Value |
|---|---|
| Report | `TUBITAK_REPORT_DRAFT_R1.docx` (new file; R0 is unchanged) |
| Base | `TUBITAK_REPORT_DRAFT_R0.docx` / `TUBITAK_REPORT_AUDIT_R0.md` (`fc26acb`) |
| Scientific source of truth | `RESEARCH_REPORT_DRAFT_R0_3.md` (`4b52a2a`, unchanged) |

**What was run and what was not.**

- **Not run:** simulations, PINN training or evaluation, benchmarks, validation scripts, numerical experiments.
- **Not modified:** source code, datasets, results, validation records, README, original figures, R0, R0.3 or the
  mapping file.
- **Run:**
  - plotting code that reads stored result files (`visualizations/tubitak_r1/make_figures_tr.py`);
  - read-only checks of `results/transport/analysis.json`;
  - reference look-ups on Crossref, PMLR, NeurIPS proceedings and the arXiv API. Network access to these four hosts
    was added to the environment for this purpose, with the user's approval.

---

## 1. Files created

| File | Content |
|---|---|
| `TUBITAK_REPORT_DRAFT_R1.docx` | revised report, built on a fresh copy of the official template |
| `TUBITAK_REPORT_AUDIT_R1.md` | this audit |
| `visualizations/tubitak_r1/make_figures_tr.py` | plotting script for the six Turkish main-report figures (reads stored data only) |
| `visualizations/tubitak_r1/*_tr.png` | six regenerated figures, 300 dpi, sized to their printed width |
| `TUBITAK_SUPPLEMENTARY_MATERIAL_R1/EK_BELGE_7_Kisaltmalar.{md,docx}` | the abbreviation list, moved out of the main report (formerly "Ek 1") |
| `TUBITAK_SUPPLEMENTARY_MATERIAL_R1/README.md` | which EK BELGE lives where |

EK BELGE 1–6 are unchanged in `TUBITAK_SUPPLEMENTARY_MATERIAL_R0/`. Their references to main-report table and figure
numbers remain correct, because R1 keeps Tablo 1–10 and Şekil 1–6.

---

## 2. Page count and formatting

### 2.1 Page count

**20 pages**, rendered with LibreOffice 24.2. Times New Roman is replaced by the metric-compatible Liberation Serif.

| Section | Starts on page |
|---|---|
| AI declaration | 1 |
| Özet, Amaç | 2 |
| Giriş | 3 |
| Yöntem | 4 |
| Proje İş-Zaman Çizelgesi | 8 |
| Bulgular | 9 |
| Sonuç ve Tartışma | 16 |
| Öneriler, Kaynaklar | 17 |
| Ekler | 20 |

Page 20 is about 60% full, so the LibreOffice slack is about 0.4 page.

**LibreOffice pagination is not identical to Microsoft Word pagination.** The page count must be checked in Word before
submission. If Word exceeds 20 pages, shorten the Giriş literature paragraphs first; do not cut limitations.

### 2.2 Formatting checks (on the DOCX)

| Check | Result |
|---|---|
| Margins | 2.5 cm on all four sides |
| Font | Times New Roman in every text run |
| Body text, tables, captions, references | 12 pt |
| Schedule table | label column 9 pt, month cells 10 pt; the template's own header cells keep 9 pt |
| Spacing and alignment | single line spacing; body paragraphs justified. The reference list is left-aligned with a hanging indent, as in R0; revert if strict justification is required |
| Colours | only `000000` / `auto`; no red instruction text and no red paragraph marks |
| Content | 6 figures, 10 numbered tables plus the schedule, 8 Word equations |
| Abstract | 182 words |
| Keywords | 5 |
| Ek 1 | removed from the main report; now EK BELGE 7, listed under "Ekler" |

### 2.3 Template headings

All 19 template labels were checked programmatically against the template: the label runs are identical in text, bold
and size, and in the original order. The labels are Proje Ana Alanı, Proje Tematik Alanı, Proje Adı (Başlığı), Yapay
Zeka Araçları Kullanım Beyanı:, Kullanılan araçlar:, Kullanım amacı:, Kullanım kapsamı:, İnsan müdahalesi:, Özet,
Anahtar Kelimeler:, Amaç, Giriş, Yöntem, Proje İş-Zaman Çizelgesi, Bulgular, Sonuç ve Tartışma, Öneriler, Kaynaklar,
Ekler.

No top-level heading was added. The bold body-text sub-headings inside Yöntem and Bulgular are as in R0.

### 2.4 Schedule table

The template's structure is kept: 1 label column, 10 month columns, the two header rows and the colours. Two layout
changes were made:

- **Month labels** now run HAZİRAN–MART, so that the schedule starts in June 2026 as the student stated. They were
  NİSAN–OCAK in the template.
- **Column widths:** the label column was widened to 4.4 cm, and the month columns were narrowed equally so that the
  total width is unchanged. This keeps the activity labels readable.

---

## 3. Student-provided information and how it is represented

| Student statement | Where it appears | Wording constraint applied |
|---|---|---|
| Project began in June 2026 | AI declaration (A); schedule (first column HAZİRAN; legend) | — |
| First activities: literature review and the first 3D Plotly visualization code | AI declaration (A); schedule rows 1–2 marked ✓ in June | — |
| The initial physical research question was defined by the student | AI declaration (A); schedule row 3 (month: Ö) | the declaration also states that the initial question was **not** originated by Claude |
| In the last two months the detailed scope, model and methodology evolved with substantial assistance from Claude | AI declaration (*Kullanım amacı*, B); schedule row 4 (months: Ö) | the report does **not** claim that the student independently developed later models, methods, code or analyses |

**Dates.** No activity date was inferred from Git commits. Commit dates are mentioned only as "Claude commits are dated
November 2026", explicitly not as the full period of use.

---

## 4. Unresolved placeholders

**In the report** (8 occurrences of `[ÖĞRENCİ TARAFINDAN TAMAMLANACAK]`, plus the Ö cells of the schedule):

| # | Location | Needed |
|---|---|---|
| P-1 | Kullanılan araçlar | Claude version as shown in the tool |
| P-2 | Kullanılan araçlar | access dates |
| P-3 | Kullanılan araçlar | other AI tools (name, version, date), if any |
| P-4 | Kullanım amacı | confirmation or correction of the purpose list |
| P-5 | Kullanım kapsamı (A) | the student's contributions beyond the four confirmed facts |
| P-6 | İnsan müdahalesi (C) | checks actually performed by the student and advisor, and changes made to AI output |
| P-7 | Yöntem, Tekrar üretilebilirlik | repository visibility and licence |
| P-8 | Schedule legend and Ö cells | months of all activities after June; whether B4 and the ζ campaign are planned |

**Required but not marked in the text:**

| # | Item |
|---|---|
| P-9 | Commit counts. The declaration says "103 non-merge commits, 79 by Claude", which describes the repository **before** the R1 commit. After the R1 commit the counts are 104 and 80. Update them at submission |
| P-10 | Confirmation that "Kuantum Teknolojileri" is a valid thematic field in the 2026–2027 call |
| P-11 | Author, school and advisor details, if the application system requires them |

**In the supplements** (unchanged from R0):

- EK BELGE 6: `[ÖĞRENCİ TARAFINDAN TAMAMLANACAK: --strict …]`
- EK BELGE 5: `[AÇIK KONU: …]`

---

## 5. AI usage declaration (R1)

The four template items are kept. The content now separates:

- **(A) Student contribution.** Only the four student-confirmed facts; everything else is a placeholder.
- **(B) Claude contribution.** Substantial assistance in the last two months: scope, model, method, coding, analysis
  and report drafting. The declaration adds that every section of this report was first drafted by Claude, and states
  explicitly that commit authorship alone does not show who conceived a task.
- **(C) Verification.**
  - The report does **not** claim independent human verification of derivations, code or results.
  - It lists the checks actually made while drafting: an automatic number comparison with R0.3, and reference look-ups
    on Crossref and publisher/conference records. It states that these checks were also made by the AI.
  - It lists what still needs student and advisor review: equations (1)–(8), the solver and analysis code, the
    statistical method, and the declaration itself.

**Not stated, because it would be invented:** the Claude version, access dates, other AI tools, percentages, and any
verification by the student.

**No eligibility determination was made.**

---

## 6. Born comparison paragraph

**What R0 said.** The R0 sentence "Zayıf saçılmada (kℓ_tr ≥ 10; ζ, W ve k taramalarının bütün noktaları) …" was
ambiguous, and the regime list did not read as one coherent paragraph.

**Rewrite.** The paragraph was rewritten in four parts: weak scattering; rough vs smooth disorder; stronger scattering;
and the order-of-magnitude caveat.

**Verification.** Every value was recomputed from `results/transport/analysis.json`, the source of R0.3 §7.2, using
`ratio_avg` and k = (k_x0² + π²)^{1/2}:

| Statement in R1 | Source values |
|---|---|
| weak scattering (kℓ_tr ≥ 10): 14 points (13 distinct settings, since W = 100 / ζ = 0.05 appears in both the W and ζ sweeps), ratio 0,950–1,121 | min 0.950 (ζ = 0.10 and 0.20), max 1.121 (k_x0 = 20) |
| rough disorder, kζ < 2: 11 points, kζ = 0,40–1,51, ratio 1,022–1,121 (estimate low by 2–12%) | kζ max 1.51 at k_x0 = 30 |
| smooth disorder, kζ ≥ 2: 3 points, kζ = 2,02–4,05, ratio 0,950–0,965 (estimate high by 3,5–5%) | 0.950, 0.965, 0.950 |
| stronger scattering, kℓ_tr < 10: 5 points, ratio 1,134–1,309; for kℓ_tr ≈ 1–3, 1,213–1,309 (W = 200, k_x0 = 12, W = 300) | 1.134 (W = 150), 1.172 (k_x0 = 16), 1.213, 1.303, 1.309 |
| monotonic increase as kℓ_tr decreases, within the W sweep and within the k sweep | W: 1.022 → 1.309; k: 1.062 → 1.303 |
| energy averaging: W = 300 ratio 1,479 (not averaged) → 1,309 (averaged) | `ratio` 1.479, `ratio_avg` 1.309 |

**Discrepancy found in R0.3 (not modified).** R0.3's Şekil 7.2 caption and its §7.2 text define rough disorder as
"kζ ≲ 1.4". The rough subset actually extends to kζ = 1.51 (the k_x0 = 30 point). R0.3's own K-2 table uses "kζ < 2",
which is correct. R1 uses kζ < 2 with the measured range 0,40–1,51.

**Values not stated verbatim in R0.3.** 1,134 and 1,479 are new in the report text; both are stored values in
`analysis.json`. 1,51 is computed from the stored k_x0 and ζ.

---

## 7. Reference verification (APA 7)

All 41 references were audited. Every entry is cited in the main report and every citation has an entry (checked by the
build). In-text author–year forms are unchanged from R0; the disambiguation "Wang, Wang vd." / "Wang, Teng vd." is
kept.

| Status | Count | References |
|---|---|---|
| **Verified on Crossref** (DOI resolves; authors, title, journal, volume, issue and pages taken from the Crossref record) | 31 | Abrahams; Anderson; Beenakker; **Billy** (all 10 authors; 453(7197), 891–894); Evers; Feit; Fisher & Lee; Jahnke; Karniadakis; Kosloff; Kramer; Kuhn (article 161); Lagendijk; Lee & Fisher; Lee & Ramakrishnan; Liu; MacKinnon; Muga; Phipson; Raissi; **Roati** (all 9 authors; 453(7197), 895–898); Sanchez-Palencia; Schwartz; Shinozuka; Strang; Wang/Sankaran; Wang/Teng; Wang/Wang; Wu; **Welch**; **North** |
| **Verified on proceedings pages** | 3 | **Rahaman** (PMLR vol. 97, pp. 5301–5310, all 8 authors); **Tancik** (NeurIPS 33, pp. 7537–7547, all 9 authors); Krishnapriyan (NeurIPS 34, pp. 26548–26560) |
| **Verified on arXiv** (authors, title, identifier) | 2 | Kingma & Ba (the arXiv comment confirms ICLR 2015, San Diego); Shah et al. (the workshop venue is taken from R0.3 and was **not** verified on a venue page) |
| **Not verifiable online in this session** (kept as in R0.3) | 5 | Akkermans & Montambaux (2007), Efron & Tibshirani (1993), Davison & Hinkley (1997), Datta (1995): books with no DOI. Holm (1979): not in Crossref, JSTOR not reachable; details from R0.3, URL formed from R0.3's JSTOR id 4615733 |

### 7.1 Specific resolutions

**Billy, Roati, Rahaman, Tancik.** The "vd." placeholders were replaced by complete author lists from the primary
records.

**Welch (1947).**

- Crossref holds two records for the same article:
  - Oxford University Press: https://doi.org/10.1093/biomet/34.1-2.28, *Biometrika* 34(1–2), 28–35;
  - JSTOR: 10.2307/2332510.
- The publisher DOI is used. The "[DOI doğrulanmadı]" note was removed, because a complete, verified entry now exists.

**North et al. (2002).**

- DOI 10.1086/341527 was verified, with issue 2 and pages 439–441.
- Crossref also lists a 2003 record with the same title and authors (*AJHG* 72(2), 498–499, 10.1086/346173). It is not
  cited; it should be checked whether it is an erratum relevant to the cited note.

**arXiv DOIs (Kingma, Shah).** The 10.48550/arXiv.… DOIs could not be checked, because DataCite was not reachable. The
verified arXiv abstract URLs are used instead.

**Details taken from R0.3 that Crossref does not contain:**

- the full initials of Muga et al. (Crossref gives single initials);
- "Jr." for Fleck;
- "Article 39" for Phipson & Smyth (Crossref has no page or article number).

**Journal titles.** All are full titles, matching the Crossref container titles.

---

## 8. Figures

### 8.1 Regenerated with Turkish labels (main report)

| Şekil | New file | Original (unchanged) | Changes |
|---|---|---|---|
| 1 | `visualizations/tubitak_r1/disorder_examples_tr.png` | `figures/disorder_examples.png` | Turkish titles and colour-bar label, decimal commas, 11 cm canvas |
| 2 | `…/transport_vs_correlation_length_tr.png` | `figures/transport_vs_correlation_length.png` | Turkish labels; kζ ticks 0,4/0,6/1/2/3/4 instead of 4×10⁻¹ notation; legend moved so it covers no data |
| 3 | `…/born_ratio_vs_kl_tr.png` | `figures/born_ratio_vs_kl.png` | Turkish labels; same error bars (1,96 × T_se / T_est_avg) |
| 4 | `…/transport_vs_strength_tr.png` | `figures/transport_vs_strength.png` | Turkish labels; three panels rearranged as (a) on top and (b), (c) below, for print readability; histogram y-axis headroom added for the legend; legends repositioned |
| 5 | `…/transport_vs_length_tr.png` | `figures/transport_vs_length.png` | Turkish labels ("Ohmic fit" → "difüzif fit"); same fits and parameters, read from `analysis.json` |
| 6 | `…/pinn_errors_tr.png` | `figures/pinn_errors.png` | Turkish labels; legend placed above the axes; run identifiers kept as data labels and explained in the caption |

**What is identical to the originals.** Data files, data selection, plotted quantities, error bars, fits, colours and
markers are identical: the same code paths as `make_figures.py`, reading the same stored files.

- The disorder field (Şekil 1) and the Born mean free paths (Şekil 2 b, Şekil 4 a) are evaluated with `disorder.py`,
  exactly as in the original script. This is deterministic function evaluation, not simulation.
- No smoothing, interpolation or new data.

**Readability at print size.** Each canvas is drawn at the width at which it is printed (10–16 cm). Text is 7.5–8 pt
(legends 6.2–6.8 pt) on paper; in R0, the multi-panel figures printed at roughly 5–6 pt.

**Aspect ratio.** Images are inserted without distortion (width only). The canvas proportions differ from the originals
because of the print-size redesign; Şekil 4 also changed from a 1×3 to a 2-row panel layout. The panel content is
unchanged.

**Font.** DejaVu Sans is used inside the figures, which supports all Turkish characters; the report text is Times New
Roman.

### 8.2 Figures still in English (supplements only; not regenerated)

| EK BELGE | Figures |
|---|---|
| EK BELGE 1 | `validation_solver.png`, `realization_convergence.png` |
| EK BELGE 2 | `transport_vs_momentum.png`, `localization_closed_box.png`, `paired_corr_length.png`, `transport_barrier.png`, `waveguide_snapshots.png` |
| EK BELGE 3 | `pinn_error_vs_time.png`, `pinn_vs_reference_W56.png` |
| EK BELGE 4 | `results/rgf_b2_validation/error_vs_dx.png` |

These can be regenerated the same way. `waveguide_snapshots.png` would require rerunning the solver (it has no stored
data), so it must stay as is.

**Already Turkish:**

- the three 3D PINN figures (EK BELGE 3), which remain there and are not in the main results;
- the `docs/turkce_ozet_figurler` figures (EK BELGE 4).

---

## 9. Scientific consistency audit (R1 vs R0.3)

| Check | Result |
|---|---|
| Numbers | all decimal values in the report body appear in R0.3 or the visualization notes, except 1,134 and 1,479, which come from `analysis.json` (§6). No unsourced number |
| Correlation-length claim | unchanged: an exploratory non-monotonic trend over seven sampled points, a coarse minimum at kζ ≈ 0,6–1,4, the exact location not determined; all four limits stated (seven points, post-hoc contrasts, no frozen campaign, ζ = 0,07 contrast not significant after correction) |
| Statistics | C1/C2/C3/C3b, Holm-adjusted p = 1,5×10⁻⁴ at bootstrap resolution, argmin shares 98,9% / 1,1%, vertex kζ = 0,95 [0,88; 1,12] (post-hoc) — unchanged |
| Born comparison | rewritten (§6); described as an order-of-magnitude comparison, not a quantitative validation; Yöntem keeps all four limits (infinite 2D, finite guide, 6 channels, rough diffusion expression) |
| PINN | secondary; errors unchanged (clean 0,489/0,487/0,567; W = 56: 0,586/0,574/0,610; etc.); Welch test: 3 seeds per group, t = 2,70, p = 0,086, not significant, exploratory |
| Validation status | A/B/C separation kept; **B4 NOT PERFORMED** ("YAPILMADI"); **W/Z NOT VALIDATED** ("DOĞRULANMADI"); **coupling-scaling NOT RUN** ("ÇALIŞTIRILMADI"); coupled-block design insufficient; R0 review pending; no independent stationary validation of disordered transport |
| Localization | "lokalizasyon kanıtı değildir" (4 places); "Anderson lokalizasyonu gösterilmemiştir"; closed-box negative result kept |
| Reproducibility | partial reproduction disclosed (Yöntem; Sonuç, limitation 9) |
| Negative results | closed box, PINN shortfall and the "düzensizlik iletimi artırır" artefact are kept; the barrier result is still only in EK BELGE 2 (open item S-9 from R0) |
| Outdated claims | ENAQT: 0 occurrences; no 3D physical simulation claimed; "3B" appears only as "3B Plotly görselleştirme" and the EK BELGE 3 visualizations |

---

## 10. Remaining student / advisor decisions

1. Fill in P-1 to P-8 and review the AI declaration together.
   - **No eligibility determination was made.** The evaluation form's ethics clause concerns work not developed by the
     project team.
2. Check the **page count in Microsoft Word** (≤ 20).
3. Update the commit counts at submission (P-9); confirm the thematic field (P-10).
4. Accept or revert the two schedule changes (month relabelling, wider label column).
5. Accept the Şekil 4 two-row layout, or request a 1×3 version.
6. Decide whether to regenerate the English supplementary figures (§8.2).
7. Check the uncited 2003 North et al. record, and the five references not verifiable online (§7).
8. Carried over from R0:
   - the keyword interpretation ("en fazla beş kelime");
   - an optional barrier sentence in Bulgular (S-9);
   - a targeted literature check for closer prior work (S-6);
   - the repository caption fix in `figures/CAPTIONS.md` (K-2), which should also correct R0.3's "kζ ≲ 1.4" wording (§6).
