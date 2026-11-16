# TÜBİTAK Report Audit — R0

**Audited files:**

- `TUBITAK_REPORT_DRAFT_R0.docx`
- `TUBITAK_SUPPLEMENTARY_MATERIAL_R0/`

**Sources, in priority order:**

1. `RESEARCH_REPORT_DRAFT_R0_3.md` (`4b52a2a`)
2. `TUBITAK_REPORT_MAPPING_R0.md` (`308fbb0`)
3. `docs/tubitak_visualization_notes.md` (`0692646`)
4. the official template `proje_yazma_sablonu_2026-2027.docx` (uploaded by the user)
5. the evaluation form `on_degerlendirme_kriterleri_2025_web.pdf` (uploaded by the user)

**What was and was not run:**

- No simulation, training, benchmark or validation script was run.
- No source, result, figure, validation record, README, R0.3 or mapping file was modified.
- The report was assembled by a script that edits a copy of the template. LibreOffice was used only to render the result
  for page counting and visual checks. The LibreOffice Math package was installed in the (ephemeral) container so that
  Word equations would render; the repository was not affected.

---

## 1. Files created

| File | Content |
|---|---|
| `TUBITAK_REPORT_DRAFT_R0.docx` | first complete report draft, built on a copy of the official template |
| `TUBITAK_REPORT_AUDIT_R0.md` | this audit |
| `TUBITAK_SUPPLEMENTARY_MATERIAL_R0/EK_BELGE_1_Cozucu_Dogrulamasi.{docx,md}` | solver validation details (R0.3 Ek A, realization convergence) |
| `TUBITAK_SUPPLEMENTARY_MATERIAL_R0/EK_BELGE_2_Ek_Tasinim_Sonuclari.{docx,md}` | full W/k/L_d tables, closed-box figure, paired-analysis figure, barrier result, illustration |
| `TUBITAK_SUPPLEMENTARY_MATERIAL_R0/EK_BELGE_3_PINN_Ayrintilari.{docx,md}` | full PINN section, training diagnostics, the three 3D visualizations (fig2 first), two further PINN figures |
| `TUBITAK_SUPPLEMENTARY_MATERIAL_R0/EK_BELGE_4_RGF_ve_WZ.{docx,md}` | RGF B1–B4 measurements and statuses, W/Z diagnostics (R0.3 §9 + Ek B) |
| `TUBITAK_SUPPLEMENTARY_MATERIAL_R0/EK_BELGE_5_Yasak_Alan_Notu.{docx,md}` | the 15% / 25% forbidden-area note (analytic estimate, not a measurement) |
| `TUBITAK_SUPPLEMENTARY_MATERIAL_R0/EK_BELGE_6_Tekrar_Uretim.{docx,md}` | reproduction commands, environment, reproduction status, metadata note, Born-ratio caption analysis |
| `TUBITAK_SUPPLEMENTARY_MATERIAL_R0/figures/` | copies of the 16 supplementary images (originals unchanged) |

**Supplementary documents:**

- The text is extracted verbatim from R0.3. Only three kinds of change were made:
  - numbered citations were converted to APA;
  - internal cross-references were remapped to the report and EK BELGE numbering;
  - "Taslak notu" boxes were removed.
- Supplementary documents keep R0.3's decimal points. The main report uses Turkish decimal commas.

**Page counts** (LibreOffice):

| Document | Pages |
|---|---|
| EK BELGE 1 | 3 |
| EK BELGE 2 | 4 |
| EK BELGE 3 | 7 |
| EK BELGE 4 | 6 |
| EK BELGE 5 | 1 |
| EK BELGE 6 | 2 |

**Build scripts.** The scripts that generated the files (`build_report.py`, `content.txt`, `refs.py`,
`build_supplements.py`) are in the session scratchpad and were **not committed**, because the task asked to commit only
the report outputs. They can be added on request.

---

## 2. Page count and formatting verification

### 2.1 Page count

**The main report is 20 pages**, rendered with LibreOffice 24.2 and Liberation Serif, which is metric-compatible with
Times New Roman. The limit is 20 pages.

| Section | Starts on page |
|---|---|
| Field, title, AI declaration | 1 |
| Özet, Amaç | 2 |
| Giriş | 3 |
| Yöntem | 4 |
| Proje İş-Zaman Çizelgesi | 8 |
| Bulgular | 9 |
| Sonuç ve Tartışma | 15 |
| Öneriler | 16 |
| Kaynaklar | 17 |
| Ekler | 19 |
| Ek 1 (abbreviations, new page) | 20 |

**The margin is tight.**

- Microsoft Word may paginate slightly differently. **The page count must be re-checked in Word.**
- If Word gives 21 pages, the first remedy is to move Ek 1 (abbreviations) into the EK BELGELER and keep the list of
  supplementary documents in "Ekler".
- Do not cut limitations to save space.

### 2.2 Formatting checks

| Check | Result |
|---|---|
| Font | Times New Roman in every text run |
| Font size | 12 pt for body text, tables, captions and references. Exceptions: the schedule table's body cells are 10 pt (so the 10 month columns fit); its header cells keep the template's own 9 pt |
| Margins | 2.5 cm on all sides, set explicitly. **The template file's own page setup is top 1.5 cm / bottom 1.25 cm**; the written rule (2.5 cm everywhere) was followed |
| Line spacing | single (1.0) for all inserted paragraphs |
| Alignment | body paragraphs justified. **Exception:** the reference list is left-aligned with a hanging indent, because justified DOI lines were badly stretched. Revert if the advisor requires strict justification |
| Red instructional text | removed. 19 instructional paragraphs were deleted, including the four "Örnek:" paragraphs and the final "NOT: Kırmızı…" note. The red explanatory runs inside the five label paragraphs (four AI items and "Anahtar Kelimeler:") were removed. Red colour on paragraph marks inherited from the template was set to black, so that text typed later in Word does not appear red. Verified: only `000000` and `auto` colours remain |
| Equations | 8 numbered Word (OMML) equations, (1)–(8) |
| Figures | 6, inserted at their native aspect ratio (width only) |
| Tables | 11 numbered tables plus the schedule table |
| Abstract | 182 words (limit 150–250) |
| Keywords | 5 (see open issue S-7) |
| Title | the exact 11-word title, unchanged |

### 2.3 Template headings

All 19 black bold labels of the template were checked **programmatically** against the template. They are present,
textually identical and in the original order:

1. Proje Ana Alanı
2. Proje Tematik Alanı
3. Proje Adı (Başlığı)
4. Yapay Zeka Araçları Kullanım Beyanı:
5. Kullanılan araçlar:
6. Kullanım amacı:
7. Kullanım kapsamı:
8. İnsan müdahalesi:
9. Özet
10. Anahtar Kelimeler:
11. Amaç
12. Giriş
13. Yöntem
14. Proje İş-Zaman Çizelgesi
15. Bulgular
16. Sonuç ve Tartışma
17. Öneriler
18. Kaynaklar
19. Ekler

The heading paragraphs' own runs and formatting are untouched; values were appended after the ":" of the three field
labels.

**No new top-level headings were added.** Inside sections, bold body-text sub-headings (not Word heading styles) were
added for readability:

| Section | Sub-headings added |
|---|---|
| Yöntem | 10 sub-headings |
| Bulgular | 6 sub-headings |
| Ekler | "Ek 1. Kısaltmalar" |

Remove them if the advisor wants plain sections.

### 2.4 Schedule table

- The template's layout is kept: "İşin Adı" column, "AYLAR" header, NİSAN–OCAK columns and the template's colours.
- The template's example row labels (Literatür Taraması, Arazi Çalışması, Verilerin Toplanması ve Analizi, Proje Raporu
  Yazımı) were replaced by the documented activities. "Arazi Çalışması" does not apply to a numerical study.
- Rows were added by copying the template's own row XML.
- Cell legend: ✓ = documented by commit dates; Ö = [ÖĞRENCİ TARAFINDAN TAMAMLANACAK]; — = no documented activity.
- The planned, not-done work (B4, the pre-registered ζ campaign) is in a row labelled "(plan)" with Ö cells only. It is
  not marked as completed.

---

## 3. Unresolved student-input placeholders

**In the main report** (8 occurrences of `[ÖĞRENCİ TARAFINDAN TAMAMLANACAK]`):

| # | Location | Needed |
|---|---|---|
| P-1 | AI declaration, *Kullanılan araçlar* | tool version as shown in the tool |
| P-2 | AI declaration, *Kullanılan araçlar* | access dates other than the documented commit range (5–16 Nov 2026) |
| P-3 | AI declaration, *Kullanılan araçlar* | any other AI tools (name, version, date), e.g. for the July–September code, translation or literature search |
| P-4 | AI declaration, *Kullanım amacı* | confirmation or correction of the purpose list |
| P-5 | AI declaration, *Kullanım kapsamı* | which parts the student designed, wrote, ran and checked; a percentage only if it can be justified |
| P-6 | AI declaration, *İnsan müdahalesi* | what the student changed, re-derived or verified; how the references were checked; the advisor's role |
| P-7 | Yöntem, Tekrar üretilebilirlik | whether the repository is public, and its licence |
| P-8 | Schedule legend and Ö cells | real months of the work; April–June 2026 activity; plans for December–January |

**Not placeholders in the text, but required before submission:**

| # | Item |
|---|---|
| P-9 | Author, school and advisor details. The template body does not ask for them; the application system may |
| P-10 | Confirmation that "Kuantum Teknolojileri" is a valid thematic field in the 2026–2027 call list (not in the provided documents) |
| P-11 | The AI commit counts (102 / 78) describe the repository **before** this report commit. After this commit the repository holds 103 non-merge commits, 79 of them by Claude. Update the counts at submission time |

**In the supplements:**

| File | Marker |
|---|---|
| EK BELGE 6 | `[ÖĞRENCİ TARAFINDAN TAMAMLANACAK: --strict ile yeniden üretim]` |
| EK BELGE 5 | `[AÇIK KONU: README'nin bu değerleri nasıl elde ettiği belgelenmelidir.]` |

---

## 4. AI usage declaration

- The declaration uses the template's exact four items.
- It states that AI was used extensively.
- It names the tool as documented in the repository: Claude, in the Claude Code environment. Commits are authored
  "Claude", with a `Co-Authored-By: Claude …` trailer.
- It gives the documented commit range and counts.
- It lists the purposes that are visible from commit content.
- It states that the first draft of **every** section of this report was prepared by AI.
- It states that the repository does not document the student's verification work, and that the report must not be
  submitted before the student and advisor review the declaration.

**Not stated, because it would be invented:** the tool version, other access dates, percentages, tasks done personally
by the student, and any claim of independent verification or understanding.

**No eligibility determination was made.** The evaluation form's ethics clause eliminates a project if the work is found
not to have been developed by the project team. The advisor must judge this; it is the main blocker (§10).

---

## 5. Scientific audit

### 5.1 Checks performed

1. **Numbers.** All 82 decimal numbers in the report body were checked automatically. Each appears verbatim (comma → dot)
   in R0.3 or in `docs/tubitak_visualization_notes.md`.
   - One rounding introduced during writing (|ΔT| = 2,0×10⁻⁵) was corrected to R0.3's 1,96×10⁻⁵.
   - Integers and exponents were checked manually against R0.3.
   - Validation-table values (4,6×10⁻¹⁴, 9,4×10⁻¹⁰, 5,2×10⁻⁴, 1,3×10⁻²) are R0.3's own rounded forms (R0.3 §6.4).
2. **Forbidden claims.** The report text was scanned; every match is a negation or a correctly hedged statement.
   - "ENAQT": 0 occurrences.
   - No 3D physical simulation is claimed; the only "3B" items are the supplementary visualizations, labelled as 2D
     fields shown as surfaces.
   - "lokalizasyon kanıtı değildir" is repeated in Özet, Bulgular and Sonuç.
   - B4 is "YAPILMADI" (H6, Tablo 3, Tablo 9, Sonuç, Öneriler).
   - W/Z is "DOĞRULANMADI"; coupling-scaling is "ÇALIŞTIRILMADI"; R0 is "BEKLEMEDE"; coupled-block is "tasarım
     yetersiz".
   - PINN appears as "ikincil", with "üstünlüğü bulunmamıştır" and "PINN'ler hakkında genel bir hüküm değildir".
3. **Main claim.** The wording is consistent with R0.3: an exploratory non-monotonic trend over seven sampled ζ values,
   with a coarse minimum at kζ ≈ 0,6–1,4 whose exact location is not determined. All four limits are stated in Bulgular
   and Sonuç; the abstract states the post-hoc and no-campaign limits.
4. **Validation language.** A (V1, V2, V6), B (V3, V4, V5, V7) and C (B2 disorder-free only; B4 not performed) are kept
   separate in Tablo 3. The Sonuç sentence follows R0.3 §12.1.
5. **Born comparison.** It is described as a theoretical approximation and an order-of-magnitude comparison, not a
   validation. All four limits (infinite-2D derivation, finite guide, 6 channels, rough diffusion expression) are in
   Yöntem.
6. **PINN statistics.** Three seeds per group, t = 2,70, p = 0,086, not significant at 5%, small sample, exploratory.
7. **Partial reproducibility.** Disclosed in Yöntem (Tekrar üretilebilirlik) and in Sonuç limitation (9).
8. **PINN 3D visualizations.** Placed only in EK BELGE 3, with: a separate 1 × 1 box, W = 56.25, ζ = 0.08, seed 0,
   T = 0.02, a 127 × 127 grid, ψ error 0.5864 vs |ψ|² error 0.6932, a single realization and a single network seed.
   The white animation is noted as exhibition-only.

### 5.2 Unresolved scientific issues

| # | Issue | Status / action |
|---|---|---|
| S-1 | **No independent stationary cross-validation of disordered transport (B4 not performed).** | disclosed; future work |
| S-2 | **The ζ result is exploratory.** Post-hoc contrasts, 7 points, no frozen campaign, ζ = 0.07 contrast not significant after correction. | disclosed |
| S-3 | The K-2 mismatch remains in the repository: `figures/CAPTIONS.md` and `make_figures.py` still say "~1.05". | the report caption uses the data ranges; the analysis is in EK BELGE 6; fix the repository files in a separate task |
| S-4 | README's 15% / 25% forbidden-area values have no documented source. | EK BELGE 5 only, labelled as an analytic estimate |
| S-5 | Partial reproduction: momentum, length, barrier and convergence sweeps and 17 PINN runs not reproduced. | disclosed; a `--strict` re-run would need a separately authorised task |
| S-6 | The research gap is phrased modestly; no targeted literature search for closer prior work was done. | student/advisor literature check recommended (indicators 2 and 7) |
| S-7 | The template says keywords should be "en fazla beş kelimeden oluşan". The report uses **5 keyword terms** (several multi-word), as instructed. If the jury reads it as 5 words in total, a single-word set is needed. | advisor decision |
| S-8 | The Born ℓ_tr minimum "kζ ≈ 0,8" comes from README §9.3 / `disorder.born_rates`. | presented as a theoretical expectation, not a finding |
| S-9 | In Bulgular, the barrier/tunnelling negative result ("düzensizlik destekli tünelleme gözlenmedi") appears only through the EK BELGE list and EK BELGE 2. | add one sentence to Bulgular if Word pagination allows |
| S-10 | R0.3 §4.5 project history is not in the report body: the template says to avoid history narration in Amaç. | the earlier artefact (the "düzensizlik iletimi artırır" claim) is mentioned in Bulgular; the AI declaration documents the timeline |

---

## 6. Reference conversion

- **All 41 references** from R0.3 §14 were converted to APA 7 author–date and listed alphabetically.
- Every in-text citation is generated from a single reference table, so the text and the list cannot diverge.
  - **Every list entry is cited in the main report** (41/41, checked by the build).
  - Every citation appears in the list.
  - The supplements cite 15 of the entries; all are in the main list.
- In-text forms:
  - "(Yazar, yıl)" for one author;
  - "(Yazar & Yazar, yıl)" for two, parenthetical;
  - "Yazar ve Yazar (yıl)" for two, narrative;
  - "Yazar vd." for three or more.
- Same-surname, same-year collisions were disambiguated: "Wang, Wang vd., 2021" and "Wang, Teng vd., 2021". The
  earlier entries Lee & Fisher 1981 and Fisher & Lee 1981 are distinct.

**Issues to resolve (no details were invented):**

| # | Issue | Entries |
|---|---|---|
| R-1 | R0.3 gives the author list only as "et al.". APA 7 requires all authors (up to 20); the list currently has "vd." | Billy 2008, Roati 2008, Rahaman 2019, Tancik 2020 |
| R-2 | Only start pages are known; APA expects page ranges | almost all journal articles (exception: Holm 1979, Welch 1947) |
| R-3 | Journal abbreviations were expanded to full titles. These are standard expansions; spot-check recommended | e.g. "Phys. Rev." → *Physical Review* |
| R-4 | **Welch (1947) has no DOI**, marked "[DOI doğrulanmadı]" in the list. Remove the marker or add a verified DOI before submission | [39] |
| R-5 | Identifiers kept as in R0.3; arXiv URLs were formed from the arXiv IDs given in R0.3 | Holm (JSTOR 4615733), North (PMID 12111669) |
| R-6 | Conference and workshop items use APA-style approximations from R0.3's data; verify proceedings details | ICML / NeurIPS / ICLR / NeurIPS workshop entries |
| R-7 | Software (Python, NumPy, SciPy, PyTorch) is not cited. Adding it would require verified entries | — |

---

## 7. Figures used and omitted

### 7.1 Main report

All files exist in the repository and were inserted at their native aspect ratio.

| No | File | Width | Source of data (as in caption) |
|---|---|---|---|
| Şekil 1 | `figures/disorder_examples.png` | 9.5 cm | generated deterministically by `disorder.py` |
| Şekil 2 | `figures/transport_vs_correlation_length.png` | 16 cm | `results/transport/corr_length/` |
| Şekil 3 | `figures/born_ratio_vs_kl.png` | 10 cm | `results/transport/analysis.json` (caption uses the data ranges, not `CAPTIONS.md`) |
| Şekil 4 | `figures/transport_vs_strength.png` | 16 cm | `results/transport/strength/` |
| Şekil 5 | `figures/transport_vs_length.png` | 16 cm | `results/transport/length/`, `analysis.json` |
| Şekil 6 | `figures/pinn_errors.png` | 11.5 cm | `results/pinn/summary.csv` |

`figures/localization_closed_box.png` was **moved to EK BELGE 2** because of the page budget. Its data table (Tablo 8)
stays in the main report.

### 7.2 Supplements

| EK BELGE | Figures |
|---|---|
| EK BELGE 1 | `validation_solver.png`, `realization_convergence.png` |
| EK BELGE 2 | `transport_vs_momentum.png`, `localization_closed_box.png`, `paired_corr_length.png`, `transport_barrier.png`, `waveguide_snapshots.png` (illustration; no stored data) |
| EK BELGE 3 | `visualizations/tubitak/fig2_pinn_dst_karsilastirma_3b.png` (main supplementary PINN comparison), `fig3_pinn_dst_mutlak_fark_3b.png`, `fig1_duzensizlik_potansiyeli_3b.png`, `figures/pinn_error_vs_time.png`, `figures/pinn_vs_reference_W56.png` |
| EK BELGE 4 | `results/rgf_b2_validation/error_vs_dx.png`, `docs/turkce_ozet_figurler/sekil1_b1_faz.png`, `sekil2_63kanal.png`, `sekil3_coupled.png` |

**Omitted:**

- the white-background animation (exhibition only);
- `results/rgf_b2_validation/error_vs_energy.png` and `order_vs_energy.png` (available in the repository).

### 7.3 Figures that still need Turkish labels

These have English axis labels, titles or legends. They were **not altered**; their captions say "Şekil etiketleri
İngilizcedir".

- `figures/disorder_examples.png`
- `figures/transport_vs_correlation_length.png`
- `figures/born_ratio_vs_kl.png`
- `figures/transport_vs_strength.png`
- `figures/transport_vs_length.png`
- `figures/pinn_errors.png`
- `figures/localization_closed_box.png`
- `figures/paired_corr_length.png`
- `figures/transport_vs_momentum.png`
- `figures/transport_barrier.png`
- `figures/validation_solver.png`
- `figures/realization_convergence.png`
- `figures/pinn_error_vs_time.png`
- `figures/pinn_vs_reference_W56.png`
- `figures/waveguide_snapshots.png`
- `results/rgf_b2_validation/error_vs_dx.png`

Regenerating them with Turkish labels is a plotting-only change from stored data and needs separate authorisation. The
`visualizations/tubitak/*.png` figures are already Turkish.

### 7.4 Readability at print size

`transport_vs_strength.png` (three panels, natively 12.9 in wide) and `transport_vs_length.png` are scaled to 16 cm.
Their axis text prints at roughly 5–6 pt, which is small. Regenerating them with larger fonts is the fix, under the same
authorisation as §7.3.

---

## 8. Evaluation-criteria weaknesses (12 indicators)

| # | Indicator | Remaining weakness in R0 |
|---|---|---|
| 1 | main-field fit | the physics fit is clear; the link to "Kuantum Teknolojileri" is one modest sentence (P-10) |
| 2 | problem and question | stated in Amaç; the gap is modest and not backed by a targeted search for direct prior work (S-6) |
| 3 | objectives linked to the problem | Tablo 1 links H1–H6 to the questions |
| 4 | measurable and achievable objectives | measurable outputs are given; H6 is not achieved. The objectives were defined during the work; the report says they are not pre-registered |
| 5 | method suitable and sufficient | weakened by B4 not being performed and the exploratory statistics (S-1, S-2); both are disclosed |
| 6 | method explained | detailed. Implementation specifics (PINN, RGF) are in EK BELGELER; the jury must open them for full detail |
| 7 | innovative approach | the novelty is modest (a numerical test of a known prediction in a finite few-channel geometry); stated without exaggeration |
| 8 | potential impact | limited to the code/method base and a possible comparison basis for experiments; no technology claims |
| 9 | clarity, correctness, reproducibility | partial reproduction (S-5); English figure labels (§7.3); small text in two figures (§7.4); reference details (§6) |
| 10 | evidence supports the objectives | H1–H5 are supported; H6 is not; the main claim is exploratory |
| 11 | comparison with standard procedures | spectral split-step (standard) vs stationary RGF (incomplete) vs PINN (worse and slower) is discussed in Sonuç; there is no comparison table, for space reasons |
| 12 | can be developed and disseminated | depends on repository visibility and licence (P-7) |
| Ethics | — | **AI authorship (§4) is the decisive open item** |

---

## 9. Content not transferred, or only condensed

| R0.3 part | Treatment |
|---|---|
| §4.5 project history | not in the body (template: no history in Amaç). The timeline is reflected in the AI declaration and schedule |
| §5.1 symbol table | omitted for space; symbols are defined in the text and in Ek 1 |
| §6.1 layout list | condensed into a sentence and Tablo 2 |
| §7.2 strong-scattering ratio table | given in prose |
| §7.3 ⟨R⟩/P_in columns, §7.4 table, §7.5 raw table | condensed in the main report; full tables in EK BELGE 2 |
| §7.7 barrier result and illustrations | EK BELGE 2 only (S-9) |
| §8 PINN | condensed; the full text is in EK BELGE 3 |
| §9 RGF | status table only; details in EK BELGE 4 |
| §10.2 observation / expectation / explanation / untested table | condensed into prose in Bulgular and Sonuç; no row was dropped |
| §10.4 repository-document inconsistency table (K-1…K-10) | not transferred (internal audit notes); K-2's scientific content is in Şekil 3 and EK BELGE 6 |
| §11 limitations | all 11 items kept, condensed into one numbered paragraph in Sonuç ve Tartışma |
| §13 reproducibility | one paragraph in Yöntem; commands in EK BELGE 6 |
| Ek E abbreviations | Ek 1 in the report |
| Ek F, Ek G and all "Taslak notu" boxes | internal; not transferred |

---

## 10. Remaining blockers before submission

1. **AI declaration and authorship (§4).** The student must fill in P-1 to P-6. The advisor must review the declaration
   and decide on eligibility under the ethics clause. No eligibility determination was made here.
2. **Word page count (§2.1).** Confirm ≤ 20 pages in Microsoft Word.
3. **Placeholders and confirmations:**
   - P-7: repository visibility and licence;
   - P-8: schedule months;
   - P-10: thematic field name;
   - P-11: commit counts.
4. **Reference completion:**
   - R-1: full author lists for 4 entries;
   - R-2: page ranges;
   - R-4: Welch DOI marker.
5. **Figure language and readability decision (§7.3, §7.4).** Needs separate authorisation to regenerate the plots.
6. **Optional:**
   - S-6: targeted literature check;
   - S-9: barrier sentence;
   - S-7: keyword interpretation;
   - S-3: repository caption fix;
   - S-5: `--strict` reproduction.
