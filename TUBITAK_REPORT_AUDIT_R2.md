# TÜBİTAK Report Audit — R2

**Files.** `TUBITAK_REPORT_DRAFT_R2.docx` (new). Base: `TUBITAK_REPORT_DRAFT_R1.docx` / `TUBITAK_REPORT_AUDIT_R1.md`
(`f670d42`). Scientific source of truth: `RESEARCH_REPORT_DRAFT_R0_3.md`.

**Not modified:** R0, R1, R0.3, source code, datasets, original figures, existing supplementary documents.

**Not run:** simulations, training, benchmarks, validation scripts.

**Network.** For the Shah et al. venue check, `ml4physicalsciences.github.io` was added, with the user's approval, to the
environment's network allowlist. It joins the four reference hosts added in R1.

---

## 1. Page count and formatting

### 1.1 Page count

**20 pages**, rendered with LibreOffice 24.2 (Liberation Serif, metric-compatible with Times New Roman).

| Section | Starts on page |
|---|---|
| Özet, Amaç | 2 |
| Giriş | 3 |
| Yöntem | 4 |
| Proje İş-Zaman Çizelgesi | 8 |
| Bulgular | 9 |
| Sonuç ve Tartışma | 16 |
| Öneriler | 17 |
| Kaynaklar | 18 |
| Ekler | 20 |

Page 20 is about 45% full. **Word pagination may differ; confirm the count in Microsoft Word.**

All 20 pages were rendered and inspected, with particular attention to pages 3–4, 8–10, 14–16 and 17–20. No overlapping
elements, orphaned captions or broken tables were found.

### 1.2 Formatting checks

| Check | Result |
|---|---|
| Template labels | **all 19 identical** to the template (text, bold, size) and in the original order |
| Margins | 2.5 cm on all sides |
| Font | Times New Roman in every text run |
| Font size | body, tables, captions and references 12 pt; schedule cells 9–10 pt and template header cells 9 pt, as in R1 |
| Colours | only black or auto |
| Spacing and alignment | single spacing; body paragraphs justified (reference list left-aligned with a hanging indent, as in R0/R1) |
| Content | 6 figures, 8 equations, 10 numbered tables plus the schedule |
| Abstract | **194 words** (R1: 182; the PINN sentence is longer, see §5) |
| Keywords | 5 |

---

## 2. Introduction layout (pages 3–4)

**Inspection.** The rendered PDF pages, not only the extracted text, were inspected at 100 dpi.

**Finding.** Sentences and paragraphs were **not** reordered in R1. The problems were of two other kinds.

1. **Line breaks inside short mathematical expressions.**
   - "kζ" was separated from "≪ 1" and from "≫ 1" at line ends: twice on page 4, in the correlation-length paragraph and
     in the Born paragraph.
   - "kℓ" was separated from "yaklaşık 1" on page 3, in the 2D localization-length paragraph.
   - "√2" was separated from "sin(nπy)" on page 5.
   - In text extraction these breaks look like interrupted sentences.
2. **Extraction artefact.** `pdftotext` places the glyphs of display equations out of order; for example, equation (2)
   is dumped above and below its line. This is not a visual problem: in the rendered page each equation follows its
   introducing sentence.
3. **Forward equation reference.** The Born paragraph in Giriş referred to "(denklem 5)", which only appears later, in
   Yöntem.

**Fixes (no scientific rewriting).**

- Non-breaking spaces are now inserted on both sides of the relation operators ≪ ≫ ≈ = < > ≤ ≥ ∼, in "√2 sin" and in
  "kℓ yaklaşık 1". Short expressions can no longer wrap mid-expression.
- The reference now reads "(Yöntem bölümündeki denklem 5)".

**Result.** Pages 3–4 were re-rendered after the fix. All four checked items now read continuously:

- the 2D localization-length paragraph;
- the kζ ≪ 1 / kζ ≫ 1 explanation;
- the Born paragraph;
- the transition into Yöntem with equation (1).

Equation placement was left as it was, because the equations themselves did not cause the problems.

---

## 3. Born comparison wording

- "W ve k taramalarının her birinde oran, kℓ_tr azaldıkça monoton olarak artar." was replaced, as requested, by
  "İncelenen W ve k taramalarındaki örneklenmiş noktalarda oran, kℓ_tr azaldıkça artan bir eğilim göstermiştir."
- Two neighbouring phrases that could read as a general law were also limited to the sampled points:
  - Bulgular: "oran 1,134–1,309'a çıkar" became "oran 1,134–1,309 aralığındadır";
  - Sonuç: "örneklenmiş noktalarda oran 1,134–1,309 aralığındadır".
- All values are unchanged (0,950–1,121; 1,022–1,121; 0,950–0,965; 1,134–1,309; 1,213–1,309; 1,479 → 1,309; kζ
  0,40–1,51 and 2,02–4,05). The order-of-magnitude caveat ("mertebe kıyası … niceliksel bir doğrulaması değildir") is
  kept.

---

## 4. Schedule

| Row | Cells | Basis |
|---|---|---|
| June | ✓ for literature review and the first 3D Plotly visualization code | student-confirmed |
| Literature review | label changed from "(sürüyor)" to "(Haziran'da başladı)" | whether it is still ongoing is not confirmed |
| Initial research question defined by the student | **Ö placeholders** (Haziran–Ağustos) | — |
| Later activities | Ö placeholders | — |
| Rapor yazımı | D (ongoing) in November | the drafting work currently being done |
| B4 / ζ campaign | Ö (yapılmadı; planned only if the student says so) | nothing planned is marked completed |

- **Research-question row.** The legend now states that the month in which the question was defined is not yet
  confirmed. The June–August window reflects only the student's statement that this question preceded the
  Claude-assisted refinement of the last two months. It is not inferred from Git. If the advisor prefers, it can be
  reduced to a single Ö cell.
- **Dates.** No date is taken from Git history.

---

## 5. PINN wording (limited to the actual experiment)

| Location | R1 | R2 |
|---|---|---|
| Özet | "… mevcut hesaplama bütçesinde klasik çözücüye yakın doğruluğa ulaşamamıştır." | "İkincil bir karşılaştırma olarak 1 × 1 kapalı kutu başlangıç değer probleminde denenen fizik-bilgili sinir ağının doğruluğu, kullanılan mimari ve hesaplama bütçesi altında klasik çözücüye göre yetersiz kalmıştır." |
| Bulgular | "yani PINN bu bütçede taşınım fiziği için yeterince doğru değildir." | "İncelenen kapalı kutu başlangıç değer probleminde, kullanılan PINN mimarisi ve hesaplama bütçesi altında elde edilen doğruluk, klasik çözücüyle karşılaştırıldığında yetersiz kalmıştır." |
| Bulgular | "Bu ileri problem için PINN'in … bir üstünlüğü bulunmamıştır." | "Bu problemde ve bu bütçede … doğruluk veya hız açısından bir üstünlüğü gözlenmemiştir; PINN ana deneylerin dalga kılavuzu geometrisinde sınanmamıştır." |
| Sonuç (comparison) | "PINN ise bu ileri problemde … daha yavaştır." | "İncelenen kapalı kutu başlangıç değer probleminde PINN … daha yavaş olmuştur." |
| Sonuç (scope) | "PINN'ler hakkında genel bir hüküm değildir." | "bütün PINN'ler veya bütün kuantum taşınımı problemleri hakkında genel bir hüküm değildir." |
| Sonuç (contribution) | "PINN'in bu bütçede yetersiz kaldığının nicel gösterimi" | limited to "incelenen kapalı kutu probleminde, kullanılan PINN mimarisi ve bütçesiyle…" |
| Sonuç (applicability) | — | "kapalı kutu problemindeki PINN denemesi için olumsuz … sonuç" |

Errors, run times (0,7–2,9 s vs 470–3100 s), the Welch result (3 seeds per group, p = 0,086, not significant) and all
PINN limitations are unchanged.

---

## 6. References

### 6.1 North et al. (2002) and the 2003 same-title record

**Evidence from Crossref:**

| DOI | Record |
|---|---|
| 10.1086/341527 | *AJHG* 71(2), 439–441, August 2002 |
| 10.1086/346173 | *AJHG* 72(2), 498–499, February 2003; same title and authors; its reference list **cites 10.1086/341527** |

- Neither record has an `update-to`, `updated-by` or erratum relation.
- The publisher page (Elsevier) is not reachable from this environment.

**Conclusion.**

- The 2003 item is a **separate record that refers back to the 2002 note**. Whether it is an erratum, a correction or a
  follow-up letter **could not be determined** from the available records.
- R0.3 intends the 2002 article: it gives vol. 71, p. 439 and PMID 12111669.
- **The citation is unchanged** (North et al., 2002; 71(2), 439–441; DOI 10.1086/341527).
- *Open:* check the 2003 item on the journal site; if it is an erratum, APA allows noting it, but no change is needed
  for the cited content.

### 6.2 Shah et al. (2022) workshop venue

**Verified on the official workshop page** (https://ml4physicalsciences.github.io/2022/):

- The page title reads "Machine Learning and the Physical Sciences, NeurIPS 2022".
- The workshop was held on 3 December 2022 at the New Orleans Ernest N. Morial Convention Center.
- The accepted-papers list contains paper **#142**, "Physics-Informed Neural Networks as Solvers for the Time-Dependent
  Schrödinger Equation", by Shah, K.; Stiller, P.; Hoffmann, N.; Cangi, A.
- The linked PDF (`files/NeurIPS_ML4PS_2022_142.pdf`) is reachable (HTTP 200, application/pdf).

**New entry.**

- Shah, K., Stiller, P., Hoffmann, N., & Cangi, A. (2022, 3 Aralık). *Physics-informed neural networks as solvers for
  the time-dependent Schrödinger equation* [Bildiri]. Machine Learning and the Physical Sciences Workshop, NeurIPS 2022,
  New Orleans, LA, ABD. https://ml4physicalsciences.github.io/2022/files/NeurIPS_ML4PS_2022_142.pdf
- The in-text form "Shah vd., 2022" is unchanged.

**Internal note.** The workshop page's metadata says "35th" NeurIPS while its body text says "36th annual conference";
the entry therefore uses "NeurIPS 2022" without an ordinal.

### 6.3 APA recheck

- All 41 entries are present and all are cited (checked by the build); no entry contains "vd." or editorial notes.
- **Verification status:**

  | Source of verification | Count |
  |---|---|
  | Crossref | 31 |
  | proceedings pages | 3 |
  | arXiv | 1 (Kingma & Ba) |
  | official workshop page | 1 (Shah et al., new in R2) |
  | not verifiable online here | 5: four books, Holm 1979 |

---

## 7. Scientific consistency (R2 vs R1 vs R0.3)

| Check | Result |
|---|---|
| Numerical values | the set of decimal numbers in the report body is **identical** to R1 (no value added or removed), and R1 was verified against R0.3 plus `analysis.json`; percentages identical |
| Confidence intervals and statistics | unchanged (Tablo 4, 5, 6, 7, 8, 10; Holm p = 1,5×10⁻⁴; argmin shares; vertex; Welch p = 0,086) |
| Validation statuses | identical counts in R1 and R2: "YAPILMADI"/"Yapılmadı" (B4), "DOĞRULANMADI" (W/Z), "ÇALIŞTIRILMADI" (coupling-scaling), "BEKLEMEDE" (R0); "B4 yapılmadı" 3× |
| Localization | "lokalizasyon kanıtı …" 4×; "Anderson lokalizasyonu gösterilmemiştir" 1×; no independent stationary validation of disordered transport claimed |
| PINN | still secondary ("ikincil" 4×); wording now limited to the closed-box problem (§5) |
| Negative findings and limitations | all kept |
| Outdated claims | ENAQT: 0 occurrences |

---

## 8. Remaining placeholders

**In the report** (8 occurrences of `[ÖĞRENCİ TARAFINDAN TAMAMLANACAK]`, plus the schedule's Ö cells):

| # | Item |
|---|---|
| 1 | Claude version |
| 2 | access dates |
| 3 | other AI tools |
| 4 | confirmation of the purpose list |
| 5 | the student's further contributions |
| 6 | checks actually performed by the student and advisor |
| 7 | repository visibility and licence |
| 8 | months of activities after June, including the month the research question was defined, and plans for B4 and the ζ campaign |

**Required but not marked in the text:**

- Update the commit counts at submission. The declaration's "103 / 79" describes the repository before the R1 commit;
  it is now 105 / 81 including the R1 and R2 commits.
- Confirm the "Kuantum Teknolojileri" thematic field.
- Check the page count in Microsoft Word.

**Still open from earlier audits:**

- the North 2003 record (§6.1);
- the 5 references not verifiable online;
- the English supplementary figures;
- the keyword interpretation;
- the barrier sentence;
- a targeted literature check;
- the `CAPTIONS.md` / R0.3 "kζ ≲ 1.4" wording.

**Eligibility.** No eligibility determination was made; the AI declaration must be reviewed by the student and advisor.
