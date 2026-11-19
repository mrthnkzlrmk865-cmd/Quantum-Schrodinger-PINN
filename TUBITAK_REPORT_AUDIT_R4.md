# TÜBİTAK Report Audit — R4

**Document audited:** `TUBITAK_REPORT_DRAFT_R4.docx`, the final readability and compliance revision of R3.

**Baselines:**

- R3 (`TUBITAK_REPORT_DRAFT_R3.docx`, commit `2c4b617`) is the scientific and textual baseline.
- `RESEARCH_REPORT_DRAFT_R0_3.md` remains the scientific source of truth.

**What was not done or changed:**

- No simulation, training, benchmark, validation test or new numerical experiment was run.
- None of these were modified: R0–R3, R0.3, the mapping, source code, datasets/results, figures, validation records,
  README, existing supplementary documents.

**New files:**

| File | Content |
|---|---|
| `TUBITAK_REPORT_DRAFT_R4.docx` | the R4 report |
| `TUBITAK_REPORT_AUDIT_R4.md` | this audit |
| `TUBITAK_SUPPLEMENTARY_MATERIAL_R4/EK_BELGE_3_PINN_Ayrintilari_R4.docx` and `.md` | EK BELGE 3 with corrected 3D captions (§7) |
| `TUBITAK_SUPPLEMENTARY_MATERIAL_R4/README.md` | which supplementary version to upload |

---

## 1. Final page count

**Result: 20 pages** in LibreOffice 24.2 (headless DOCX → PDF conversion, `pdfinfo`).

**Rendering check.** Every page was rendered to an image and inspected. No page is half-empty: the white space at the
bottom of a page is at most 15% (page 11, before a section heading). The text on page 20 ends a few lines above the
bottom margin.

**Layout changes from R3.** The R4 edits added about 8 lines, so the following changes were made. No font size was
reduced and no scientific content was removed.

| Change | Detail |
|---|---|
| Figure widths | Şekil 4: 13.2 → 12.6 cm; Şekil 5: 13.6 → 13.0 cm; Şekil 6: 10.0 → 9.4 cm |
| Tablo 11 | reflowed into two balanced column pairs, so it no longer has 5 empty half-rows |
| Ekler list | shorter descriptions of each EK BELGE; every listed item is kept |
| Wording | three small tightenings with no change of meaning: limitation (5) of the length scan; Öneriler item 3; the Yöntem definition of "realizasyon", which is now defined at its first use instead (§4) |

**Pagination caveat.** Microsoft Word may paginate differently. **Check the page count in Word** before submission (§9,
blocker 1).

---

## 2. Official template compliance

| Check | Result |
|---|---|
| Template labels | the 19 template labels are unchanged (checked by the build against a snapshot of the official template) |
| Section order | unchanged from R3 |
| Title and fields | title, Fizik and Kuantum Teknolojileri unchanged |
| Font | Times New Roman only |
| Sizes | body, tables and captions 12 pt; the template's own schedule table 9–10 pt (as in R2 and R3); empty spacer paragraphs 6 pt |
| Margins | 2.5 cm on all sides |
| Spacing and alignment | single line spacing; body paragraphs justified (references left-aligned with a hanging indent, as in R3) |
| Colour | all text black |
| Abstract | unchanged from R3: 213 words (limit 150–250), no abbreviations |
| Keywords | unchanged from R3: 5 (limit ≤ 5) |
| Schedule table | 1 + 10 columns, as in the template; the month labels are adapted (see §6, an open issue) |

---

## 3. Abbreviations and symbols table review (Tablo 11)

**Removed:**

| Entry | Reason |
|---|---|
| B1–B4 | internal test identifiers. They are no longer used in the main text, except "B4", which is explained in the text at first use (Tablo 1) and in Tablo 9, and the pointer to EK BELGE 4 in Ekler |
| W/Z | internal identifier; removed from the main text (§5) |
| L-BFGS | implementation term; it occurs once and is explained there ("sınırlı bellekli bir yarı-Newton yöntemi olan L-BFGS") |
| CPU | implementation term; explained once in the text ("yalnızca ana işlemci, CPU; ekran kartı kullanılmamıştır"); elsewhere the text says "ana işlemci" / "işlemci" |
| RGF | not essential to the main narrative; the abbreviation was removed from the main text, which now says "özyinelemeli Green fonksiyonu yöntemi" |

**Kept:**

| Entry | Reason |
|---|---|
| PR (katılım oranı) | used meaningfully: equation (7), the closed-box section, and the column headings of Tablo 8 |
| CAP | used in the text and in the symbol W_CAP |
| DST, PINN | used throughout the main text |
| GA, sd | used in table headings |
| L2 | used in Tablo 10 and in the PINN text |
| ÜYZ | used in the AI declaration |
| 1B/2B/3B | used in the text |

**Corrected explanations:**

- **L2:** now reads exactly "İkinci norm; bağıl hata ölçümünde kullanılır." It is no longer described as "karesel norm".
- **kℓ_tr:** "dalga sayısı ile taşınım ortalama serbest yolunun çarpımı; büyük değer: daha zayıf saçılma". This matches
  how the report uses it (kℓ_tr ≥ 10 is weak scattering, kℓ_tr < 10 stronger scattering).
- **⟨T⟩:** "örneklerin ortalama iletimi".
- **n:** "örnek (realizasyon) sayısı".
- **χ²_red:** "(uydurma kalitesi)" added.
- **DST:** the redundant "(klasik çözücü)" was dropped.

**Final content (23 entries).**

- **9 abbreviations:** 1B/2B/3B, CAP, DST, GA, L2, PINN, PR, sd, ÜYZ.
- **14 symbols:** ħ; ψ, |ψ|²; V, W; ζ; k, k_x0; E₀; L_d, L_y; T, R, P_in; ⟨T⟩; ℓ, ℓ_tr; kℓ_tr; ξ, N; n; χ²_red.
- Every requested essential symbol is present: ħ, ψ, ζ, k, W, L_d, T and R, ℓ and ℓ_tr, ξ, ⟨T⟩, DST, PINN.
- No symbol was added to lengthen the table. Each entry is used in the main report.

---

## 4. Terminology and remaining unexplained terms

**Edits made:**

| Item | R3 | R4 |
|---|---|---|
| realizasyon (first use, Tablo 1) | "her biri 48 realizasyon" (defined only later, in Yöntem) | "her biri 48 realizasyon (rastgele üretilmiş düzensizlik örneği)"; the Yöntem definition now reads "… her alana bir düzensizlik realizasyonu (örnek) denir" |
| later uses | "realizasyon" 28× | "örnek" used where it fits: "İletimin örnekten örneğe değişimi", the bootstrap explanation, Tablo 11. "realizasyon" now appears 22×, mostly in fixed expressions and in table titles |
| bootstrap explanation | "eldeki realizasyonlardan yerine koyarak çok sayıda yeni örnek çekilip …" | "eldeki örneklerden yerine koyarak çekilen çok sayıda yeni veri kümesi oluşturulup …"; this avoids using "örnek" for two different things |
| spectral bias (Giriş) | "ağların yavaş değişen bileşenleri daha kolay öğrenmesi (spektral yanlılık)" | "sinir ağlarının bir fonksiyonun yavaş değişen kısımlarını hızlı salınan kısımlarından daha önce ve daha kolay öğrenme eğilimi (spektral yanlılık)" |
| "fit" (6×) | English | "uydurma" / "eğri uydurma" (Tablo 1, Yöntem, Tablo 7, Şekil 5, Bulgular) |
| "quasi-Newton" | English | "yarı-Newton" |
| "CPU" | 5× | explained once; elsewhere "ana işlemci" / "işlemci" |
| "coupled-block", "coupling-scaling" | English, in Tablo 9 | removed from the main text (§5) |
| "V7 dalga sayısı M" (Tablo 3) | ambiguous: in the report "dalga sayısı" otherwise means k | "düzensizlikteki kosinüs sayısı M" (M is the number of cosine terms in equation 3) |

**Kept unchanged** (standard scientific terms, each explained at first use): Born yaklaşımı, Anderson lokalizasyonu,
korelasyon uzunluğu, saçılma, ortalama serbest yol, bootstrap, Holm düzeltmesi, post-hoc (explained once as "sonuçlar
incelendikten sonra seçilmiştir"), Strang adım-bölme, Fisher–Lee bağıntısı, katılım oranı.

**Remaining expressions a reader may still find hard.** These are not errors; the student should be ready to explain
them.

1. **Şekil 6 row labels** (clean, s, d1, nofourier, rad, T0.01/T0.04, long). They are in English because they are
   printed inside the existing figure file, which was not to be modified. The caption explains each one in Turkish.
2. **"commit"** in the AI declaration and the schedule text. These sections were left unchanged (§8). A Turkish gloss
   such as "sürüm kaydı (commit)" could be added by the student.
3. **"post-hoc", "bootstrap", "Holm düzeltmesi", "(b+1)/(B+1)".** These are standard statistical terms, explained in the
   text. The student should be able to say why p = 1,5×10⁻⁴ is at the resolution limit of B = 20 000.
4. **kζ, kℓ_tr, τ_s, τ_tr, I₀/I₁ (modified Bessel functions).** Defined in Yöntem. They are mathematical background
   above high-school level.
5. **"yarı-1B lokalizasyon uzunluğu", var(ln T), ⟨ln T⟩.** Explained briefly in Bulgular.
6. **PINN training terms:** kalıntı, Fourier öznitelikleri, Adam, L-BFGS, kalıntı-uyarlamalı örnekleme, nedensel
   ağırlıklandırma. Each has a short explanation in Yöntem.

---

## 5. Internal validation identifiers

**Removed from the main text:** B1, B1-F, B2, B3, B3X, B3Y, W/Z, coupled-block, coupling-scaling, R0, RGF, and the
V1–V7 labels in Tablo 3.

**Remaining:**

- **"B4"**, the name of the test that was not performed. It is explained where it first appears and in every place it
  is used.
- **One pointer in Ekler:** "B1–B4 test kimlikleri" for EK BELGE 4, and "V1–V7 olarak numaralanan testler" for
  EK BELGE 1. This lets a reader match the main report with the supplements.

**Tablo 9** ("Bağımsız durağan çözücüyle yapılan testlerin durumu") now describes each test in plain words. Every status
is unchanged:

| R4 row (plain description) | R3 identifier | Status (unchanged) |
|---|---|---|
| Düzensizlik içermeyen kılavuz (6 yapılandırma) | B1 | passed in 5 configurations; phase criterion (10⁻¹⁰) **BAŞARISIZ** in the longest configuration (B1-F) |
| Düzensizlik içermeyen bariyer | B2 | **GEÇTİ** |
| Sabit düzensizlik örnekleri: two independent references | B3, B3X | **BAŞARISIZ** in both (ızgara inceltme ölçütü) |
| İnce ızgaralı tekrar | B3Y | **GEÇTİ** |
| Düzensiz taşınımın çapraz doğrulaması (B4) | B4 | **YAPILMADI** |
| Sayısal yuvarlama teşhis çalışmaları | W/Z; coupled-block; coupling-scaling; R0 | eşdeğer yazım **DOĞRULANMADI**; one experiment's design found insufficient; the other **ÇALIŞTIRILMADI**; independent review **BEKLEMEDE** |

**Correction made while rewriting the table.**

- "A–E geçti; B1-F" in R3 refers to **six configurations, A–F**, not to five criteria. This is confirmed by
  `docs/RGF_B1_VALIDATION.md`: "all B1 criteria pass … in 5 of 6 tested configurations".
- The R4 wording therefore says "5 yapılandırmada geçti; en uzun yapılandırmada faz ölçütü … BAŞARISIZ".

**Status counts.**

- Every status word except BAŞARISIZ appears the same number of times as in R3: YAPILMADI 1, DOĞRULANMADI 1,
  ÇALIŞTIRILMADI 1, BEKLEMEDE 1, GEÇTİ 2.
- **"BAŞARISIZ" appears 2× (R3: 3×)** because the two failed disorder tests (B3, B3X) now share one row ("her iki
  karşılaştırmada BAŞARISIZ").
- All three failures are still stated. Sonuç, limitation (8), names them in words: "bağımsız durağan çözücü
  testlerinden üçü … başarısızdır".

**What the reader can still see clearly:**

- the solver's own validation and convergence checks passed (Tablo 3);
- the independent barrier test passed;
- three independent stationary tests failed;
- B4 was not performed;
- disordered transport has **not** been independently cross-validated. This is stated in bold in Yöntem and in the
  Tablo 9 section, and repeated in Sonuç and Öneriler.

All identifiers and measurements remain in EK BELGE 4 (unchanged).

**Schedule row.** "RGF B1–B3Y (kısmen; B4 yapılmadı)" is renamed "Bağımsız durağan çözücü testleri (kısmen; B4
yapılmadı)". Its cells are unchanged.

---

## 6. Work-time schedule compliance

**Comparison with the official 2026–2027 template** (`proje_yazma_sablonu_2026-2027.docx`):

| Item | Official template | R3 / R4 |
|---|---|---|
| Structure | "İşin Adı" + "AYLAR" header over 10 month columns | same structure (1 + 10 columns) |
| Month columns | NİSAN, MAYIS, HAZİRAN, TEMMUZ, AĞUSTOS, EYLÜL, EKİM, KASIM, ARALIK, OCAK | HAZİRAN, TEMMUZ, AĞUSTOS, EYLÜL, EKİM, KASIM, ARALIK, OCAK, ŞUBAT, MART (since R1) |
| Example rows | Literatür Taraması, Arazi Çalışması, Verilerin Toplanması ve Analizi, Proje Raporu Yazımı | the project's own 10 rows |
| Explanatory note | none | a note above the table: start in June 2026, the meaning of the symbols, which months still need to be filled in |

**What the template says.** "aşağıda sunulan **tabloya benzer** bir iş-zaman çizelgenizin olması gerekir" (you must
have a schedule *similar to* the table below). The template also calls for a schedule "ilk haftadan itibaren …
projenizin süresi boyunca" (from the first week, for the duration of the project).

**Assessment.**

- This wording suggests that the table is a model and that its months may follow the actual project period.
- However, nothing in the template or the provided evaluation criteria says so explicitly. The month labels may also
  have been chosen to match the competition calendar.
- **This could not be determined from the available documents.**
- Because the question is uncertain, R4 does **not** change the R3 month columns again. It keeps the R1–R3 adaptation
  and records it here as an open issue.

**Open issue S-1 — decision for the student and advisor.** Choose one option:

- **(a) Keep the adapted columns, HAZİRAN–MART.** These match the documented start in June 2026.
- **(b) Restore the template's NİSAN–OCAK columns.** NİSAN and MAYIS would then be marked "—" (the project had not yet
  started). The current ŞUBAT and MART cells, which are only "Ö" placeholders, would be dropped. That is a small change in
  `build_report.py` (the month labels and cell lists) and alters no fact.

If the call text or the regional coordinator gives a rule, follow it.

**Dates.** No activity date was invented.

- Only June 2026 (literature review; first 3D Plotly code) carries ✓, based on the student's statement.
- November 2026 report writing is marked "D" (devam ediyor).
- All other activity months remain "Ö" or "—", as in R3.

---

## 7. 3D visualizations

- **Main report.** The 3D PINN–DST visualizations are **not** in the main report. The main report has only the six 2D
  figures (Şekil 1–6). The 3D images are mentioned only in the Ekler list ("3B PINN–DST karşılaştırma görselleri", in
  EK BELGE 3) and in one sentence of the PINN section that points to EK BELGE 3.
- **EK BELGE 3 (R0 version).** Its section text already said "3B yüzeyler 2B alanların gösterimidir, üç boyutlu
  fiziksel simülasyon değildir". But the captions of Şekil EK3.1 and EK3.2 said only "3B yüzey yalnızca bir gösterim
  biçimidir", and the EK3.2 caption referred to "Şekil 2" instead of Şekil EK3.1.
- **New `EK_BELGE_3_PINN_Ayrintilari_R4`.** Each of the three 3D captions now states: "Bu görsel, iki boyutlu bir sayısal
  alanın 3B yüzey gösterimidir; … Üç boyutlu bir fiziksel sistemin simülasyonu değildir." The cross-reference is
  corrected to Şekil EK3.1.
- **Nothing else in EK BELGE 3 was changed:** the comparison with R0 gives 3 changed lines, the same 5 images, and every
  number identical. The R0 file is unchanged.
- **Upload EK BELGE 3 from `TUBITAK_SUPPLEMENTARY_MATERIAL_R4/`.**

---

## 8. AI declaration placeholders

The AI declaration is **identical to R3** (checked paragraph by paragraph). No student-specific information was
invented or filled in, and no independent human verification is claimed.

**Remaining student inputs** (8 `[ÖĞRENCİ TARAFINDAN TAMAMLANACAK]` placeholders in the report, the same as R3):

| # | Placeholder | Where |
|---|---|---|
| 1 | Claude version (as shown in the tool) | Kullanılan araçlar |
| 2 | access dates / period | Kullanılan araçlar |
| 3 | other AI tools, if any (name, version, date) | Kullanılan araçlar |
| 4 | confirmation or correction of the list of AI uses | Kullanım amacı |
| 5 | the student's further contributions | Kullanım kapsamı (A) |
| 6 | checks actually performed by the student and by the advisor (human verification, advisor review) | İnsan müdahalesi (C) |
| 7 | repository visibility and licence | Yöntem, Tekrar üretilebilirlik |
| 8 | months of activities after June, the month the research question was defined, and whether B4 and the ζ campaign are planned | schedule note and "Ö" cells |

**Facts in the declaration that the student must update or confirm (not marked as placeholders):**

- **Commit counts.** The declaration says "103 birleştirme dışı commit'in 79'unun yazarı Claude'dur". The text was left
  unchanged, because the declaration's facts were not to be edited. Current counts:
  - before the R4 commit: 106 non-merge commits, 82 by Claude;
  - after the R4 commit: 107 / 83.
  - **The student must update the numbers at submission time.**
- The statement that the student defined the initial research question, and the description of the student's own work.
  These are based on what the student told us; only the student can confirm them.

---

## 9. Scientific consistency with R3 and R0.3

| Check | Result |
|---|---|
| Decimal numbers | the set of decimal numbers in the document XML is **identical** to R3: 215 distinct values, none added or removed. R3 was identical to R2, and R2 was verified against R0.3 and `analysis.json` |
| Equations | 8 equations, text **identical** to R3 |
| References | list **identical** to R3: 41 entries, all cited ("cited refs: 41 uncited: []") |
| Abstract, keywords, AI declaration, field labels, schedule text | identical to R3 |
| Schedule cells | identical to R3; only one row label renamed (§5) |
| Validation statuses | unchanged (§5) |
| Localization wording | "lokalizasyon kanıtı değildir" 3×, "lokalizasyonun kanıtı değildir" 1×, "Anderson lokalizasyonu" 6×, including "Bu çalışmada Anderson lokalizasyonu gösterilmemiştir" (all as in R3) |
| Qualifications | "keşifsel" 8×, "anlamlı değildir" 4×, "post-hoc" 1×, "mertebe" 15× (all as in R3) |
| Outdated claims | ENAQT: 0 occurrences |
| Limitations | all 11 items kept; item (8) now names the three failed tests in words |

**Claims that do not appear** (checked by reading every changed passage):

- 3D physical simulations;
- demonstrated Anderson localization;
- experimental verification;
- completion of B4;
- validation of W/Z.

## 10. Confirmation that no numerical results or scientific claims changed

None of the following was changed:

- numerical results, confidence intervals and statistical tests;
- equations and the physical model;
- conclusions and limitations;
- validation statuses;
- references.

**All edits are of these kinds:**

- terminology and explanations (§4);
- replacing internal identifiers with plain descriptions, with the same statuses (§5);
- the abbreviations table (§3);
- figure widths and table column widths (§1);
- EK BELGE 3 captions (§7).

**One factual precision was added.** The B1 result is now described accurately as 5 of 6 *configurations* passing (§5).
This is the meaning of R3's "A–E geçti" and agrees with `docs/RGF_B1_VALIDATION.md` and R0.3.

---

## 11. Remaining submission blockers

| # | Blocker | Owner |
|---|---|---|
| 1 | **Page count in Microsoft Word.** R4 is 20 pages in LibreOffice with little room to spare. If Word shows 21 pages, reduce one figure by about 0.5 cm or shorten the Ekler list. Do not reduce the body font | student |
| 2 | **AI declaration placeholders 1–6** (§8) | student and advisor |
| 3 | **Commit counts** in the AI declaration (103 / 79 → current values at submission) | student |
| 4 | **Repository visibility and licence** (placeholder 7) | student |
| 5 | **Schedule:** the months after June and the month the research question was defined (placeholder 8), and **open issue S-1** (adapted versus template month columns, §6) | student and advisor |
| 6 | **Thematic field.** Confirm that "Kuantum Teknolojileri" is valid in the 2026–2027 call list | student and advisor |
| 7 | **Advisor review** of the full report and of the AI declaration. **No eligibility determination was made** | advisor |
| 8 | **Supplementary uploads.** EK BELGE 1, 2, 4, 5, 6 from R0; EK BELGE 3 from R4; EK BELGE 7 from R1 (or drop it, since Tablo 11 now covers it, and remove "EK BELGE 7. Kısaltmalar." from Ekler). The EK BELGE 1–6 figures still have English labels | student |
| 9 | **References still open:** the North (2003) same-title record (erratum or not); 5 entries not verifiable online here (four books, Holm 1979) | student or advisor (library check) |
| 10 | **Oral preparation** for the terms in §4 ("remaining expressions") and the validation statuses in §5 | student |
