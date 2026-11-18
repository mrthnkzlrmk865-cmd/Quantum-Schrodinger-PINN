# TÜBİTAK Report Audit — R3

**Document audited:** `TUBITAK_REPORT_DRAFT_R3.docx` (readability revision of R2).

**Purpose of R3.** Make the report readable for a Turkish high-school jury without weakening its scientific rigor:

- explain specialised terms on first use;
- use shorter sentences;
- rewrite the abstract;
- reorder the Introduction;
- lead each result with its main finding;
- add a table of abbreviations and symbols.

**Baselines:**

- R2 (`TUBITAK_REPORT_DRAFT_R2.docx`, commit `09e0d0b`) is the numeric and textual baseline.
- `RESEARCH_REPORT_DRAFT_R0_3.md` remains the scientific source of truth.

**What was not done or changed:**

- No simulation, training, benchmark, validation script or new numerical experiment was run.
- None of these were modified: R0, R1, R2, R0.3, the mapping, source code, datasets/results, original figures, validation
  records, README, supplementary files.
- Only `TUBITAK_REPORT_DRAFT_R3.docx` and this audit are new.

---

## 1. Page count

**Result: 20 pages** in LibreOffice 24.2 (headless DOCX → PDF conversion; `pdfinfo` reports `Pages: 20`).

**Why layout work was needed.**

- The readability rewrite added about 420 words (6 668 → 7 091 words in the document's paragraphs), plus the new
  abbreviations table.
- The first R3 build was 23 pages.

**Layout changes used to return to 20 pages.** No font was reduced, and no scientific content was removed.

| Change | Detail |
|---|---|
| Body paragraphs | first-line indent 0.6 cm and 0 pt space-after (R2: no indent, 2–4 pt space-after); paragraphs stay visually separated by the indent |
| Empty spacer paragraphs | template spacers and spacers after tables set to a 6 pt line height (they were 12 pt) |
| Schedule table | changed from a floating table (`w:tblpPr`) to an inline table so it can split across pages; its first column is set to 5.6 cm. Its cells, rows and month labels are unchanged |
| Figure widths, reduced slightly | Şekil 1: 8.6 cm; Şekil 2: 15.2 cm; Şekil 3: 9.5 cm; Şekil 4: 13.2 cm; Şekil 5: 13.6 cm; Şekil 6: 10.0 cm |
| Figure positions | Şekil 3–6 each moved next to the paragraph that introduces them, to remove half-empty pages; numbering and captions unchanged |
| Table column widths | Tables 3, 9, 10 and 11 rebalanced to remove wrapped lines |
| Ekler | the appendix list is one compact paragraph; its opening clause "Ana rapor sayfa sınırı nedeniyle" was removed |

**Pagination caveat.**

- The count was measured in LibreOffice. Microsoft Word can paginate slightly differently because of font metrics and
  equation rendering.
- Page 20 is full, so a Word rendering could spill a few lines onto page 21.
- **The student must open the file in Word and confirm the page count** before submission (see §13).

---

## 2. Template compliance

| Check | Result |
|---|---|
| Template labels | the 19 template labels are unchanged (the build script compares each label against a snapshot of the official template and stops on any difference) |
| Section order | unchanged: Proje Ana Alanı, Proje Tematik Alanı, Proje Adı, Yapay Zekâ beyanı, Özet, Anahtar Kelimeler, Amaç, Giriş, Yöntem, Proje İş-Zaman Çizelgesi, Bulgular, Sonuç ve Tartışma, Öneriler, Kaynaklar, Ekler |
| Title | "Korelasyon Uzunluğunun İki Boyutlu Kuantum Dalga Paketi Taşınımına Etkisinin Sayısal İncelenmesi" (unchanged) |
| Fields | Fizik / Kuantum Teknolojileri (unchanged) |
| Font | Times New Roman only (all runs) |
| Body size | 12 pt; tables and captions are also 12 pt |
| Other sizes | the template's own schedule table (9–10 pt, as in R2) and the 6 pt empty spacer paragraphs |
| Margins | 2.5 cm on all sides |
| Line spacing | single (1.0) |
| Alignment | body paragraphs justified; references are left-aligned with a hanging indent (APA), as in R2; one table note is left-aligned |
| Colour | all text black (the template's red instruction colour removed, as in R2) |
| Template instruction paragraphs | 19 removed (as in R2) |
| New element | "Kısaltmalar ve Semboller" heading and Tablo 11, placed after Kaynaklar and immediately before Ekler. It is not a template label, so the template's label set is unchanged |

---

## 3. Abstract word count

**R3: 213 words** (R2: 194). Required range: 150–250.

**What was rewritten:**

- No abbreviations are used. The 2-or-more-capital-letter tokens DST, CAP, PINN, RGF etc. no longer appear.
- Correlation length is explained as "tepe ve çukurlarının tipik genişliği".
- Transmission is explained as "düzensiz bir bölgeyi geçme olasılığı".
- The Born prediction is stated in words ("dalganın ilk yönünü kaybetmeden önce aldığı yol").
- Realizations are described as "rastgele potansiyelin bağımsız örneği".
- kζ is written as "dalga sayısı ile korelasyon uzunluğunun çarpımı".

**Qualifications kept in the abstract:**

- exploratory analysis;
- post-hoc selection;
- no pre-planned validation experiment;
- 7 sampled correlation lengths;
- coarse minimum;
- no independent comparison for disordered transport;
- no Anderson localization;
- PINN secondary and limited to the closed box and its architecture and budget.

**Numbers in the abstract:** 0,494 / 0,670 / 0,734 and 24–48, identical to R2.

---

## 4. Keyword count

**5 keywords, unchanged from R2:** Kuantum taşınımı; Korelasyonlu düzensizlik; Korelasyon uzunluğu; Zamana bağlı Schrödinger
denklemi; Born yaklaşımı.

---

## 5. Abbreviations and symbols table review

**Format.** Tablo 11, "Raporda kullanılan kısaltmalar ve semboller": 4 columns (Kısaltma | Açıklama | Sembol | Açıklama),
14 rows.

**Abbreviations (14 entries, alphabetical):** 1B/2B/3B, B1–B4, CAP, CPU, DST, GA, L-BFGS, L2, PINN, PR, RGF, sd, ÜYZ, W/Z.

**Symbols (14 entries):** ħ; ψ, |ψ|²; V, W; ζ; k, k_x0; E₀; L_d, L_y; T, R, P_in; ⟨T⟩; ℓ, ℓ_tr; kℓ_tr; ξ, N; n;
χ²_red.

**Checks:**

- Every abbreviation in the table occurs in the report text or tables outside Tablo 11.
  - Occurrence counts: 1B 3, 2B 7, 3B 4, B1 7, B4 14, CAP 1, CPU 4, DST 4, GA 6, L-BFGS 1, L2 6, PINN 39, PR 7, RGF 4,
    sd 3, ÜYZ 4, W/Z 5.
- No abbreviation was invented. Each one was already used in R2 or is the standard short form of a term R2 already used.
- "CAP" is now expanded at its first occurrence in Yöntem: "karmaşık soğurucu potansiyel (CAP)". In R2 it appeared only as
  the subscript of W_CAP.
- "B1–B4" says "(B4 yapılmadı)"; W/Z says "RGF'nin eşdeğer yazımı (W'den farklı)". This avoids confusion with the
  disorder strength W.
- **Not included:** abbreviations used only once inside a figure file name or a source path (e.g. `results/...`), and the
  template's own labels. Statistical symbols that are explained where they are used (B for the number of bootstrap
  resamples, p) are also not included.

**Redundancy note.**

- EK BELGE 7 ("Kısaltmalar", an existing supplementary file) now overlaps with Tablo 11.
- It was not changed, because the existing supplementary files were not to be modified, and it is still listed in Ekler.
- The student may drop EK BELGE 7 from the upload, or keep it as a longer glossary. If dropped, remove "EK BELGE 7.
  Kısaltmalar." from the Ekler paragraph.

---

## 6. Readability edits made

| Section | Edit |
|---|---|
| Özet | fully rewritten (§3) |
| Amaç | shorter sentences; the research question and the gap stated separately; the duplicated gap paragraph merged |
| Giriş | reordered (see below); each paragraph opens with its main idea; one sentence on "lazer ışığının rastgele benek desenleri" removed as a side remark (its Kuhn and Sanchez-Palencia citations remain in the same paragraph) |
| Yöntem | still technical, but each term is explained at first use (§7); all equations (1)–(8) and every parameter are unchanged |
| Bulgular | every result block starts with a bold "**Ana bulgu.**" sentence, followed by the numbers, their meaning, and what they do **not** establish ("Bu sonuç ne gösterir, ne göstermez?") |
| Sonuç ve Tartışma | condensed into bold-led paragraphs: research-question answer; goals and theory comparison; contribution and comparison with existing approaches; limitations (11 numbered items, all kept); applicability |
| Kısaltmalar ve Semboller | new (§5) |
| Ekler | one compact paragraph |

**New Giriş order:**

1. phenomenon (wave scattering in a random medium);
2. effects of disorder (diffusive spreading versus Anderson localization; mean free path; scaling theory);
3. correlation length;
4. the Born prediction;
5. geometry and limitations (stationary methods versus the time-dependent approach, finite guide);
6. PINN background;
7. research question.

**Unchanged sections:** the AI declaration, field labels, keywords, schedule text, schedule table and reference list are
**identical** to R2 (checked paragraph by paragraph and cell by cell).

---

## 7. Examples of terminology explanations

Standard scientific terms are **kept** and explained; none was replaced by an informal word.

| Term | R2 | R3 (first use) |
|---|---|---|
| realizasyon | "düzensizlik realizasyonu" without definition | "Aynı W ve ζ ile farklı tohumlarla üretilen her alana bir düzensizlik realizasyonu denir." |
| ansambl | "ansambl ortalaması" (3×) | replaced by the precise equivalent "realizasyonlar üzerinden ortalama"; ⟨T⟩ is defined in Tablo 11 |
| bootstrap | "bootstrap güven aralıkları" | "bootstrap yöntemiyle, yani eldeki realizasyonlardan yerine koyarak çok sayıda yeni örnek çekilip ortalamanın bu örnekler arasındaki dağılımına bakılarak" |
| post-hoc | "eşleştirilmemiş sonuçlar görüldükten sonra (post-hoc) seçilmiştir" | "sonuçlar incelendikten sonra seçilmiştir (post-hoc); … önceden planlanmış bir doğrulama testi değildir" |
| Holm correction | named only | "birden fazla karşılaştırmada yanlışlıkla 'anlamlı' sonuç bulma olasılığını denetlemek için … Holm düzeltmesi" |
| CAP | "soğurucu sınır" | "dalgayı geri yansıtmadan yutması amaçlanan bir karmaşık soğurucu potansiyel (CAP)" |
| DST | "ayrık sinüs dönüşümüyle (DST)" | the same, plus "(raporda 'DST çözücüsü' veya 'klasik çözücü')" |
| Strang splitting | "Strang adım-bölme" | "her zaman adımı, kinetik ve potansiyel enerjinin etkilerini ayrı ayrı içeren alt adımlara bölünür" |
| Fisher–Lee | "Fisher–Lee bağıntısıyla hesaplanır" | "iletimi Green fonksiyonundan veren Fisher–Lee bağıntısı" |
| RGF | "RGF" | "özyinelemeli Green fonksiyonu (RGF) yöntemi" |
| B4 | "B4" | "zamana bağlı ve durağan çözümün karşılaştırılması (B4)" in Tablo 1; "Yapılmadı" kept |
| W/Z | "W/Z teşhisi" | "yöntemin cebirsel olarak eşdeğer bir yazımı (W/Z)"; Tablo 9: "yöntemin sayısal yuvarlama davranışının incelenmesi" |
| PINN | "fizik-bilgili sinir ağları" | "bir diferansiyel denklemin çözümünü yapay bir sinir ağıyla temsil eder ve ağı denklemin seçilen noktalarda sağlanacağı biçimde eğitir" |
| spectral bias | "spectral bias" | "ağların yavaş değişen bileşenleri daha kolay öğrenmesi (spektral yanlılık)" |
| L-BFGS | "200 L-BFGS iterasyonu" | "sınırlı bellekli bir quasi-Newton yöntemi olan L-BFGS'nin 200 iterasyonu" |
| kℓ_tr | symbol only | "kℓ_tr ≫ 1 zayıf, kℓ_tr ∼ 1 güçlü saçılmayı gösterir"; Tablo 11 "saçılma ölçüsü (büyük: zayıf)" |
| χ²_red | symbol only | "değerin 1'e yakın olması iyi uyumu, 1'den belirgin biçimde büyük olması verilerle tam uyuşmamayı gösterir" |
| PR | "katılım oranı" | "olasılık yoğunluğunun etkin olarak kapladığı alan", with how it behaves for localized and spreading packets |
| relative L2 error | "bağıl L2 hatası" | "PINN ile klasik çözüm arasındaki farkın büyüklüğünün (L2 normu) klasik çözümün büyüklüğüne oranı"; a value of 1 equals the error of predicting zero everywhere |

**One duplicated phrase removed.** In Yöntem, "(mertebe kıyası)" was removed after "büyüklük mertebesi düzeyinde",
because it repeated the same phrase. The Bulgular sentence "Bu karşılaştırma bir **mertebe kıyasıdır**" is kept.

---

## 8. Terms still needing the student's own explanation

A jury may ask about these. The report defines them, but the student should be able to explain each one aloud:

1. **Realizasyon / tohum (seed).** "Realizasyon" first appears in Tablo 1 (Amaç) and is defined only later, in Yöntem.
   The student should explain it on first mention in a presentation.
2. **Eşleştirilmiş tasarım (paired design).** Why the same seed across ζ values reduces noise but makes the points
   dependent.
3. **Bootstrap, Holm correction and (b+1)/(B+1) p-values.** Why p = 1,5×10⁻⁴ is at the resolution limit of
   B = 20 000 resamples.
4. **kℓ_tr, kζ, ℓ, ℓ_tr, τ_s, τ_tr.** The difference between the scattering rate and the transport rate.
5. **Born approximation and the Born + diffusion estimate.** Why it is only a qualitative, order-of-magnitude comparison
   (unbounded 2D medium versus a 6-channel guide).
6. **CAP and spurious reflection** at k_x0 = 12.
7. **Strang splitting and DST.** Why the kinetic step is exact in the sine basis.
8. **Anderson localization, ξ ∼ Nℓ (quasi-1D), var(ln T).** Why L_d/ξ ≤ 2,6 and χ²_red > 1 prevent a localization
   claim.
9. **Participation ratio (PR)** in the closed-box test.
10. **RGF, Fisher–Lee, the B1–B4 stages, B3X/B3Y, W/Z, coupled-block, coupling-scaling, R0.** What each status means,
    especially "B4 yapılmadı" and "W/Z doğrulanmadı".
11. **PINN terms:** residual, Fourier features, Adam / L-BFGS, residual-adaptive sampling, causal weighting, spectral bias.
12. **Pre-registered (önceden kaydedilmiş) campaign:** designed but not frozen and not run.

---

## 9. Consistency with R2 and R0.3

| Check | Result |
|---|---|
| Decimal numbers | the set of decimal numbers in the document XML is **identical** to R2: 215 distinct values, none added or removed. R2 was verified against R0.3 and `analysis.json` |
| Equations | 8 OMML equations, text **identical** to R2 |
| Tables | R2's 11 tables kept with identical values; 1 new table (Tablo 11) |
| Figures | 6 figures, same files (`visualizations/tubitak_r1/*_tr.png`), same numbering and captions; only widths and positions changed |
| References | reference list **identical** to R2: 41 entries, all cited (build check "cited refs: 41 uncited: []") |
| AI declaration, schedule, keywords, field labels | identical to R2 |
| Validation statuses | identical counts in R2 and R3: "YAPILMADI" 1, "DOĞRULANMADI" 1, "ÇALIŞTIRILMADI" 1, "BEKLEMEDE" 1, "BAŞARISIZ" 3, "GEÇTİ" 2; "B4 yapılmadı" 3× (as in R2: Tablo 3, Tablo 11, the schedule row); "Yapılmadı" for H6/B4 in Tablo 1 and "YAPILMADI" for B4 in Tablo 9 |
| Localization wording | "lokalizasyon kanıtı değildir" 3× and "lokalizasyonun kanıtı değildir" 1× (R2: 4× in total); "Anderson lokalizasyonunun kanıtı değildir" in the abstract; "Bu çalışmada Anderson lokalizasyonu gösterilmemiştir" kept |
| Outdated claims | ENAQT: 0 occurrences |
| Wording counts that fell | "keşifsel" 8× (R2: 10×) and "post-hoc" 1× (R2: 6×). R3 uses the Turkish explanation "sonuçlar incelendikten sonra seçilmiştir" / "sonradan seçilmiştir" (3×) and defines "post-hoc" once. Every place where R2 qualified a result as exploratory or post-hoc still carries the qualification |

---

## 10. Confirmation that numbers and claims are unchanged

None of the following was changed:

- numerical values, confidence intervals and statistical tests;
- equations, the model and its parameters;
- conclusions and validation statuses;
- schedule dates and placeholders;
- the AI-declaration facts;
- references.

**All qualifications are present in R3:**

| Qualification | Where |
|---|---|
| Exploratory analysis | Özet; İstatistiksel analiz; Bulgular intro; Sonuç |
| Post-hoc comparisons | Tablo 5 text; Sonuç limitations (6) |
| Only 7 ζ points | Özet; Bulgular; Sonuç |
| Coarse minimum, exact position not determined | Bulgular "Bu sonuç ne gösterir, ne göstermez?"; Şekil 2 caption; Sonuç |
| ζ = 0,07 contrast not significant after correction | Tablo 5 note; Bulgular; Sonuç |
| No localization claim | Özet; Bulgular (length scan, closed box); Sonuç |
| No independent stationary validation (B4 not performed) | Özet; Yöntem; Tablo 1, 3 and 9; Sonuç |
| PINN secondary and limited to the closed box, architecture and budget | Amaç; Bulgular (the R2 sentence kept verbatim); Sonuç |
| Born comparison qualitative / order of magnitude only | Yöntem; Bulgular; Şekil 3 caption; Sonuç |

The R2 Born sentence ("İncelenen W ve k taramalarındaki örneklenmiş noktalarda oran, kℓtr azaldıkça artan bir eğilim
göstermiştir.") and the R2 PINN sentence are both kept verbatim.

---

## 11. Remaining placeholders

**In the report:** 8 occurrences of `[ÖĞRENCİ TARAFINDAN TAMAMLANACAK]`, the same as R2, plus the schedule's "Ö" cells.

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

- **Commit counts in the AI declaration.** The declaration still says "103 / 79". This fact was deliberately left
  unchanged. Counts at the time of writing:
  - before the R3 commit: 105 non-merge commits, 81 by Claude;
  - after the R3 commit: 106 / 82.
  - The student should update the numbers at submission time.
- Confirm the "Kuantum Teknolojileri" thematic field.
- Confirm the page count in Microsoft Word (§1, §13).

---

## 12. References requiring attention

Unchanged from R2; no reference was edited in R3.

- **North et al. (2002).** A 2003 record with the same title (DOI 10.1086/346173, *AJHG* 72(2), 498–499) cites the
  2002 note.
  - Whether it is an erratum could not be determined.
  - The citation stays as 2002 (71(2), 439–441; DOI 10.1086/341527).
  - The student or advisor should check the journal site.
- **5 entries not verifiable online here:** four books and Holm (1979). They are formatted in APA 7 but must be checked
  against a library record.
- **Shah et al. (2022):** verified on the official ML4PS 2022 page in R2; nothing is open.
- **All other entries:** 31 via Crossref, 3 via proceedings pages, 1 via arXiv. All verified in R1/R2.

---

## 13. 20-page limit check

**Checked: 20 pages in LibreOffice.**

**How it was checked:**

- `soffice --headless --convert-to pdf`, then `pdfinfo`.
- Every page was rendered to an image and inspected on contact sheets.
- The blank space at the bottom of each page was measured. Pages 1–20 are filled to the bottom margin, with at most 16%
  white at the bottom (page 5, before Şekil 1). No page is half-empty.

**Margin of safety: none.** Page 20 ends close to the bottom margin.

**Before submission:**

1. Open `TUBITAK_REPORT_DRAFT_R3.docx` in Microsoft Word and check the page count.
2. If Word shows 21 pages, the simplest fixes that keep the formatting rules are:
   - remove "EK BELGE 7. Kısaltmalar." from Ekler if EK BELGE 7 is dropped (§5);
   - reduce one or two figure widths by about 0.5 cm.
3. Do not reduce the body font.

---

## 14. Sections the student must read carefully

1. **Yapay Zekâ Araçları Kullanım Beyanı.** Unchanged from R2. All placeholders, the commit counts, and the description
   of the student's own contribution must be the student's own, verified statements.
2. **Özet.** Fully rewritten. The student should be able to say every sentence in their own words.
3. **Giriş.** Reordered. Check that the line of argument matches how the student understands the project.
4. **Bulgular, "Ana bulgu" sentences.** These are the shortest statements of the results. Each must be read together with
   the "ne göstermez" sentences that follow it.
5. **Sonuç ve Tartışma, Sınırlamalar (11 items).** The jury will test these.
6. **Tablo 9 and the RGF paragraph.** Validation statuses (B1-F başarısız, B3/B3X başarısız, B3Y geçti, B4 yapılmadı,
   W/Z doğrulanmadı, coupling-scaling çalıştırılmadı, R0 beklemede).
7. **Tablo 11.** Check that every explanation matches how the student will use the term.
8. **İş-Zaman Çizelgesi.** Fill in the "Ö" cells only with dates the student can support.
9. **Ekler.** Make sure each EK BELGE listed is actually uploaded.

**Eligibility.** No eligibility determination was made. The AI declaration must be reviewed by the student and the
advisor.
