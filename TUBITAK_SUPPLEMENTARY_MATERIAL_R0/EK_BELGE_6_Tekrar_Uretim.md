# EK BELGE 6. Tekrar Üretim Komutları, Yazılım Ortamı ve Şekil Açıklaması Notu

*Bu ek belge, “EK BELGE 6. Tekrar Üretim Komutları, Yazılım Ortamı ve Şekil Açıklaması Notu” başlıklı TÜBİTAK 2204-A proje raporunun (Korelasyon Uzunluğunun İki Boyutlu Kuantum Dalga Paketi Taşınımına Etkisinin Sayısal İncelenmesi) ek belgesidir. İçerik, projenin araştırma taslağından (RESEARCH_REPORT_DRAFT_R0_3.md) değiştirilmeden aktarılmıştır; kaynakların tam künyesi ana raporun Kaynaklar bölümündedir. Sayılarda ondalık ayırıcı olarak nokta kullanılmıştır. Hiçbir yeni hesaplama yapılmamıştır.*


**Depo.** `https://github.com/mrthnkzlrmk865-cmd/Quantum-Schrodinger-PINN`

- Dal: `claude/quantum-schrodinger-audit-iy6kn4`.
- Revizyon bilgisi: başlık bloğundaki tablo.
- Fiziksel şekiller `make_figures.py` ile `1c72630` commit'inde üretilmiştir.

**Ortam.** Kesin sürümler `requirements-lock.txt`'tedir: Python 3.13, NumPy 2.5.3, SciPy 1.18.1, PyTorch 2.15.0 (CPU).

**Komutlar** (README §11; **bu taslak için çalıştırılmamıştır**):

```bash
python -m pytest -q tests
python validate_solver.py
python experiment_disorder.py --sweep all
python experiment_localization.py --n-real 16
python experiment_pinn.py
python paired_analysis.py
python make_figures.py
```

**Köken bilgisi.** Her sonuç yapılandırmasını, tohumlarını, yapılandırma ve kod parmak izlerini, git commit'ini ve
kütüphane sürümlerini kaydeder (`provenance.py`). Rapora girecek sonuçlar için `--strict` önerilir (README §11).

**Tekrar üretim durumu** (`docs/REPRODUCIBILITY_AUDIT.md` §5):

| Sonuç | Durum |
|---|---|
| V1–V7, W taraması, ζ taraması, kapalı kutu | bit düzeyinde aynı |
| PINN `clean_fourier_uniform_s0` (2 iş parçacığı) | ağırlıklar bit düzeyinde aynı |
| momentum, length, barrier, convergence taramaları; diğer 17 PINN koşusu | yeniden çalıştırılmadı (**UNKNOWN**) [ÖĞRENCİ TARAFINDAN TAMAMLANACAK: `--strict` ile yeniden üretim] |

**Meta veri notu.** 44 transport `.npz` dosyasının 22'sinde `config_json` T = 1.0 yazar, gerçek t_end ise 1.25–2.5'tir.
Sayılar bundan etkilenmez (`docs/REPRODUCIBILITY_AUDIT.md` §2).

---



### Born oranı şeklinin açıklaması

**Şekil açıklaması tutarsızlığının incelenmesi.** `figures/CAPTIONS.md` ve aynı metni üreten
`make_figures.py::fig_born_ratio`, oranın zayıf saçılmada "~1.05'e yaklaştığını" söyler. Şekil `ratio_avg` alanını çizer.
`analysis.json`'daki saklı değerler:

| Alt küme | Nokta sayısı | `ratio_avg` |
|---|---|---|
| pürüzlü düzensizlik, kℓ_tr ≥ 10 (kζ < 2; W, ζ ve k taramaları) | 11 | 1.022–1.121; ortalama ve medyan ≈ 1.08 |
| bunların en zayıf saçılmalı üçü (kℓ_tr ≈ 42–167: W = 50, k_x0 = 30, W = 25) | 3 | 1.058, 1.062, 1.022 |
| yumuşak düzensizlik (kζ ≥ 2; ζ = 0.10, 0.15, 0.20) | 3 | 0.950, 0.965, 0.950 |

**Değerlendirme.** Açıklama yanlış değil, **eksik ve yuvarlanmıştır**:

- "~1.05" yalnızca pürüzlü düzensizlik alt kümesinin en zayıf saçılmalı noktalarını kabaca tanımlar.
- Yumuşak düzensizlik noktalarını ayrıca "oran < 1" diye anar.
- kℓ_tr ≥ 10'daki tam aralığı (0.95–1.12) vermez. Pürüzlü noktaların ortalaması (≈ 1.08) da 1.05'ten büyüktür.


Ana rapordaki Şekil 3 açıklaması veriden alınan aralıkları kullanır.
