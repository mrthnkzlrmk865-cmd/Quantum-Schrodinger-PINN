# Nihai Rapor Öncesi Bilimsel Durum ve Sonuç Envanteri (Rapor Hazırlık Denetimi)

**Denetimin türü.** READ-ONLY. Hiçbir çözücü, driver, referans hesabı, manifest üretimi veya test çalıştırılmadı.
Mevcut dosyaların hiçbiri değiştirilmedi. Bu dosya, denetim sırasında oluşturulan tek yeni dosyadır.

- **Tarih:** 2026-11-15.
- **Kanıt düzeyi.** Aşağıdaki her sayısal değer, belirtilen sonuç dosyasından okunarak doğrulandı. Okunamayan veya kaynağı
  bulunamayan değerler **UNKNOWN** olarak işaretlendi.
- **Yeniden hesaplananlar.** Yalnızca iki basit aritmetik kontrol yapıldı ve ikisi de ham CSV/JSON dosyalarından alındı:
  - PINN için Welch t testi;
  - Born oranlarının min–maks aralığı.

  Bu kontroller fiziksel bir hesap değildir.

---

## 1. Depo durumu

| öğe | değer | kaynak |
|---|---|---|
| aktif dal | `claude/quantum-schrodinger-audit-iy6kn4` | `git branch --show-current` |
| HEAD | `26174afa7826ecba0ce65ea615aef580837b3559` ("docs: Turkish technical project summary…") | `git rev-parse HEAD` |
| `origin/<dal>` | aynı (`26174af`) | `git fetch`; `git rev-parse` |
| `origin/main` | `dd0007a852f00731ec1c98e07667610c0de0908d`, PR #25 birleştirmesi; ebeveynleri `b0b3dac` ve `26174af` | `git log -1 origin/main` |
| main ile fark | main 1 commit ileride (yalnızca birleştirme commit'i); **dosya ağacı HEAD ile birebir aynı** (`git diff HEAD origin/main` boş) | `git diff --quiet` |
| working tree | temiz (denetimden önce 0 değişiklik) | `git status --porcelain` |
| önceki raporlarla karşılaştırma | önceki aşamalarda raporlanan commit'ler (`aa06437`, `b0b3dac`, `fac5f5e`, `1668540`, `26174af`) geçmişte sırasıyla mevcut. O aşamalardan sonra **yeni bir kod veya sonuç commit'i yok** | `git log --oneline -8` |

---

## 2. Projenin bilimsel amacı: iki ayrı araştırma hattı

### 2.1 Ana fiziksel araştırma hattı (zamana bağlı, TDSE)

**Kaynak:** `README.md` §1–§3, `docs/PHASE2_METHOD_SELECTION.md` (giriş).

> **Araştırma sorusu** (README §1). Gauss-korelasyonlu rastgele potansiyelin gücü W ve korelasyon uzunluğu ζ, 2B sert
> duvarlı bir dalga kılavuzundaki sonlu düzensiz bölgeden geçen kuantum dalga paketinin iletim ve yansımasını nasıl
> belirler? Gözlenen eğilimler pertürbatif (Born) taşınım kuramıyla tutarlı mı? Erişilebilir sistem boyutlarında
> lokalizasyon sonlu-boyut etkilerinden ayırt edilebilir mi?
>
> **İkincil (yöntemsel) soru.** Aynı düzensiz 2B TDSE için PINN çözümü, doğrulanmış klasik çözücüye göre ne kadar doğrudur
> ve onu ne sınırlar?

| Konu | Gerçek tanım | Kaynak |
|---|---|---|
| Model | iψ_t = −½∇²ψ + Vψ, ħ = m = 1, sert duvarlar (ψ = 0) | README §3; `reference_solver.py` |
| Dalga kılavuzu | L_y = 1, L_x = 2.9 + L_d; mod-1 Gauss paket, k_x0 = 20, σ_x = 0.2; uçlarda CAP; k ≈ 20.25, E₀ ≈ 205, **6 açık kanal** | README §3; `config.py` (`waveguide_config`) |
| Düzensizlik | sürekli Gauss-korelasyonlu alan, ⟨VV⟩ = W² e^{−r²/2ζ²}; M = 512 mod (Shinozuka–Deodatis); bölge 1.9 < x < 1.9 + L_d | README §4; `disorder.py` |
| Gözlemlenebilirler | T, R, P_in = 1 − T − R; t_end anında | README §8; `transport.py` |
| Klasik referans | DST-I Strang split-step (2. derece Δt, spektral uzay), CAP W₀ = 300 | README §5; `reference_solver.py` |
| PINN | 4×64 tanh MLP, Fourier öznitelikleri, sert IC/BC; Adam + L-BFGS | README §6; `model.py`, `pinn_core.py`, `train.py` |
| Kapalı kutu | L × L, L ∈ {1, 2}; MSD ve katılım oranı (PR) | README §8; `experiment_localization.py` |

**ENAQT / Lindblad / HEOM.**

- Eski sürümün "açık kuantum sistemi / ENAQT" motivasyonu modelle uyuşmuyordu (`docs/AUDIT.md` 2.4, [M]) ve **kaldırıldı**
  (README "Project history").
- Güncel modelde çevre, dephasing, Lindblad veya HEOM **yoktur**; dinamik üniterdir.
- Nihai raporda bu yöntemler **kullanılıyormuş gibi sunulmamalıdır**.

### 2.2 RGF / W-Z sayısal doğrulama hattı (zamandan bağımsız)

**Kaynak:** `docs/PHASE2_METHOD_SELECTION.md` §2 (Task B), §3 (Task C, B1–B4).

**Amaç:** TDSE sonuçlarını, **bağımsız bir sabit-enerji yöntemiyle** (mod uzayında RGF + Fisher–Lee) doğrulamak.

**Doğrulama zinciri:**

- **B1:** temiz kılavuz.
- **B2:** y'den bağımsız bariyer; 1B ODE ile karşılaştırma.
- **B3 / B3X / B3Y:** sabit düzensizlik realizasyonları.
- **B4:** durağan çözüm ile zamana bağlı çözümün karşılaştırması.

  Spec §3 kuralı: Phase 2 taşınım sonuçları **yalnızca B1–B4'ün hepsi geçerse** bağımsız olarak doğrulanmış sayılır
  (`docs/RGF_B2_VALIDATION.md` başlığı).

**W/Z alt hattı (`experiments/wz_rgf_diagnostic/`).**

- B1 konfigürasyon F'nin faz kriteri başarısızlığını teşhis etmek için başlatıldı.
- Ardından cebirsel olarak eşdeğer W/Z formülasyonunun sayısal davranışı incelendi.
- Bu hat **üretim RGF çözücüsünün yuvarlama davranışıyla ilgilidir.** Ana fiziksel sorunun (W, ζ, k, L etkisi) bir parçası
  değildir ve hiçbir taşınım sonucunu doğrudan etkilemez:
  - TDSE taşınım verileri dx = 1/64 ile **ayrı bir çözücüyle** üretildi.
  - B1-F, dx = 1/1024 ve L_d = 3 olan bir RGF konfigürasyonudur.

**İlişkinin doğru sunumu.** Fiziksel hat, bir doğrulama hattı tarafından desteklenmesi **planlanan** bir sonuç setidir.
Doğrulama hattı henüz bu desteği sağlamamıştır (B4 yapılmadı, B1-F FAIL). Nihai raporda iki hat **tek bir araştırma sorusu
olarak sunulmamalıdır.**

---

## 3. Ana fiziksel araştırma sonuçları

Bu bölümün ortak özellikleri:

- **Çözücü.** Tüm transport verileri `experiment_disorder.py` (DST çözücü) ile üretildi; ham veriler
  `results/transport/<sweep>/<point>.npz`, istatistikler `summary.csv` dosyalarında.
- **Güven aralıkları.** %95 bootstrap aralıklarıdır (`T_ci_boot_lo/hi`).
- **Eşleştirilmiş tasarım.** Realizasyon i her noktada seed i'yi kullanır. Bu nedenle bir sweep içindeki noktalar
  bağımsız değildir (README §8).
- **Denetimde doğrulananlar.** Aşağıdaki tüm değerler `summary.csv` ve `analysis.json` dosyalarından okunarak README ile
  karşılaştırıldı; **uyuşmazlık bulunmadı.**

### 3.1 Çözücü doğrulaması (V1–V7)

| Kontrol | Değer (kaynak: `results/validation/validation_summary.json`, `*.csv`) |
|---|---|
| V1 tam özmod | max \|hata\| = 4.58e-14 |
| V2 serbest Gauss | rel. L2 ≤ 9.39e-10 |
| V3 Δt (PINN problemi, W = E₀) | Δt = 1e-5'te 5.29e-7 |
| V3 Δt (dalga kılavuzu, W = 200) | Δt = 2e-4'te \|ΔT\| = 1.96e-5, \|ΔR\| = 2.25e-5; 2. derece |
| V4 uzaysal | dx = 1/64'te \|ΔT\| ≤ 1.4e-5 (README §7; JSON'daki V4 satırlarıyla tutarlı) |
| V5 CAP sahte geri dönüşü | k = 20'de 1.0e-4; **k = 12'de 1.3e-2, k = 10'da 3.1e-2** (README §7, `v5_cap_reflection.csv`) |
| V7 mod sayısı M | M'ye anlamlı bağımlılık yok |

- **Sınıflandırma: A (RAPORA HAZIR).** Sonuçlar ölçülmüş ve bit düzeyinde yeniden üretilmiştir
  (`docs/REPRODUCIBILITY_AUDIT.md` §5).
- **Uyarı.** CAP sapması, k_x0 = 12 noktasında R'yi yaklaşık 0.01 yukarı kaydırır.

### 3.2 Düzensizlik gücü W (`results/transport/strength/summary.csv`)

- **Parametreler:** W ∈ {0, 25, 50, 75, 100, 150, 200, 300}, ζ = 0.05, L_d = 1, k_x0 = 20, n = 48. W = 0 için n = 1
  (deterministik).
- **⟨T⟩ [95% CI]:**

  | W | ⟨T⟩ | 95% CI |
  |---|---|---|
  | 25 | 0.938 | [0.935, 0.941] |
  | 50 | 0.792 | [0.783, 0.801] |
  | 75 | 0.629 | [0.614, 0.645] |
  | 100 | 0.494 | [0.475, 0.513] |
  | 150 | 0.316 | [0.295, 0.337] |
  | 200 | 0.222 | [0.202, 0.242] |
  | 300 | 0.122 | [0.104, 0.141] |

- P_in(t_end) ≤ 3.2e-3.
- **Yorum.** ⟨T⟩, W ile monoton azalır. Bu bir **gözlemdir**; bu geometri ve paket için geçerlidir.
- **Sınıflandırma: A.**

### 3.3 Born + difüzyon tahmini ile karşılaştırma (`results/transport/analysis.json`)

- **Oran aralıkları.** kℓ_tr ≥ 10 olan noktalarda ⟨T⟩_sim/tahmin (enerji ortalamalı, `ratio_avg`) **0.950–1.121**
  aralığındadır; ortalamasız oran 0.907–1.131'dir. kℓ_tr ≈ 1–3'te oran 1.21–1.31'dir (W200: 1.213, W300: 1.309,
  k12: 1.303).
- **Yorum.** kℓ azaldıkça oran sistematik olarak artıyor. **Mekanizma belirlenmedi** (README §9.2). ℓ/(L+ℓ) biçimi bir
  mertebe tahminidir.
- **Sınıflandırma: B (SINIRLAMAYLA HAZIR).** Born oranları sınırsız 2B içindir; kanal yapısını içermez. π/2 öneki
  konvansiyona bağlıdır, bu yüzden yalnızca eğilim yorumlanmalıdır.

### 3.4 Korelasyon uzunluğu ζ (`results/transport/corr_length/summary.csv`; `results/analysis/paired_corr_length.json`)

- **Parametreler:** ζ ∈ {0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2}, W = 100, n = 48.
- **⟨T⟩ [95% CI]:**

  | ζ | kζ | ⟨T⟩ | 95% CI |
  |---|---|---|---|
  | 0.02 | 0.4 | 0.670 | [0.656, 0.684] |
  | 0.05 | 1.0 | 0.494 | [0.475, 0.513] |
  | 0.2 | 4.0 | 0.734 | [0.687, 0.777] |

- ⟨T⟩, ζ'ye göre monoton değildir ve minimumu kζ ≈ 1'dedir.
- sd(T), ζ ile artar: 0.050 → 0.159.
- **Eşleştirilmiş kontrastlar** (post-hoc; seçimi sonuçlar görüldükten sonra yapılmıştır):

  | kontrast | fark | %95 CI |
  |---|---|---|
  | C1 | −0.176 | [−0.198, −0.154] |
  | C2 | −0.240 | [−0.284, −0.193] |
  | C3 (eğrilik) | +0.046 | [0.023, 0.068] |

  Holm-düzeltilmiş bootstrap p = 1.5e-4, bu değer bootstrap çözünürlük sınırındadır.
- **Kaynak notu.** Paired analizin girdisi `results/reproduction/transport/corr_length` dizinidir, `results/transport`
  değildir. İkisinin bit düzeyinde aynı olduğu `docs/REPRODUCIBILITY_AUDIT.md` §5'te belgelenmiştir.
- **Sınıflandırma: B.**
  - Minimum yalnızca 7 noktalı ızgara çözünürlüğünde belirlenmiştir; ζ = 0.07 ile farkın düzeltilmiş aralığı sıfırı
    içerir (README §9.3).
  - Kontrastlar post-hoc'tur ve bağımsız bir tekrar yapılmamıştır.

### 3.5 Gelen momentum k (`results/transport/momentum/summary.csv`)

- **Parametreler:** k_x0 ∈ {12, 16, 20, 25, 30}, W = 100, n = 32; her nokta için temiz referans da var.
- **⟨T⟩:**

  | k_x0 | ⟨T⟩ | 95% CI |
  |---|---|---|
  | 12 | 0.251 | [0.222, 0.280] |
  | 16 | 0.360 | — |
  | 20 | 0.505 | — |
  | 25 | 0.691 | — |
  | 30 | 0.835 | [0.826, 0.843] |

- Temiz kılavuzda k12 için T = 0.984 ve R = 0.014'tür; bu R değeri CAP sapmasıdır.
- **Sınıflandırma: B.** k = 12 noktası CAP sapması nedeniyle yaklaşık 0.01 yanlıdır.

### 3.6 Bölge uzunluğu L ve lokalizasyon (`results/transport/length/summary.csv`, `analysis.json["length"]`)

- **Parametreler:** L_d ∈ {0.25, 0.5, 1, 2, 3}, W ∈ {100, 200}, n = 24.
- **⟨T⟩:**

  | W | L_d = 0.25 | 0.5 | 1 | 2 | 3 |
  |---|---|---|---|---|---|
  | 100 | 0.827 | 0.694 | 0.518 | 0.320 | 0.232 |
  | 200 | 0.602 | 0.414 | 0.249 | 0.104 | 0.056 |

- **Fitler** (`analysis.json["length"]`):

  | W | Ohm fiti: ℓ, χ²_red | Üstel fit: ξ, χ²_red |
  |---|---|---|
  | 100 | 1.000 ± 0.026, 3.60 | 2.01, 6.83 |
  | 200 | 0.242, 16.1 | 1.147 ± 0.038, 5.71 |

- **var(ln T)** L = 3'te: W = 100 için 0.042, W = 200 için 0.133.
- **Yorum** (README §9.5).
  - W = 100: difüzif biçim daha iyi uyar.
  - W = 200: lokalizasyonun başlangıcıyla **tutarlı, ama kanıt değil.**
    - L_d/ξ ≤ 2.6'dır.
    - var(ln T) ≪ |⟨ln T⟩|; enerji ortalaması dalgalanmaları bastırır.
    - Hiçbir fit χ²_red ≈ 1 vermez.
    - Alanın yaklaşık %15'i klasik olarak yasaktır (V > E₀).
- **Sınıflandırma: B.** Bu sweep'te n = 24, yani diğerlerinden düşük. "Lokalizasyon" ifadesi kullanılmamalı; yalnızca
  "üstel benzeri azalma" denmeli.

### 3.7 Kapalı kutu: yayılma ve sonlu-boyut testi (`results/localization/summary.csv`)

- **PR/alan** (L = 1 / L = 2; ± SE, n = 16):

  | W | L = 1 | L = 2 |
  |---|---|---|
  | 0 (n = 1) | 0.375 | 0.366 |
  | 100 | 0.415 ± 0.007 | 0.425 ± 0.022 |
  | 200 | 0.334 ± 0.020 | 0.312 ± 0.031 |
  | 300 | 0.267 ± 0.023 | 0.228 ± 0.032 |

- **Sonuç.** PR oranları 3.4–4.1'dir, yani alanla yaklaşık orantılıdır. Bu boyutlarda **2B lokalizasyon saptanmadı.**
- **Sınıflandırma: A** (negatif sonuç olarak). W = 300'deki düşüşün klasik bir katkısı vardır: alanın yaklaşık %25'inde
  V > E₀.

### 3.8 Düzensizlik içinde bariyer (`results/transport/barrier/summary.csv`)

- **Parametreler:** V_b ∈ {150, 250}, W ∈ {0, 50, 100, 200}, n = 24.
- **⟨T⟩:**

  | V_b | W = 0 | W = 50 | W = 100 | W = 200 |
  |---|---|---|---|---|
  | 150 | 0.708 | 0.587 | 0.402 | 0.203 |
  | 250 | 0.355 | 0.315 | 0.247 | 0.160 |

- **Sonuç.** Düzensizliğin yardımıyla tünelleme gözlenmedi; bu sonuç bu parametrelerle sınırlıdır.
- **Sınıflandırma: A** (sınırlı kapsamlı bir negatif sonuç).

### 3.9 Realizasyon yakınsaması (`results/transport/convergence/summary.csv`)

- n = 192: ⟨T⟩ = 0.509 [0.500, 0.519], SE = 0.0048.
- 48 realizasyonluk ortalama 0.494'tür ve tam ortalamayla 1.5 SE içinde uyumludur.
- **Sınıflandırma: A.**

### 3.10 PINN ile klasik çözücü karşılaştırması (`results/pinn/summary.csv`)

- **Kurulum.** 18 koşu, yalnızca CPU (GPU yok); her koşu aynı potansiyelle DST referansına karşı değerlendirildi.
  Konfigürasyonlar: W ∈ {0, 28, 56, 112}, ζ = 0.08 (W = 56 için {0.04, 0.16} de var), T = 0.02.
- **t = T anındaki rel. L2(ψ):**

  | konfigürasyon | rel. L2(ψ) |
  |---|---|
  | temiz, 3 seed | 0.4885, 0.4867, 0.5665 |
  | W = 56, 3 seed | 0.5864, 0.5744, 0.6100 |
  | W = 112 | 0.7133 |
  | Fourier özniteliği olmadan, temiz | 1.152 |
  | T = 0.04 | 0.986 |
  | T = 0.01 | 0.184 |
  | 3 kat bütçe | 0.281 |

- **Maliyet.** Eğitim yaklaşık 470–3100 s, referans çözücü 0.7–2.9 s.
- **Welch testi** (denetimde yeniden hesaplandı): temiz ile W56 arasında t = 2.70, p = 0.086 (n = 3 + 3). README ile aynı
  sonuç.
- **Yorum.**
  - **PINN bu bütçede taşınım fiziği için yeterince doğru değildir.** En kolay ayarda bile son hata yaklaşık %50'dir;
    klasik çözücünün hatası yaklaşık 1e-6'dır.
  - Hata zaman penceresiyle büyür.
  - Düzensizlik etkisi bir eğilimdir ve %5 düzeyinde **anlamlı değildir.**
- **Sınıflandırma: B.** Sonuç yalnızca bu CPU bütçesi, tek mimari ve en fazla 3 seed için geçerlidir.
  - PINN eğitimi PyTorch iş parçacığı sayısına bağlıdır: 2 iş parçacığı bit düzeyinde yeniden üretilebilir; 4 iş
    parçacığında son hata 0.48850'den 0.48987'ye çıkar (`docs/REPRODUCIBILITY_AUDIT.md` §5).
  - PINN dalga kılavuzu geometrisinde **çalıştırılmamıştır**; yalnızca 1×1 kutu problemi incelenmiştir.

### 3.11 Tekrar üretilebilirlik (`docs/REPRODUCIBILITY_AUDIT.md`)

**Bit düzeyinde yeniden üretildi:**

- V1–V7;
- strength (8 × 48) ve corr_length (7 × 48) sweep'leri;
- localization (8 × 16);
- 2 iş parçacığıyla PINN `clean_fourier_uniform_s0`.

**Yeniden çalıştırılmadı:** momentum, length, barrier ve convergence sweep'leri ile diğer 17 PINN koşusu. Bunlar için
tekrar üretilebilirlik **UNKNOWN**'dur. Altyapı (provenance ve önbellek) bunu mümkün kılar ama yapılmamıştır.

**Meta veri kusuru.** 44 transport `.npz` dosyasının 22'sinde `config_json` T = 1.0 yazıyor, oysa gerçek t_end
1.25–2.5'tir. Bilimsel sayılar bundan etkilenmez, çünkü t_end `times` dizisinden alınır (Audit §2).

---

## 4. RGF ve W/Z araştırmaları

| İş | Sonuç dosyası | Kriter | Durum (aynen) | Ölçülen | Test ETMEDİĞİ | Sınıf |
|---|---|---|---|---|---|---|
| **B1** uzun cihaz (6 konfig.) | `results/rgf_b1_validation/summary.json`, `phase_audit_hp_summary.json`; `docs/RGF_B1_VALIDATION.md` (`bce6587`, denetim `d5a97a3`) | faz ≤ 1e-10 (geçici kriter, değiştirilmedi) | A–E geçti; **B1-F FAIL** | F: 1.76851694417934e-10; 100 basamaklı referansa göre 1.5424928001550615e-10 | sıfırdan farklı potansiyel; fiziksel doğrulama (B1 bir iç tutarlılık testidir) | B |
| **B2** bariyer, RGF vs 1B ODE | `results/rgf_b2_validation/summary.json`; `docs/RGF_B2_VALIDATION.md` (`bc8f144`) | \|ΔT\| ≤ 2e-4 ve diğerleri | **tüm bağlayıcı kriterler PASS** | max \|ΔT\| = 3.59e-5 / 5.95e-5; TD vs ODE 3.2e-4 / 7.0e-4 (≤ 1e-3) | kanal bağlaşımı (y'ye bağımlı V) | A |
| **B3** sabit realizasyon | `results/rgf_b3_validation/summary.json` (`d3f9291`) | 10 kriter grubu | **FAIL: dx_1_512_vs_1_1024** (9/10 geçti) | S uyumu 4.9e-12; dx: 5.97e-4 > 2e-4 | toplu istatistik | B |
| **B3X** 6 konfigürasyon | `results/rgf_b3x_validation/summary.json` (`5d907cd`) | 14 kriter | **FAIL: C11_dx_refinement** (13/14) | S uyumu ≤ 5.3e-12; C11: 9.17e-4 | süreklilik limiti, istatistik | B |
| **B3Y** ince ızgara | `results/rgf_b3y_validation/summary.json` (sürücü `1a5688b`, sonuç `36055a7`) | max\|ΔS\| ≤ 1e-10 ve diğerleri | **PASS** (`overall_status`; `stopped_at = null`) | 4 yöntem ≤ 9.1e-11; dx değişimleri 4 yöntemde de aynı → sonlu-fark kafesinin özelliği | sürekli denkleme yakınsama; başka seed ve enerjiler | B |
| **B4** durağan vs TD (sabit realizasyon) | **yok** | `PHASE2_METHOD_SELECTION.md` §3 B4 | **YAPILMADI** (README: "B4 open") | — | — | C |
| **W/Z-1** scalar, tek nokta | `experiments/wz_rgf_diagnostic/REPORT.md`, `results/worst_point.json`; `ERRATUM_worst_point.md` | önceden kayıtlı kural (`PREDICTIONS.md`, `3ce732d`) | **INCONCLUSIVE** | faz −1.19e-15 (orijinal +1.5425e-10); genlik +4.82e-14 (orijinal −1.08e-15) | yalnızca 1×1 kanal-5; kapalı kanallar çalışmadı (erratum E1–E3) | B |
| **W/Z-2** 63 kanal | `amplitude_roundoff_report.md`, `results_amplitude_roundoff/amplitude_roundoff_results.json` (koşu 3: `e274064`) | `PREDICTIONS_amplitude_roundoff.md` | H_A: **SUPPORTED, kanal 5'te bağımsız olarak ayırt edilemedi** | W/Z faz ≤ 3.24e-15; orijinal 1.0e-11–1.54e-10; W/Z genlik max 1.93e-13 (üretim 5.7e-14; ≈3.4×); \|P\| − 1 = 1.37e-17 → 4.757e-14 | bağlı bloklar, düzensizlik, N-ölçeklenmesi; tek enerji | B |
| **W/Z-3** coupled-block | `results_coupled_blocks/coupled_blocks_report.md` ve `.json`/`.csv` (koşu `63e046d`) | `PREDICTIONS_coupled_blocks.md` §8 | **DESIGN INSUFFICIENT** (W/Z FAIL değil) | 76 execution / 72 konfigürasyon; referans ≤ 2.0e-61 / 3.2e-61; uç ≤ 2.5e-16; 480/480 mutant; N8 eğimi 0.78–0.79 < 0.8 | uzun cihaz, eşik yakını, düzensizlik ansamblı | B |
| **E** coupling-scaling v2 tasarımı ve manifestleri | `PREDICTIONS_coupling_scaling_v2.md` (Rev 3.7, `0f502c9`); `results_coupling_scaling_v2/manifest_{XP,HO}.json`, `manifest.sha256` (`fac5f5e`) | §2.4, §8.2 | manifestler **GENERATED AND FROZEN**; ana deney **NOT RUN** | XP 60 (24/24/12), HO 120 (48/48/24); redraw 11 / 10 (hepsi enerji); max binary64 cond 382.1 | — (dizinde **yalnızca** üç manifest dosyası var; results/points/pairs dosyası yok) | C |
| **F** T1–T6 yazılım doğrulaması | resmi koşu kaydı (`1668540`); `R0_REVIEW_PACKET.md` I.7 | bayraksız resmi koşu | **PASS** (yazılım/protokol) | 54/507/28/468/3528/36, 0 başarısız; C4 doğrulayıcı değil; pytest 691 passed + 2 xfailed; altyapı testleri 122 | fiziksel doğruluk; bağımsız inceleme | D |
| **G** uygulama düzeyi R0 | `R0_REVIEW_PACKET.md` (H boş); `REVIEW_coupling_scaling_v2.md` **yok** | spec §4.3 | **PENDING** | — | — | C/D |

**Notlar.**

1. **B1-F** FAIL olarak kalır. Kriter değiştirilmedi. Bu belge B1-F için bir mekanizma iddiası eklemez; yalnızca
   `REPORT.md` §4'teki tek noktalı tam ayrıştırmayı aktarır: temsil +1.6108e-10, yaklaşık %104.
2. **W/Z genel doğrulaması: NOT VALIDATED.** Hiçbir deney bunu göstermez.
3. **Coupled-block DESIGN INSUFFICIENT** sonucu, önceden kayıtlı tasarımın karar üretememesidir. Bu W/Z'nin başarısız
   olduğu anlamına gelmez (rapor §9.1).
4. **Ana coupling-scaling deneyi çalıştırılmadı.** Driver kapılıdır: G-1 "duplicate set" tanımsız, R0 imzası yok,
   `previous_files.sha256` yok.
5. **Testlerin yeniden çalıştırılması.** Bu denetimde yeniden çalıştırılmadı (solver çalıştırma yasağı). En son resmi
   koşu, yalnızca belge commit'inden (`26174af`) önceki `1668540` üzerindedir.

---

## 5. Raporlaştırma sınıflandırması (özet)

| Bulgu | Sınıf | Gerekçe |
|---|---|---|
| DST çözücü doğrulaması (V1–V7) | **A** | ölçülmüş, bit düzeyinde yeniden üretilmiş; düşük k'da CAP sapması belirtilmeli |
| ⟨T⟩(W) monoton azalma | **A** | n = 48, CI'lar ayrık, yeniden üretilmiş |
| Kapalı kutuda 2B lokalizasyon saptanmaması | **A** | negatif sonuç; boyut sınırı açık |
| Bariyer + düzensizlik, tünelleme yok | **A** | sınırlı parametrelerde negatif sonuç |
| Realizasyon yakınsaması | **A** | doğrudan ölçüm |
| ⟨T⟩(ζ) minimumu kζ ≈ 1'de | **B** | ızgara çözünürlüğü kaba; eşleştirilmiş analiz post-hoc |
| Born oranının kℓ ile artması | **B** | mekanizma belirsiz; tahmin bir mertebe formudur |
| ⟨T⟩(k) artışı | **B** | k = 12'de CAP sapması; tekrar üretilmedi |
| L-bağımlılığı; W = 200'de üstel benzeri azalma | **B** | n = 24; L_d/ξ ≤ 2.6; χ²_red > 1; lokalizasyon kanıtı değil; tekrar üretilmedi |
| PINN başarısızlık modları | **B** | bütçe, tek mimari, ≤ 3 seed; iş parçacığı bağımlılığı |
| B1 (A–E geçti, F FAIL) | **B** | iç tutarlılık testi; F FAIL olarak raporlanmalı |
| B2 PASS | **A** | bağımsız ODE'ye karşı, kanal bağlaşımı yok (kapsam açıkça belirtilmeli) |
| B3 / B3X FAIL (dx), B3Y PASS | **B** | B3Y, dx değişiminin kafes özelliği olduğunu gösterir; süreklilik limiti açık |
| B4 (TD vs durağan) | **C** | tasarlandı, sonuç yok |
| W/Z scalar (INCONCLUSIVE) ve 63 kanal | **B** | tek nokta, temiz cihaz; erratum düzeltmeleri esas alınmalı |
| Coupled-block (DESIGN INSUFFICIENT) | **B** | sonuç gerçek; karar üretmedi |
| Coupling-scaling v2 ana deneyi | **C** | manifest var, sonuç yok |
| T1–T6 / pytest PASS | **D** | yalnızca yazılım/protokol; bağımsız inceleme yok |
| Clean-room transcription, oracle, sertifika kodu | **D** | aynı yapay zekâ ekosistemi; R0 PENDING |
| Planlanan fiziksel kampanya (`PHYSICAL_CAMPAIGN_PLAN.md` rev. 3) | **C** | dondurulmadı, hash'lenmedi, çalıştırılmadı; kök entropi ve C0 seed kaydedilmedi (`CAMPAIGN_STATUS_CONTRACT.md`) |
| B1 için "korelasyon 0.91–0.93, eğim 0.99–1.05" | **E** | hiçbir commit edilmiş dosyada yok (önceki görev sırasında aranıp bulunamadı) |
| "Disorder increases transmission (0.424 vs 0.277)" (eski README) | **E** | `docs/AUDIT.md` 2.3'te çürütüldü; README §9.9 düzeltmesi esas alınmalı |
| ENAQT / fotosentez / Lindblad / HEOM çerçevesi | **E** | model bunları uygulamıyor (`AUDIT.md` 2.4) |
| Eski README'nin 1B/periyodik model iddiası | **E** | `AUDIT.md` 1.1: kod 2B ve Dirichlet |

---

## 6. Rapor için kullanılabilecek grafikler

Ana fiziksel hat grafiklerinin tamamı `make_figures.py` ile `1c72630` (2026-11-05) commit'inde üretildi. Açıklamaları
`figures/CAPTIONS.md` dosyasındadır.

| Dosya | Gösterdiği | Kaynak veri | Parametreler | Uygunluk | Eksik açıklama |
|---|---|---|---|---|---|
| `figures/transport_vs_strength.png` | T, R, P_in vs W; Born tahmini | `results/transport/strength/summary.csv` | ζ = 0.05, L_d = 1, k = 20, n = 48 | uygun | eşleştirilmiş tasarım nedeniyle noktalar bağımsız değil |
| `figures/transport_vs_correlation_length.png` | T, R vs kζ; Born ℓ_tr | `corr_length/summary.csv` | W = 100, n = 48 | uygun | minimumun kaba çözünürlüğü |
| `figures/paired_corr_length.png` | eşleştirilmiş kontrastlar | `results/analysis/paired_corr_length.json` | 48 seed, B = 20000 | uygun (post-hoc etiketiyle) | veri yolu `results/reproduction/...` |
| `figures/transport_vs_momentum.png` | T vs k; temiz değerler | `momentum/summary.csv` | W = 100, n = 32 | uygun | k = 12'de CAP sapması |
| `figures/transport_vs_length.png` | ⟨T⟩ vs L_d; Ohm ve üstel fitler | `length/summary.csv`, `analysis.json` | W = 100/200, n = 24 | uygun | "lokalizasyon kanıtı değil" uyarısı şart |
| `figures/born_ratio_vs_kl.png` | sim/tahmin vs kℓ_tr | `transport/analysis.json` | üç sweep | uygun | açıklamada "oran ≈ 1.05'e yaklaşır"; README'de 0.95–1.12 aralığı (bkz. §7 çelişki K-2) |
| `figures/transport_barrier.png` | bariyer + düzensizlik | `barrier/summary.csv` | V_b 150 ve 250, n = 24 | uygun | — |
| `figures/realization_convergence.png` | n ile yakınsama | `convergence/` | n = 192 | uygun | — |
| `figures/localization_closed_box.png` | PR/alan vs t | `results/localization/` | L = 1/2, n = 16 | uygun | klasik yasak bölge katkısı |
| `figures/validation_solver.png` | V2–V4 yakınsama | `results/validation/*.csv` | — | uygun | — |
| `figures/disorder_examples.png` | V_dis örnekleri | `disorder.py` ile deterministik (seed 0) | W = 100, üç ζ | uygun (illüstrasyon) | — |
| `figures/waveguide_snapshots.png` | \|ψ\|² anlık görüntüleri | **saklı veri yok**; çözücü `make_figures.py` içinde çalıştırılıyor | W = 150, seed 0, dx = 1/64 | yalnızca illüstrasyon | tek realizasyon; veri dosyası yok |
| `figures/pinn_errors.png` | tüm PINN koşularının son hataları | `results/pinn/summary.csv` | 18 koşu | uygun | ≤ 3 seed; iş parçacığı bağımlılığı |
| `figures/pinn_error_vs_time.png` | hata vs t | `results/pinn/<run>/eval_timeseries.csv` | — | uygun | — |
| `figures/pinn_vs_reference_W56.png` | PINN vs referans | `results/pinn/W56_fourier_uniform_s0/` | W = 56, ζ = 0.08 | uygun (tek örnek) | temsili değil, tek koşu |
| `results/rgf_b2_validation/{error_vs_dx,error_vs_energy,order_vs_energy}.png` | B2 hata ve yakınsama | `results/rgf_b2_validation/` | B2 | ek / ekler için uygun | — |
| `docs/turkce_ozet_figurler/sekil{1,2,3}_*.png` | B1 faz; 63 kanal ε_φ/ε_A; coupled-block e_glob | `phase_audit_hp_summary.json`; `amplitude_roundoff_results.json`; `coupled_blocks_comparisons.csv` | — | RGF/W-Z bölümü için uygun | — |

Yeni grafik üretilmedi.

---

## 7. Bilimsel ve matematiksel eksikler, çelişkiler

### 7.1 Eksikler

| # | Eksik | Neden önemli | Rapordan önce çözülmeli mi? | Çözülmezse raporda nasıl açıklanmalı |
|---|---|---|---|---|
| M-1 | **B4 yok** (TD vs durağan, sabit realizasyon) | TDSE taşınım sonuçlarının bağımsız sabit-enerji doğrulaması tamamlanmadı | Hayır, ancak "bağımsız doğrulandı" iddiası için **evet** | "TDSE sonuçları yakınsama ve analitik testlerle doğrulandı; bağımsız durağan çapraz kontrol (B4) yapılmadı." |
| M-2 | Lokalizasyon için yetersiz boyut (L_d ≤ 3, 6 kanal, kutu ≤ 2×2) | W = 200 yorumu ve 2B negatif sonuç boyutla sınırlı | Hayır | README §9.5–§9.6 ve §10 dilini koru: "tutarlı, kanıt değil" |
| M-3 | Enerji ortalaması (σ_k ≈ 3.5) | ln T dağılımı ve sabit-enerji niceliklerine erişim yok | Hayır | sınırlamalarda açıkça belirt |
| M-4 | Örnek sayısı: length/barrier için n = 24, momentum için n = 32 | farklar ≳ 0.03 çözülebilir; ince eğilimler belirsiz | Hayır | CI'ları göster; nokta bazlı CI olduklarını belirt |
| M-5 | Eşleştirilmiş tasarım: noktalar bağımsız değil | çoklu karşılaştırma ve eğri fitleri için | Hayır | README §8 dilini koru; fit χ²'lerinin korelasyonu yok saydığını belirt |
| M-6 | Born sapmasının mekanizması belirsiz | ana bulgulardan biri açıklanmamış | Hayır | "gözlem; aday nedenler; test edilmedi" |
| M-7 | Kısmi tekrar üretim: momentum, length, barrier, convergence ve 17 PINN koşusu yeniden çalıştırılmadı | nihai rapor için `--strict` koşusu önerilmişti (README §11) | **Önerilir** (çözücü çalıştırmak gerekir; bu denetimin kapsamı dışında) | "yalnızca X sweep'leri bit düzeyinde yeniden üretildi" |
| M-8 | 22/44 `.npz` dosyasında eski t_end meta verisi | tekrar üretimde yanıltıcı olabilir | Hayır (sayılar etkilenmez) | ekte belirt |
| M-9 | PINN: ≤ 3 seed, CPU bütçesi, tek mimari, iş parçacığı bağımlılığı; dalga kılavuzu geometrisinde yok | genelleme yapılamaz | Hayır | "bu bütçede PINN yetersiz kaldı; genel PINN iddiası yok" |
| M-10 | B1-F FAIL (geçici kriter) | RGF uzun cihazlarda faz doğruluğu sınırında | Hayır | FAIL'i aynen raporla; TDSE verilerini etkilemediğini açıkla |
| M-11 | B3/B3X dx kriteri FAIL; süreklilik limitine yakınsama bilinmiyor | RGF transmisyonunun dx = 1/512'de 2e-4 düzeyinde yakınsamadığı (W = 200) | Hayır | B3Y'nin bunu kafes özelliği olarak ayırdığını belirt |
| M-12 | W/Z bağımsız doğrulama eksik; R0 PENDING; G-1 tanımsız; A-1 açık | W/Z hakkında hiçbir genel sonuç verilemez | W/Z sonucu iddia edilecekse **evet** | "W/Z NOT VALIDATED; ana deney NOT RUN" |
| M-13 | T1–T6 yalnızca yazılım testi; aynı yazar ekosistemi | yazılım geçişi fiziksel doğrulama değil | Hayır | "software validation only" ayrımını koru |
| M-14 | Fiziksel kampanya (rev. 3) dondurulmadı ve çalıştırılmadı; kök entropi ve C0 seed yok | ζ etkisinin önceden kayıtlı doğrulayıcı testi yok | ζ hakkında doğrulayıcı iddia için **evet** | "ζ bulguları keşif niteliğinde; doğrulayıcı kampanya planlandı, yapılmadı" |
| M-15 | Matematiksel türetimler: Born oranları ve ℓ/(L+ℓ) önekleri konvansiyona bağlı; W/Z türetimi ve oracle formülleri bağımsız olarak doğrulanmadı | raporun matematiksel bölümünün sağlamlığı | Hayır (fizik için); W/Z için R0 | önek bağımlılığını ve R0 durumunu açıkça belirt |

### 7.2 Çelişkiler (uzlaştırılmadan raporlanıyor)

| # | Çelişki | Kaynak A | Kaynak B |
|---|---|---|---|
| K-1 | README'nin depo yapısı satırı (satır 596) "B3 and B3X: dx criterion fails; B4 open" diyor ve **B3Y'den (PASS) hiç söz etmiyor**; ayrıca README'nin `docs/` listesinde B3Y yok | `README.md` §12 | `docs/RGF_B3Y_VALIDATION.md` (overall PASS); `results/rgf_b3y_validation/summary.json` |
| K-2 | Born oranı: grafik açıklaması "zayıf saçılmada ~1.05'e yaklaşır" diyor; README "kℓ_tr ≳ 10 için 0.95–1.12" diyor | `figures/CAPTIONS.md` (`born_ratio_vs_kl.png`) | `README.md` §9.2; `analysis.json` (ratio_avg 0.950–1.121) |
| K-3 | Dondurulmuş spec başlığı "nothing … implemented or run; no matrix generated" diyor; oysa altyapı yazıldı ve manifestler üretildi | `PREDICTIONS_coupling_scaling_v2.md` başlığı | `fac5f5e`; `R0_REVIEW_PACKET.md` Part I |
| K-4 | `docs/rgf_b3y_protocol.json` içindeki `status` metni "no B3Y numerical comparison has been run" diyor (tarihsel), oysa kampanya tamamlandı | protokol JSON | `RGF_B3Y_VALIDATION.md` "Status (reconciled)" tablosu; bu belge çelişkiyi zaten açıklıyor |
| K-5 | B1 için görev metninde geçen korelasyon/eğim değerleri | önceki görev metni | depoda kaynak yok → **E** |
| K-6 | Paired analizin girdi dizini `results/reproduction/...`; README ise analizi `results/transport` sonuçları üzerinden sunuyor | `paired_corr_length.json` `data` alanı | `REPRODUCIBILITY_AUDIT.md` §5: bit düzeyinde aynı; çelişki yok, ama raporda belirtilmeli |

---

## 8. Nihai rapor için önerilen yapı

**İlke.** İki hat **ayrı sorular** olarak sunulmalıdır.

- **Bölüm I** fiziksel çalışmadır.
- **Bölüm II**, Bölüm I'in sayısal güvenilirliğini destekleme amacıyla başlatılan ve henüz tamamlanmamış bir doğrulama
  çalışmasıdır.
- W/Z, Bölüm II'nin içinde bir **teşhis alt çalışması** olarak yer almalıdır.

**Önerilen başlıklar:**

1. Özet. Ana bulgular ve sınırlamalar; durum tablosu.
2. Giriş ve araştırma sorusu. Fiziksel soru ve ikincil yöntemsel soru (PINN).
3. Fiziksel model. TDSE, dalga kılavuzu geometrisi, Gauss-korelasyonlu düzensizlik, Born kıyas formülleri (kuram olarak
   işaretlenmiş).
4. Sayısal yöntem ve doğrulama. DST split-step, CAP, V1–V7; CAP sapması; tekrar üretilebilirlik durumu.
5. Sonuçlar.
   - 5.1 W
   - 5.2 Born karşılaştırması
   - 5.3 ζ, eşleştirilmiş analiz dahil (post-hoc etiketiyle)
   - 5.4 k
   - 5.5 L ve lokalizasyon soruları ("kanıt değil")
   - 5.6 kapalı kutu, negatif sonuç
   - 5.7 bariyer
   - 5.8 realizasyon istatistikleri
6. PINN değerlendirmesi. Yöntemsel sonuç: bu bütçede yetersiz; başarısızlık modları; istatistiksel güç sınırı.
7. Bağımsız sabit-enerji doğrulama hattı (RGF).
   - Motivasyon (Phase 2).
   - Durum: B1 (F FAIL), B2 PASS, B3/B3X dx FAIL, B3Y PASS, **B4 yapılmadı**.
   - Sonuç: fiziksel sonuçlar **henüz bağımsız olarak doğrulanmadı**.
8. RGF yuvarlama teşhisi ve W/Z (ek veya ayrı bir teknik not olabilir).
   - B1-F ayrıştırması;
   - scalar ve 63-kanal deneyleri (erratum dahil);
   - coupled-block (DESIGN INSUFFICIENT);
   - coupling-scaling v2 tasarımı ve manifestleri (NOT RUN);
   - R0 PENDING.
9. Sınırlamalar. README §10, ayrıca M-1 … M-15.
10. Sonuç ve açık sorular. Kuram ile proje sonucu ayrımı (README §9.11 tablosu).
11. Tekrar üretilebilirlik ve veri erişimi. Commit'ler, sonuç dosyaları, `--strict` önerisi.
12. Kaynaklar. `docs/LITERATURE.md`'deki doğrulanmış liste.

**Bilimsel katkı (abartısız).** Projenin katkısı şudur:

- Bu belirli geometri (6 kanallı kılavuz, Gauss-korelasyonlu düzensizlik) için ⟨T⟩(W, ζ, k, L_d)'nin güven aralıklı
  sayısal haritası;
- Born + difüzyon tahminine göre kℓ'ye bağlı sistematik sapmanın gözlenmesi;
- kζ ≈ 1 civarında bir ⟨T⟩ minimumu;
- erişilebilir boyutlarda 2B lokalizasyonun saptanmaması;
- PINN'in bu bütçede taşınım için yetersiz kaldığının nicel gösterimi;
- üretim RGF çözücüsünün uzun cihazlardaki faz yuvarlamasının belgelenmesi.

**Yeni kuram, yeni yöntem veya lokalizasyon keşfi iddia edilmemelidir.**

---

## 9. Kontrol edilen dosyalar

- **Belgeler:**
  - `README.md` (tamamı);
  - `docs/AUDIT.md`, `REPRODUCIBILITY_AUDIT.md`, `PHASE2_METHOD_SELECTION.md` (§1–§3);
  - `ODE1D_VALIDATION.md`, `RGF_B1/B2/B3/B3X/B3Y_VALIDATION.md`, `RGF_B3Y_TIER3_AUDIT.md`;
  - `PHYSICAL_CAMPAIGN_PLAN.md`, `CAMPAIGN_STATUS_CONTRACT.md`, `LITERATURE.md`;
  - `figures/CAPTIONS.md`.
- **Sonuç dosyaları:**
  - `results/transport/*/summary.csv` ve `results/transport/analysis.json`;
  - `results/analysis/paired_corr_length.json`;
  - `results/localization/summary.csv`, `results/pinn/summary.csv`, `results/validation/validation_summary.json`;
  - `results/ode1d_validation/summary.json`;
  - `results/rgf_b{1,2,3,3x,3y}_validation/summary.json`.
- **W/Z hattı:** `experiments/wz_rgf_diagnostic/` altındaki raporlar, sonuçlar ve manifest dizini.
- **Kod:** `make_figures.py` (kaynak eşlemesi).
- **Git:** dal, HEAD, `origin/main` ve son commit'ler.

---

## 10. Güncel durum (değiştirilmeden)

| Konu | Durum |
|---|---|
| Ana TDSE taşınım sonuçları | mevcut, kısmen bit düzeyinde yeniden üretilmiş; **bağımsız durağan doğrulama (B4) yok** |
| PINN | bu bütçede taşınım için **yetersiz** |
| B1-F | **FAIL** |
| B2 | PASS |
| B3 / B3X | FAIL (dx kriteri) |
| B3Y | PASS |
| B4 | **yapılmadı** |
| Coupled-block W/Z | **DESIGN INSUFFICIENT** |
| W/Z genel doğrulaması | **NOT VALIDATED** |
| Main coupling-scaling experiment | **NOT RUN** |
| Implementation-level R0 | **PENDING** |
| XP/HO manifestleri | **GENERATED AND FROZEN** (`fac5f5e`) |
| T1–T6 yazılım doğrulaması | PASS (yalnızca yazılım/protokol) |
| Fiziksel doğrulayıcı kampanya | önerildi (rev. 3); **dondurulmadı, çalıştırılmadı** |
