# TÜBİTAK 2204-A Görselleştirme Notları (3B Plotly yüzey gösterimleri)

> **Kapsam.** Bu belge, `visualizations/tubitak/` altındaki görselleri ve bunların rapora nasıl entegre edileceğini
> açıklar. Resmî rapor henüz düzenlenmemiştir.
>
> - Bu görseller **yalnızca mevcut proje verisinden** üretilmiştir.
> - Hiçbir simülasyon çalıştırılmamış, hiçbir PINN eğitilmemiş veya yeniden değerlendirilmemiştir.
> - Hiçbir sayısal sonuç değiştirilmemiştir.
> - Fiziksel model **iki boyutludur**. 3B yüzey, 2B bir alanın (x, y) üzerindeki değerinin yükseklik olarak
>   gösterilmesidir; üç boyutlu bir fiziksel simülasyon değildir.

---

## 1. Oluşturulan dosyalar

| Dosya | İçerik |
|---|---|
| `visualizations/tubitak/fig1_duzensizlik_potansiyeli_3b.{png,svg,html}` | Şekil 1: V(x, y) |
| `visualizations/tubitak/fig2_pinn_dst_karsilastirma_3b.{png,svg,html}` | Şekil 2: DST ve PINN \|ψ\|², yan yana |
| `visualizations/tubitak/fig3_pinn_dst_mutlak_fark_3b.{png,svg,html}` | Şekil 3: \| \|ψ_PINN\|² − \|ψ_DST\|² \| |
| `visualizations/tubitak/wavefunction_animation_white.html` | beyaz arka planlı, iyileştirilmiş animasyon (yeni sürüm) |
| `visualizations/tubitak/preview_wavefunction_animation_white_t_end*.png` | animasyonun t = T karesinin önizlemeleri (tam ve yakın z ölçeği) |
| `visualizations/tubitak/figure_values.json` | şekillerde kullanılan değerlerin kaydı (tepe değerleri, hatalar, ızgaralar) |
| `visualizations/tubitak/make_tubitak_3d.py` | HTML şekilleri ve animasyonu üreten betik (yalnızca veri okur) |
| `visualizations/tubitak/export_static.js` | HTML'den PNG/SVG dışa aktarımı (headless Chromium + `Plotly.toImage`) |

**Değiştirilmeyen dosyalar.** Kök dizindeki özgün `wavefunction_animation.html` (ve onu açan `index.html`) aynen
korunmuştur; özgün sürüm olarak kullanılmaya devam eder. Hiçbir veri, sonuç, şekil veya kaynak kod dosyası
değiştirilmemiştir.

**Dosya özellikleri:**

- **PNG:** 3 kat ölçek; Şekil 1 ve 3 için 2700 × 2280 px, Şekil 2 için 4500 × 2280 px. Beyaz arka plan.
- **SVG:** Plotly, 3B (WebGL) sahneyi SVG içine **raster görüntü** olarak gömer. Yalnızca başlık, renk çubuğu ve
  açıklamalar vektördür. Baskı için PNG önerilir.
- **HTML:** Plotly.js dosyanın içine gömülüdür (5–8 MB). İnternet bağlantısı olmadan da açılır; sergi standı için
  uygundur.
- **Yazı tipi:** "Times New Roman"; sistemde yoksa metrik olarak eşdeğer Liberation Serif kullanılır. PNG'ler Liberation
  Serif ile üretilmiştir.
- **Ondalık ayırıcı:** Türkçe virgül.

**Yeniden üretim** (depo kökünden):

```bash
python visualizations/tubitak/make_tubitak_3d.py
node visualizations/tubitak/export_static.js
```

---

## 2. Veri kökeni (Phase 1)

### 2.1 Kaynak dosyalar ve diziler

| Şekil | Kaynak dosya | Dizi / alan | Not |
|---|---|---|---|
| Şekil 1 | `results/pinn/W56_fourier_uniform_s0/eval_fields.npz` | `V` (127 × 127), `X`, `Y` | — |
| Şekil 2 | aynı dosya | `psi_ref[-1]`, `psi_pinn[-1]` (karmaşık, 127 × 127) | gösterilen ρ = \|ψ\|² |
| Şekil 3 | aynı dosya | \|ψ_pinn[-1]\|² − \|ψ_ref[-1]\|² farkının mutlak değeri | — |
| Hata değerleri | `results/pinn/W56_fourier_uniform_s0/eval_summary.json` | `rel_l2_final` = 0,58638 (ψ); `rel_l2_rho_final` = 0,69322 (\|ψ\|²) | — |
| Parametreler | `results/pinn/W56_fourier_uniform_s0/config.json` | — | — |
| Animasyon | kök dizindeki `wavefunction_animation.html` (commit `8f013cf`) | gömülü 31 kare × {PINN, DST} \|ψ\|², 59 × 59 ızgara | yeniden hesaplanmamış, yalnızca yeniden biçimlendirilmiştir |

**Mevcut 4 panelli şeklin kaynağı.** `figures/pinn_vs_reference_W56.png`, `make_figures.py` (satır 432–456) tarafından
aynı `eval_fields.npz` dosyasından üretilmiştir:

- panel 1: `V`;
- panel 2 ve 3: `|psi_ref[-1]|²` ve `|psi_pinn[-1]|²`;
- panel 4: **işaretli** fark PINN − referans.

Yeni Şekil 3 ise istendiği gibi **mutlak** farkı gösterir.

**Kaynak dosyaların SHA-256 önekleri** (görseller üretilirken):

| Dosya | SHA-256 öneki |
|---|---|
| `eval_fields.npz` | `61f7fcbb5b40f9ec` |
| `eval_summary.json` | `f29375949d913e8f` |
| `config.json` | `0d24c05954b4e491` |
| `wavefunction_animation.html` | `94ba62bdc0c73da0` |

### 2.2 Parametreler, ızgara ve normalizasyon

**Problem.** 1 × 1 kapalı kutu, dört duvar sert (ψ = 0). Bu, ana taşınım deneylerinin dalga kılavuzu değil, **PINN test
problemidir**.

**Fiziksel parametreler** (`config.json`):

| Parametre | Değer |
|---|---|
| birimler | ħ = m = 1; uzunluklar boyutsuz |
| düzensizlik gücü W (`disorder_strength`) | 56,25 (= E₀/2) |
| korelasyon uzunluğu ζ (`disorder_corr_length`) | 0,08 |
| düzensizlik tohumu | 0 |
| mod sayısı M | 512 |
| başlangıç paketi | Gauss; merkez (0,35; 0,5); σ = 0,07; k_x0 = 15 |
| son zaman | T = 0,02 |

**PINN ağı.** 4 × 64, Fourier öznitelikli; ağ tohumu 42.

**Izgara.** 127 × 127 iç nokta; x, y ∈ [1/128, 127/128], aralık 1/128.

**Zaman.** `evaluate.py` 21 zaman noktasının her 5.'sini saklar: t = 0; 0,005; 0,01; 0,015; 0,02. Şekil 2–3, son
kareyi (t = T = 0,02) kullanır.

**Normalizasyon:**

- ∑\|ψ_ref\|² h² = 1,000.
- ∑\|ψ_PINN\|² h² = 0,988 (t = T'de). PINN'in normu tam korunmaz; `norm_max_dev` = 0,039.
- Şekillerde \|ψ\|² değerleri **dönüştürülmeden** gösterilmiştir. Yeniden normalize etme yapılmamıştır.

### 2.3 Aynı realizasyon ve aynı zaman noktası kontrolleri

Bu kontroller yalnızca okuma ve deterministik fonksiyon değerlendirmesidir; zaman ilerletme yoktur.

1. **V aynı realizasyondur.** `potential.potential_numpy(X, Y, cfg)`, `config.json` ile yeniden değerlendirildi.
   Saklanan `V` ile en büyük fark **0,0**'dır. Bu, deterministik potansiyel formülünün değerlendirilmesidir;
   simülasyon değildir.
2. **ψ_ref ve ψ_pinn aynı zaman noktasındadır.** Saklanan alanlardan hesaplanan bağıl L2 hataları (0,586381… ve
   0,693223…), `eval_summary.json`'daki `rel_l2_final` ve `rel_l2_rho_final` ile 10⁻¹² içinde aynıdır. Betik bunu
   `assert` ile denetler.
3. **Animasyon aynı koşuya aittir.** İki ızgaranın ortak 9 noktasında (x, y ∈ {0,25; 0,5; 0,75}) karşılaştırma
   yapıldı:
   - Animasyondaki PINN \|ψ\|² değerleri, `eval_fields.npz`'deki değerlerle **bit düzeyinde aynıdır** (t = 0; 0,01;
     0,02 kareleri).
   - Animasyondaki DST değerleri ~10⁻⁵ bağıl düzeyde uyuşur. Fark, DST çözümünün animasyon için 59 × 59 ızgarada
     ayrıca hesaplanmış olmasından gelir.
   - Kare zamanları t_i = i·T/30'dur (i = 0…30). Bu, özgün kare etiketlerinin 31'iyle birebir uyuşur ve üç karede veriyle
     doğrulanmıştır.
4. **Eksik veri yoktur.** Ekran görüntüsünden veri çıkarılmamıştır.

**Sınırlama (animasyon kökeni).** Animasyonun 59 × 59 ızgara ve 31 kareyle hangi komutla üretildiği depoda
kayıtlı değildir. `visualize.py` 89 × 89 ızgara ve 41 kare kullanır. Animasyonun içeriği yukarıdaki karşılaştırmalarla
W56 koşusuna bağlanmıştır.

---

## 3. Türkçe şekil açıklamaları (rapor için önerilen metin)

**Şekil 1. İki boyutlu düzensizlik potansiyelinin 3B yüzey gösterimi.**

> PINN test probleminde kullanılan Gauss-korelasyonlu rastgele potansiyel V(x, y) (W = 56,25, ζ = 0,08, düzensizlik
> tohumu 0; 1 × 1 kapalı kutu, 127 × 127 ızgara, boyutsuz birimler). Fiziksel model iki boyutludur; düşey eksen
> yalnızca V'nin değerini gösterir, üçüncü bir uzay boyutu değildir. Gösterilen, tek bir düzensizlik
> realizasyonudur; ansambl ortalaması değildir. Bu potansiyel, ana taşınım deneylerinin dalga kılavuzu geometrisine ait
> değildir. Veri: `results/pinn/W56_fourier_uniform_s0/eval_fields.npz` (`V`).

**Şekil 2. PINN ve DST referans dalga fonksiyonu karşılaştırması.**

> Aynı düzensizlik realizasyonunda (W = 56,25, ζ = 0,08, tohum 0) t = T = 0,02 anındaki olasılık yoğunluğu
> \|ψ(x, y)\|². Sol panel referans çözümdür (DST tabanlı split-step çözücü); sağ panel PINN tahminidir. İki panelde
> ızgara, z aralığı, renk ölçeği ve kamera açısı aynıdır. Farklar gizlenmemiştir: PINN'in tepe değeri 9,31, referansın
> tepe değeri 6,66'dır ve yoğunluğun uzaysal dağılımı belirgin biçimde farklıdır. Son anda bağıl L2 hatası ψ için
> 0,586, \|ψ\|² için 0,693'tür. Fiziksel model iki boyutludur; 3B yüzey yalnızca bir gösterim biçimidir. Bu tek bir
> realizasyon ve tek bir ağ tohumudur; ansambl ortalaması veya temsili bir PINN performansı değildir. PINN bu çalışmada
> ikincil bir yöntemsel karşılaştırmadır. Veri: `results/pinn/W56_fourier_uniform_s0/eval_fields.npz`,
> `eval_summary.json`.

**Şekil 3. PINN tahmininin DST referansından mutlak farkı.**

> Şekil 2'deki iki yoğunluk arasındaki **mutlak** fark \| \|ψ_PINN\|² − \|ψ_DST\|² \|, t = T = 0,02 (işaretli fark
> değildir; işaret bilgisi bu gösterimde yoktur). DST çözümü referans alınmıştır. En büyük mutlak fark 4,02'dir.
> Kaynak veride saklanan bağıl L2 hatası, son anda ψ için 0,586, ρ = \|ψ\|² için 0,693'tür. Fiziksel model iki
> boyutludur; 3B yüzey yalnızca bir gösterimdir. Tek bir düzensizlik realizasyonu (W = 56,25, ζ = 0,08, tohum 0) ve tek
> bir ağ tohumu gösterilmektedir; ansambl ortalaması değildir. Veri:
> `results/pinn/W56_fourier_uniform_s0/eval_fields.npz`, `eval_summary.json`.

**Animasyon (sergi / ek belge).**

> Aynı realizasyonda (W = 56,25, ζ = 0,08, tohum 0) dalga paketinin olasılık yoğunluğunun t = 0'dan T = 0,02'ye zamanla
> değişimi: solda DST referans çözümü, sağda PINN tahmini (31 kare, 59 × 59 ızgara). Fiziksel model iki boyutludur; 3B
> yüzey bir gösterim biçimidir. Farklar gizlenmemiştir. Tek bir realizasyondur.

**Sayısal ifade notu.** Kullanıcı görevindeki "≈ 0,59" değeri **ψ'nin** bağıl L2 hatasıdır (0,586). Şekil 2–3
**yoğunluğu** (\|ψ\|²) gösterdiği için yoğunluk hatası (0,693) da yazılmıştır. İkisi karıştırılmamalıdır.

---

## 4. Animasyon (Phase 4)

**Özgün animasyon korunmuştur:**

- `wavefunction_animation.html` (kök dizin; `index.html` buna yönlendirir);
- koyu tema, tek sahne;
- PINN renkli ve DST gri-saydam üst üste çizilir.

**Yeni sürüm:** `visualizations/tubitak/wavefunction_animation_white.html`.

- Veri, özgün dosyada gömülü kareler okunarak alınmıştır. Yeniden simülasyon veya PINN değerlendirmesi yoktur.
- Görsel değişiklikler:
  - beyaz arka plan;
  - yan yana iki sahne, panel başlıklarıyla açıkça adlandırılmış ("DST referans çözümü" / "PINN tahmini");
  - aynı z aralığı ve renk ölçeği, tek ortak renk çubuğu;
  - kare başına güncellenen kutulu zaman göstergesi (t = 0,0000 … 0,0200);
  - alt kısımda parametre bilgisi (W, ζ, tohum, k_x0, T, "tek realizasyon", "2B model, 3B yüzey gösterimi") ve son
    andaki hata değerleri;
  - Türkçe etiketler ve Türkçe ondalık virgül.
- **Oynat/Duraklat ve zaman kaydırıcısı korunmuştur.** Ek bir düzeltme yapıldı: özgün kare adları 3 ondalığa
  yuvarlandığı için çakışıyordu (ör. iki kare "t=0.001"). Yeni sürümde her karenin adı benzersizdir ve kaydırıcı
  etiketleri 4 ondalıklıdır.
- **z ölçeği düğmesi:**
  - **"tam (0–68)":** varsayılan; t = 0'daki dar paketin tepe değerini (~65) de içerir.
  - **"yakın (0–20,7)":** t ≥ T/2 karelerinin en büyük değerine göre seçilmiştir; geç anlardaki farkları görünür kılar.
    Bu ölçekte **erken kareler (t < T/2) z ekseninde kesilir**. Veri değişmez, yalnızca görüntü aralığı değişir.

---

## 5. Rapor ve sergi entegrasyonu (Phase 6)

**Önemli bağlam.** Bu üç şekil, **PINN test problemine** (1 × 1 kapalı kutu, W = 56,25, ζ = 0,08, t ≤ 0,02) aittir.
Raporun **ana sorusunu** (ζ'nin dalga kılavuzundaki iletime etkisi) yanıtlayan veriler değildir. Bu nedenle ana rapordaki
yerleri sınırlı olmalıdır (bkz. `TUBITAK_REPORT_MAPPING_R0.md` §3, §5: PINN ikincil sonuçtur ve sayfa bütçesi dardır).

| Görsel | Önerilen yer | Gerekçe |
|---|---|---|
| **Şekil 2** (PINN ve DST yan yana) | **Ana rapor, Bulgular, PINN alt bölümü** (en fazla bir PINN alan görseli). Sayfa yetmezse EK BELGE 3 | Aynı ölçekte yan yana gösterim, PINN'in DST'den farkını en dürüst ve en anlaşılır biçimde gösterir. Ana rapordaki `figures/pinn_errors.png` ile birlikte kullanılabilir; ikisinden birini seçmek gerekirse, hatanın bütün yapılandırmalar için nicel özeti olan `pinn_errors.png` ana metinde kalmalı, Şekil 2 EK BELGE 3'e gitmelidir |
| **Şekil 3** (mutlak fark) | **EK BELGE 3** (PINN ayrıntıları) | Şekil 2'nin tamamlayıcısıdır; ana soruya katkısı yoktur |
| **Şekil 1** (V(x, y)) | **EK BELGE 3**, Şekil 2–3'ün bağlamı olarak | Ana raporda düzensizliği göstermek için dalga kılavuzu geometrisindeki ve üç ζ değerini karşılaştıran `figures/disorder_examples.png` daha uygundur; Şekil 1 PINN kutusuna aittir |
| **Animasyon (beyaz)** | **Sergi standı** (dokunmatik ekran / dizüstü, çevrimdışı açılır) | Etkileşimli; zaman kaydırıcısı ve yan yana karşılaştırma ziyaretçiye PINN–DST farkını doğrudan gösterir |
| **Şekil 2 (PNG)** | **Sergi posteri** | Statik baskı için en açıklayıcı görsel; yanına "2B model, 3B gösterim; tek realizasyon; DST referans" notu konmalıdır |

**Raporda kullanırken:**

- Şeklin içindeki başlık, Word'deki "Şekil n." açıklamasıyla aynı bilgiyi tekrarlar. Raporda istenirse başlıksız sürüm
  üretilebilir: `base_layout` içindeki başlık satırı kaldırılır.
- Genişlik önerileri:
  - Şekil 2: tam sayfa genişliği (16 cm). Bu genişlikte eksen yazıları ≈ 7–9 pt olur.
  - Şekil 1 ve 3: ≈ 11–12 cm.

---

## 6. Sınırlamalar

1. **2B model.** 3B yüzeyler 2B alanların gösterimidir. Üç boyutlu fiziksel bir simülasyon yapılmamıştır ve bu
   görseller öyle sunulmamalıdır.
2. **Tek realizasyon, tek ağ tohumu.** Hiçbir görsel ansambl ortalaması veya temsili PINN performansı değildir.
   - PINN hatasının yapılandırmalar arası özeti `results/pinn/summary.csv` ve `figures/pinn_errors.png`'dedir.
   - Temiz ile W = 56 karşılaştırması (grup başına 3 tohum; Welch p = 0,086) istatistiksel olarak anlamlı değildir.
3. **PINN test problemi ≠ ana deney.** Görseller, ana sonuç olan ⟨T⟩(ζ) taramasının dalga kılavuzu geometrisine ait
   değildir.
4. **PINN normu.** PINN normu tam korunmaz (t = T'de 0,988). Yoğunluk farkının bir kısmı bundan gelir. Veriler
   yeniden normalize edilmemiştir.
5. **Animasyon ızgarası.** Animasyonun DST verisi 59 × 59 ızgarada ayrıca hesaplanmıştır ve üretim komutu kayıtlı
   değildir (§2.3). 127 × 127 verisiyle ~10⁻⁵ bağıl düzeyde uyuşur.
6. **SVG.** 3B sahne SVG'de rasterdır (§1).
7. **Yazı tipi.** PNG'ler Liberation Serif ile üretilmiştir; Times New Roman ile metrik olarak eşdeğerdir, glif
   biçimleri ise çok az farklıdır.
8. **Görsel seçimler.** Renk ölçekleri, kamera açısı ve ışıklandırma yalnızca görsel tercihlerdir.
   - Şekil 1: simetrik `RdBu_r`, 0'da beyaz.
   - Şekil 2: `Viridis`.
   - Şekil 3: `YlOrRd`.
   - Veriye enterpolasyon veya yumuşatma uygulanmamıştır. Plotly yüzeyleri ızgara noktaları arasında doğrusal
     üçgenleme ile çizer.
9. **DST referansının kapsamı.** DST çözümü bu kapalı-kutu probleminde referanstır (ayrıklaştırma hatası ≈ 10⁻⁶; R0.3
   Ek A). Bu durum, düzensiz taşınım sonuçlarının bağımsız durağan doğrulaması anlamına gelmez; B4 yapılmamıştır.
