---
title: "Quantum-Schrodinger-PINN"
subtitle: "Düzensiz dalga kılavuzlarında Recursive Green's Function (RGF) çözücüsünün yuvarlama hatası ve W/Z formülasyonu: bulgular, yazılım doğrulaması ve bağımsız R0 incelemesi talebi"
lang: tr
---

::: {.kapak}
::: {.kapak-ust}
Bağımsız Teknik İnceleme İçin Proje Özeti
:::

# Quantum-Schrodinger-PINN {.unnumbered .unlisted}

::: {.kapak-alt}
Düzensiz dalga kılavuzlarında Recursive Green's Function (RGF) çözücüsünün yuvarlama hatası ve W/Z formülasyonu:
bulgular, yazılım doğrulaması ve bağımsız R0 incelemesi talebi
:::

| | |
|---|---|
| **Hazırlayan** | [Ad Soyad] |
| **Kurum / danışman** | [doldurulacak] |
| **Tarih** | 14 Kasım 2026 |
| **Belge türü** | Bağımsız Teknik İnceleme İçin Proje Özeti |
| **Depo** | `github.com/mrthnkzlrmk865-cmd/Quantum-Schrodinger-PINN` |
| **Dal / incelenen commit** | `claude/quantum-schrodinger-audit-iy6kn4` / `1668540f2fd25d9903ad005d31b06f4c340798f6` |
| **Dondurulmuş şartname** | `PREDICTIONS_coupling_scaling_v2.md`, Revizyon 3.7 (`0f502c9`) |

::: {.uyari}
**Bu belge bir başarı ilanı değildir.**

- B1-F: **FAIL**.
- W/Z: **NOT VALIDATED**.
- Ana coupling-scaling deneyi: **NOT RUN**.
- Uygulama düzeyi R0: **PENDING**.
- Bu belgede hiçbir hakem onayı bulunmamaktadır.
:::
:::

# Yönetici özeti

Proje, iki boyutlu düzensiz bir dalga kılavuzunda zamandan bağımsız Schrödinger denkleminin saçılma (S) matrisini
kanal uzayında (*mode space*) özyinelemeli Green fonksiyonu (*Recursive Green's Function*, RGF) yöntemiyle hesaplıyor.
Ardından bu hesabın kayan nokta (binary64) yuvarlama hatalarını önceden kaydedilmiş (*preregistered*) protokollerle
inceliyor. Temel bulgular:

1. **B1 uzun cihaz testi.** Altı temiz kılavuz konfigürasyonundan beşi geçti. En uzun cihazda (F: N = 3483 dilim) faz
   büyüklüğü 1.76851694417934e-10 çıktı ve önceden belirlenmiş 1e-10 kriterini aştı. 100 basamaklı kafes referansına göre
   gerçek RGF faz hatası 1.5424928001550615e-10'dur. Kriter gevşetilmedi; **B1-F = FAIL**.
2. **W/Z formülasyonu.** Aynı ayrık problemin küçük-büyüklüklü değişkenlerle yazılmış cebirsel bir eşdeğeridir.
   - Tek bir temiz noktada (B1-F'nin en kötü satırı) faz hatasını ~1e-15 düzeyine indirdi.
   - Buna karşılık genlik hatasını büyüttü. Tüm açık kanallarda maksimum ~1.9e-13 ölçüldü; bu, orijinal yöntemin 3.4
     katıdır.
   - Önceden kaydedilmiş kurala göre sonuç **INCONCLUSIVE** oldu.
3. **Birbiriyle değişmeyen (non-commuting) bloklarla tasarım deneyi.** 76 execution / 72 benzersiz konfigürasyon
   çalıştırıldı.
   - Referanslar ve W/Z hataları ε düzeyinde kaldı.
   - Önceden kaydedilmiş bir ölçekleme kuralı ise tetiklendi; dört eğim 0.78–0.79 çıktı ve alt sınır 0.8'in altında kaldı.
   - Sonuç **DESIGN INSUFFICIENT** oldu. Bu, W/Z'nin başarısız olduğu anlamına gelmez.
4. **Coupling-scaling v2 deneyi.** Altyapısı yazıldı ve XP/HO manifestleri üretilip donduruldu (`fac5f5e`). T1–T6 resmi
   yazılım testleri geçiyor. Ana deney **çalıştırılmadı**. Uygulama düzeyi bağımsız R0 incelemesi henüz yapılmadı ve
   şartnamedeki bir tanım boşluğu (G-1) çözülmeden driver çalışmayı reddediyor.

| Konu | Durum |
|---|---|
| B1-F | **FAIL** |
| Coupled-block W/Z | **DESIGN INSUFFICIENT** |
| W/Z genel doğrulaması | **NOT VALIDATED** |
| Main coupling-scaling experiment | **NOT RUN** |
| Implementation-level R0 | **PENDING** |
| XP/HO manifests | **GENERATED AND FROZEN** |
| T1–T6 software validation | **PASS** (yalnızca yazılım ve protokol testi) |

: Tablo 1. Güncel bilimsel durum. Durumlar birbirinin yerine kullanılamaz: yazılım testlerinin geçmesi, bilimsel
doğrulama anlamına gelmez.

# Bilimsel problem ve motivasyon

Düzensiz (*disordered*) ortamlarda dalga taşınımı, Anderson lokalizasyonu [1] ve mezoskopik iletkenlik [2–4] gibi
olguların temelidir. Projenin taşınım hesapları, iki uçlu (*two-probe*) bir dalga kılavuzunun S matrisine dayanır. Bu
matris, cihazın Green fonksiyonundan Fisher–Lee bağıntısıyla [3] elde edilir.

Uzun cihazlarda ve eşik enerjilerine (*threshold*) yakın kanallarda, dilim başına küçük yuvarlama hataları birikebilir.
Bu birikim, önceden belirlenmiş doğruluk kriterlerini aşabilir.

Projenin incelediği soru şudur: **Mevcut RGF çözücüsünün ölçülen faz hatası nereden kaynaklanıyor? Cebirsel olarak eşdeğer
bir yeniden yazım olan W/Z formülasyonu bu hatayı gerçekten azaltıyor mu ve başka bir hata türünü büyütüyor mu?**
Bu soru doğrulama protokolleri, önceden kaydedilmiş karar kuralları ve bağımsız inceleme kapılarıyla (R0) ele alınıyor.

Faz ve genlik hataları **ayrı ayrı** değerlendirilir, çünkü farklı mekanizmalardan kaynaklanırlar:

- **Faz hatası.** Eşik yakınında açık kanalın dalga sayısı k küçüktür. Dilim başına bir köşegen yuvarlaması, faz
  hatasını yaklaşık ε/sin(k dx) ile büyütür.
- **Genlik hatası.** Böyle bir büyütme içermez; ancak dilim başına O(ulp(1)) birikerek N ile doğrusal büyüyebilir
  (bölüm 6).

Bir yöntemin fazı iyileştirmesi, genliği iyileştirdiği anlamına gelmez.

# Matematiksel model: RGF ve W/Z

Aşağıdaki tanımlar `stationary.py`, `wz_rgf.py`, `PREDICTIONS_coupled_blocks.md` §1 ve dondurulmuş şartnamenin §N–§2
bölümlerinden alınmıştır.

**Ayrıklaştırılmış model.** Dalga kılavuzu x yönünde dx aralıklı N dilime bölünür. Her dilimde enine $N_{y}$ kanal
(enine özfonksiyon) bulunur ve enine enerjiler $E_{\perp,n}$'dir. Dilim içi ve dilimler arası bloklar şöyledir:

$$t = \frac{1}{2\,dx^2},\qquad H_{jj} = 2t\,I + \mathrm{diag}(E_\perp) + V_j,\qquad H_{j,j\pm1} = -t\,I .$$

Burada $V_{j}$, j. dilimin potansiyelinin kanal uzayındaki matrisidir. Sistem matrisi blok üç köşegenlidir
(*block-tridiagonal*):

$$M = E - H_D - \Sigma_L \oplus \Sigma_R,\qquad G = M^{-1}.$$

**Açık sınır koşulları ve lead öz-enerjisi.** Sol ve sağ uçlara yarı sonsuz temiz kılavuzlar (*lead*) bağlıdır. Etkileri
yalnızca 1. ve N. dilime eklenen öz-enerji (*self-energy*) terimleriyle hesaba katılır:

$$\Sigma = \mathrm{diag}(-t\lambda_n) = -t\Lambda,\qquad \lambda_n + \lambda_n^{-1} = 2(1-\delta_n),\qquad \delta_n=(E-E_{\perp,n})\,dx^2 .$$

Burada $λ_{n}$ geciktirilmiş (*retarded*) köktür. Açık kanalda |$λ_{n}$| = 1 olur, kapalı kanalda |$λ_{n}$| < 1 olur.

**Orijinal RGF** (soldan sağa; `stationary.rgf_endpoint_blocks`). $A_{j}$ = E − $H_{jj}$ olmak üzere:

$$D_1 = A_1 - \Sigma_L,\qquad D_j = A_j - t^2 g_{j-1}\ (j>1),\qquad g_j = D_j^{-1}$$

Son dilimde ayrıca $Σ_{R}$ çıkarılır. Uç bloklar ($G_{11}$, $G_{N1}$, $G_{1N}$, $G_{NN}$) zincir güncellemeleriyle elde edilir.

**W/Z değişkenleri** (`wz_rgf.py`):

$$\Delta_j=(E\,I-\mathrm{diag}(E_\perp)-V_j)\,dx^2,\qquad A_j = t\,(2\Delta_j-2I)\ \ (2t\,dx^2=1\ \text{olduğundan tam}),$$

$$W_j = D_j/t + I\ \ (1\le j<N),\qquad Z_N = D_N/t,\qquad \hat G = t\,G .$$

**Cebirsel ilişki.** $D_1 = A_1 - \Sigma_L$ ve $\Sigma = -t\Lambda$'dan $W_1 = 2\Delta_1 + (\Lambda - I)$ çıkar. $j > 1$ için $t\,g_{j-1} = (W_{j-1} - I)^{-1}$
eşitliği ve (I − W)⁻¹W = −I − (W − I)⁻¹ özdeşliği kullanılarak şu elde edilir:

$$W_j = 2\Delta_j + (I-W_{j-1})^{-1}W_{j-1},\qquad Z_N = 2\Delta_N + (\Lambda - I) + (I-W_{N-1})^{-1}W_{N-1}.$$

Ĝ blokları $P_{j}$ = (I − $W_{j}$)⁻¹ ile güncellenir; $P_{j}$'nin sol ve sağ zincirlerdeki çarpım sırası önemlidir. Fisher–Lee
bağıntısı açık kanallarda $t_{RL}$ = 2i√s $\hat G_{N1}$ √s şeklindedir ($s_n = \sin\theta_n = \Gamma_n/(2t)$).

Bu yeniden yazımın sayısal amacı şudur:

- Orijinal yazımda köşegen terim $A_{jj}$ ≈ −2t ≈ −10⁶ (dx = 1/1024) mertebesindedir. Eşik yakınında asıl fiziksel bilgi bu
  büyük sayının küçük bir farkındadır (E − $E_{\perp}$).
- W/Z, bu farkı doğrudan $Δ_{j}$ = (E − $E_{\perp}$)dx² büyüklüğünde taşır.

::: {.not}
**Önemli uyarılar.**

1. W/Z ile orijinal RGF, **tam aritmetikte** aynı ayrık problemi çözer. Bu cebirsel eşdeğerlik, binary64'te sayısal
   üstünlük ya da fiziksel doğrulama anlamına gelmez.
2. Her çözücü, kendi binary64 girdilerinin tanımladığı **temsil edilen kafesi** (*represented lattice*) çözer:
   - L: tam kafes;
   - $R_{W}$: W/Z'nin gördüğü kafes ($Δ_{\mathrm{fl}}$);
   - $R_{\mathrm{orig}}$: üretim kodunun gördüğü kafes.

   Hata ayrıştırmaları bu kafeslerin farkını açıkça ayırır (bölüm 5–6).
3. Kontrollerin anlamları:
   - U = ‖S†S − I‖: akı korunumu (unitarity).
   - Q = ‖S − Sᵀ‖: karşılıklılık (*reciprocity*). Gerçek simetrik (R sınıfı) potansiyelde Q = 0 beklenir. Karmaşık
     Hermitsel (H sınıfı) potansiyelde S simetrik olmak zorunda değildir; Q orada referansla karşılaştırılır.
   - Hermitsellik (*Hermiticity*): potansiyel matrisi için yapısal bir koşuldur.

   §9.4'teki bulgu, S tabanlı testlerin bazı sıralama hatalarını **göremeyebileceğini** gösterir (bölüm 7).
:::

# B1 uzun cihaz bulguları

**Kaynak.** `docs/RGF_B1_VALIDATION.md` (2026-11-06, temel commit `bce6587`); yüksek hassasiyetli faz denetimi
`audit_rgf_phase.py` (kaynak commit `d5a97a3`); sonuçlar `results/rgf_b1_validation/summary.json` ve
`phase_audit_hp_summary.json`.

B1, potansiyelsiz (V = 0) temiz kılavuzdur. Kafes üzerinde $t_{nn}$ = $λ_{n}$^{N−1} analitik cevabı vardır. Bu nedenle B1 bir
**iç tutarlılık testidir, bağımsız fiziksel doğrulama değildir**. Her konfigürasyonda 50 enerji noktası ve her iki yön
kullanılmıştır.

| Konfig. | $N_{y}$ | dx | $L_{d}$ | N | (i) B1 faz büyüklüğü | (ii) RGF vs 100 basamak | B1 faz kriteri |
|---|---|---|---|---|---|---|---|
| A | 63 | 1/512 | 1 | 718 | 1.878e-11 | 2.360e-11 | geçti |
| B | 127 | 1/512 | 1 | 718 | 1.878e-11 | 2.360e-11 | geçti |
| C | 63 | 1/1024 | 1 | 1435 | 7.289e-11 | 6.358e-11 | geçti |
| D | 127 | 1/1024 | 1 | 1435 | 7.289e-11 | 6.358e-11 | geçti |
| E | 63 | 1/512 | 3 | 1742 | 4.560e-11 | 5.765e-11 | geçti |
| **F** | **63** | **1/1024** | **3** | **3483** | **1.76851694417934e-10** | **1.5424928001550615e-10** | **FAIL** |

: Tablo 2. B1 faz sonuçları. Değerler 50 enerji üzerinden maksimumdur; (ii) 424 satır üzerinden alınmıştır. Kriter 1e-10
olarak önceden belirlenmiştir ve **değiştirilmemiştir**.

![Şekil 1. B1 konfigürasyonlarında maksimum faz hatası. Veri: `results/rgf_b1_validation/phase_audit_hp_summary.json`. Ölçülmemiş veya tahmini değer içermez.](turkce_ozet_figurler/sekil1_b1_faz.png)

**Hata kaynağı.** En kötü nokta F'dedir: E = 125.692…, kanal n = 5, sin θ ≈ 2.1e-3. Bu noktada faz hatası 110 basamakta
tam bir teleskopik ayrıştırmayla bileşenlerine ayrıldı (`experiments/wz_rgf_diagnostic/REPORT.md` §4; kalıntı 0):

| bileşen | katkı |
|---|---|
| saklanan köşegenin temsili, $R_{\mathrm{dev}}$ − L | +1.6108e-10 (≈ %104) |
| özyineleme yuvarlaması | −6.88e-12 (≈ −%4) |
| lead uyumsuzluğu | +5.3e-14 |
| **toplam** | **+1.5425e-10** |

Yani bu noktada hatanın neredeyse tamamı, **saklanan köşegen $A_{j}$ değerinin binary64 yuvarlamasından** gelmektedir.

Ek gözlemler:

- B1 raporu §6 tüm konfigürasyonlarda bir **üst sınır** ilişkisi ölçmüştür: |Δφ| ≤ 0.89·N·ε/sin(k dx). Bu oranın medyanı
  0.21–0.23'tür.
- Sabit dx'te hata N ile doğrusal büyür.
- dx yarıya indirildiğinde ise sabit bir çarpanla ölçeklenmez.

::: {.not}
**Doğrulanamayan iddia.** Görev tanımında şu ifadeler yer alıyor:

- tahmin edilen yuvarlama hatası ile gözlenen işaretli faz hatası arasında "korelasyon ≈ 0.91–0.93";
- "eğim ≈ 0.99–1.05".

Bu değerler depodaki **hiçbir commit edilmiş rapor, sonuç dosyası veya commit mesajında bulunamamıştır**, bu yüzden bu
belgede kanıt olarak kullanılmamıştır. Depodaki kanıtlar, yukarıdaki tek noktalı tam ayrıştırma ile tüm satırlardaki sınır
ilişkisidir. Kalan hatanın tamamı, tüm satırlar için tek bir nedene bağlanmış değildir. Köşegen yuvarlaması
**kanıtlanmış tek neden olarak sunulmamaktadır**.
:::

# Scalar ve 63-kanal genlik yuvarlaması deneyleri

## İlk scalar W/Z deneyi (tek nokta, tek kanal)

**Kaynak.** `experiments/wz_rgf_diagnostic/REPORT.md`; referans commit `7422718`; tahminler `PREDICTIONS.md`
(`3ce732d`); düzeltme notu `ERRATUM_worst_point.md`.

Deney B1-F'nin en kötü satırında yapıldı: E = 125.69231487269683, kanal 5.

| ölçüt | orijinal (üretim) | W/Z (δ tabanlı lead) |
|---|---|---|
| faz $ε_{φ}$ (L'ye göre) | +1.5425e-10 | −1.19e-15 |
| genlik $ε_{A}$ | −1.08e-15 | +4.82e-14 (44 kat daha kötü) |

Önceden kaydedilmiş karar kuralı (faz en az 10 kat iyileşmeli **ve** genlik 10 kattan fazla kötüleşmemeli) uygulandığında
sonuç **INCONCLUSIVE** çıktı.

Erratum, ilk raporu üç noktada düzeltir:

1. W/Z yalnızca **1×1 skaler yolda** (kanal 5) çalıştırıldı.
2. Kapalı kanal dalı **hiç çalışmadı**.
3. U ve Q, 10×10'luk tam açık kanal S'si değil, **2×2'lik kanal-5** büyüklükleridir.

Bu deney genel bir W/Z üstünlüğü göstermez.

## Tam 63-kanal genlik yuvarlaması deneyi

**Kaynak.** `amplitude_roundoff_report.md`; temel commit `74227189`; yayımlanan çıktılar 3. koşudan (`e274064`);
`results_amplitude_roundoff/amplitude_roundoff_results.json`.

- **Kurulum.** $N_{y}$ = 63, yani 5 açık ve 58 kapalı kanal; N = 3483; tam açık S matrisi 10×10. Temiz cihaz olduğu için tüm
  bloklar köşegendir.
- **Köşegen dışı terimler.** $W_{j}$, $P_{j}$, Ĝ blokları ve $Z_{N}$'nin köşegen dışı elemanları **tam olarak 0** çıktı; S
  bloklarında kanal karışımı da 0'dır.
- **Faz.** W/Z'nin faz hatası tüm açık kanallarda en fazla **3.24e-15**'tir (kendi kafesine göre). Orijinal yöntemin L'ye
  göre faz hatası 1.0e-11 ile **1.54e-10** arasındadır.
- **Genlik.** W/Z'nin maksimum |$ε_{A}$| değeri **1.93e-13**, üretimin 5.7e-14'tür; oran ≈ **3.4**. Kanal 5'teki "44 kat"
  oranı kanala özgüdür.
- **Mekanizma (kanal 5).**
  - Tüm çözme yöntemleri, doğru yuvarlanmış aynı P = fl($λ_{W}$) değerini döndürür.
  - Bu değerin modül kusuru |P| − 1 = **+1.37e-17**'dir.
  - 3482 adım boyunca birikince **4.757e-14** olur; bu, $ε_{A}$ = 4.82e-14'ün %98.7'sidir.
- **Yineli iyileştirme** (*iterative refinement*). Residual'ı yarıya indirdi (2.2e-16 → 1.1e-16), ancak genlik hatasını
  **azaltmadı**: max |$ε_{A}$| 1.93e-13'ten 2.34e-13'e çıktı.
- **Diğer ölçütler.** U ≤ 5.24e-13, Q = 0. Kapalı kanallarda $\hat G_{N1}$ göreli hatası ≤ 4.8e-13'tür.

![Şekil 2. 63 kanallı tek noktada açık kanal başına $|ε_{φ}|$ ve $|ε_{A}|$. W/Z kendi kafesi $R_{W}$'ye göre, orijinal RGF L'ye göre ölçülmüştür. Bu noktada $R_{W}$ − L faz farkı ≤ 1.7e-15'tir. Veri: `amplitude_roundoff_results.json` (`runs["63x63:solve"]`, `production_original`).](turkce_ozet_figurler/sekil2_63kanal.png)

**Sınırlama.** Bu deney tek enerji ve tek temiz cihazla yapılmıştır. Şunları **test etmez**:

- birbiriyle değişmeyen (*non-commuting*) bağlı bloklar; temiz cihazda bütün çarpımlar köşegen olduğundan çarpım sırası
  hiç sınanmamıştır;
- düzensizlik (*disorder*);
- N ile ölçeklenme.

# Coupled-block (birbiriyle değişmeyen blok) deneyi

**Kaynak.** `PREDICTIONS_coupled_blocks.md` (tasarım `7d69260`); `results_coupled_blocks/coupled_blocks_report.md`
(koşu `63e046d`, 2026-11-10).

**Kurulum.** Küçük, deterministik, bilerek bağlı ($V_{j}$ köşegen dışı) sistemler kullanıldı:

- boyutlar N ≤ 4, $N_{y}$ ≤ 3;
- sınıflar R (gerçek simetrik), H (karmaşık Hermitsel) ve D (köşegen kontrol);
- **76 execution, 72 benzersiz konfigürasyon**; 4 tekrarlanan execution bit düzeyinde aynı çıktı.

| kontrol | en kötü değer | eşik |
|---|---|---|
| (a) `mp.lu_solve` vs `mp.inverse` | 0 | 1e-55 |
| (b) 100 basamak vs 60 basamak | 2.0e-61 | 1e-55 |
| (c) mp residual | 3.2e-61 | 1e-55 |
| W/Z uç bloklar ($e_{\mathrm{glob}}$) | ≤ 2.5e-16 | sezgisel $τ_{\mathrm{end}}$ |
| W/Z prefix blokları | ≤ 3.1e-16 | $τ_{\mathrm{pre}}$ |
| $W_{j}$ / $P_{j}$ | ≤ 4.8e-16 / ≤ 5.2e-16 | $τ_{W}$ / $τ_{P}$ (≤ %1'i) |
| enstrümanlı kopya = `wz_endpoint_blocks` (bit düzeyi) | 76/76 | tümü |
| mutantlar: algılanabilir vakalar algılandı | 480/480 | $D_{\mathrm{mut}}$ ≥ 1e-6 |

: Tablo 3. Coupled-block deneyinin özet ölçümleri. Görev metnindeki "residual ≤ 3e-61" ifadesinin rapordaki kesin
karşılığı 3.2e-61'dir.

![Şekil 3. Coupled-block deneyinde W/Z (solver A) uç blok hatası, sezgisel sınır $τ_{\mathrm{end}}$ ve 1e-10 kriteri ile birlikte. Veri: `coupled_blocks_comparisons.csv` (`row = "solver A"`).](turkce_ozet_figurler/sekil3_coupled.png)

**Karar: DESIGN INSUFFICIENT.**

- Sıfır geçersizlik (*invalid*), sıfır başarısızlık (*fail*) ve sıfır anomali vardır.
- Tetiklenen tek kural N8 (coupling ölçekleme eğimi) oldu: dört sıralama mutantında log-log eğim 0.78–0.79 çıktı ve
  önceden kayıtlı alt sınır 0.8'in altında kaldı.
- Raporun §9.2 bölümü sonradan (*post-hoc*) şu açıklamayı verir:
  - s = 1'de komütatörler henüz doğrusal rejimde değildir;
  - bu durum tasarımın kendi ön-uçuş (*preflight*) tablosunda zaten görülüyordu;
  - dolayısıyla bu bir tasarım hatasıdır.
- **Bu sonuç W/Z FAIL anlamına gelmez.** Önceden kayıtlı tasarımın, karar için yeterli kanıt üretemediği anlamına gelir.

Ek yapısal bulgu (§9.4): N = 2 ve tek açık kanallı R sınıfında M6/M7 mutantları $\hat G_{N1}$'i ~0.15 kadar değiştirir. Buna
rağmen S, U ve Q değişmez. Bu durumda yalnızca Green bloklarının doğrudan ters matris referansıyla karşılaştırılması
hatayı yakalar.

# Coupling-scaling v2 altyapısı ve manifestler

Bu deney (şartname Rev. 3.7) bir sorunun önceden kayıtlı cevabını arar: sıralama mutantlarının ayrımı, zayıf bağlaşımda
($V_{j}$(s) = s·$V_{j}$(1)) tam olarak belirlenen asimptotik mertebeyi double precision'da da izliyor mu? Ayrıca doğru W/Z kodu
81 noktalı ızgarada 1e-10 kriterini sağlıyor mu? Deney altyapısı şu bileşenlerden oluşur:

- **Oracle** (`wz_mutants_oracle.py`): her mutantın tam değeri, §4.2 formüllerinden ve doğrudan tersler kullanılarak
  hesaplanır.
- **Clean-room mpmath transcription** (`wz_mutants_mp.py`).
- **Sertifika ve covering modülleri.** Doğrulanmış (*validated*) ters norm sınırları kesin rasyonel aritmetikte
  hesaplanır; bunlar §2.4.1 sertifikası ve [0,1] üzerinde Neumann/Taylor yarıçapı $r_{0}$'dır.
- **§4.7 doğrulanmış referanslar** (`coupling_scaling_v2_refs.py`): ω/(1 − ω) ≤ 1e-55 kesin karşılaştırması; 60 → 120
  basamak yükseltmesi; aksi halde INVALID.
- **Manifest üreteci ve driver.** Driver yazılmış, ancak kapılıdır (*gated*) ve **hiç çalıştırılmamıştır**.

**Manifestler** (commit `fac5f5e3aae6f9ae9c151a649e9279bcd01fed95`; üretici kaynak commit'i `aa06437`). Manifest yalnızca
§8.2'nin izin verdiği çekiliş, filtre ve sertifika adımlarıyla üretildi; hiçbir çözücü veya ana referans çalıştırılmadı.

| küme | kök entropi | beklenen (R/H/D) | üretilen | redraw (enerji/koşullama/sertifika) | maks. deneme no. |
|---|---|---|---|---|---|
| XP | 202611110001 | 60 (24/24/12) | 60 (24/24/12) | 11 (11/0/0) | 2 |
| HO | 202611110002 | 120 (48/48/24) | 120 (48/48/24) | 10 (10/0/0) | 2 |

: Tablo 4. Manifest sayıları. Redraw üst sınırı 10.000 denemedir; hiçbir indekste yaklaşılmadı.

Sertifika ve doğrulama sonuçları:

- Her matris [0,1] üzerinde sertifikalandırıldı. s = 0 her yerde sertifikalıdır, dolayısıyla $r_{0}$ her konfigürasyonda
  tanımlıdır.
- Sertifikasız nokta yoktur. Sertifikalı nokta sayısı XP'de 18.690, HO'da 37.303'tür.
- En büyük binary64 ızgara koşul sayısı **382.1**'dir (sınır 1e3).
- Ayrı bir tekrar oynatma betiği şunları manifestten yeniden türetti ve **0 sorun** buldu:
  - PCG64 çekilişleri, enerji kuralı ve $V_{j}$(1);
  - sıralama ve tekillik;
  - [0,1] örtüsünün kesin rasyonel kontrolü.
- Temiz bir süreçte yapılan tam yeniden üretim, üç dosyada da **bayt düzeyinde aynı** çıktı verdi.

| dosya (`experiments/wz_rgf_diagnostic/results_coupling_scaling_v2/`) | SHA-256 |
|---|---|
| `manifest_XP.json` | `e51ffe43a1970cf8fa5cf4ddce14547e9d5a24598a8248086ed88cba12a90923` |
| `manifest_HO.json` | `9b3ecc924ab204a4ef742b42240065e70778126bd9bd9584638e17879ae034d2` |
| `manifest.sha256` | `d7be20bc64b96324ed7fb16e35e428ee164672bedb7356f1b08e79b41c65b969` |

: Tablo 5. Manifest hash'leri; bu belge hazırlanırken depodaki dosyalardan yeniden hesaplanmıştır. **Manifestlerin
dondurulmuş olması, ana deneyin çalıştırıldığı anlamına gelmez.**

# Yazılım doğrulaması

Resmi bayraksız T1–T6 koşturucusu `python3 -m cs2_t6_harness.run --out <dizin>` (dizin: `experiments/wz_rgf_diagnostic`)
ve tam pytest, commit `1668540` üzerinde çalıştırıldı. Bu belgenin hazırlığı sırasında, çözücü çalıştırmama kuralı
gereğince bu testler **yeniden çalıştırılmamıştır**; değerler o koşunun kayıtlarından alınmıştır.

| Aşama | Kontrol | Başarısız | Sonuç |
|---|---:|---:|---|
| 1 Ortam / import denetimi | 54 | 0 | PASS |
| 2 Referans öz-testleri (KAT ve öz-ispatlar) | 507 | 0 | PASS |
| 3 Ön koşullar (dondurulmuş girdiler) | 28 | 0 | PASS |
| 4 T1–T5 (§2.4.2) | 468 | 0 | PASS |
| 5 T6.1–T6.8 | 3528 | 0 | PASS |
| C4 sayısal çapraz kontrol | — | — | CONSISTENT, engelleyici değil |
| 6 Mutasyon paketi | 36 | 0 | PASS (36/36 algılandı) |

: Tablo 6. T1–T6 resmi koşu (çıkış kodu 0). Ek olarak tam pytest **691 passed, 2 xfailed**, altyapı testleri **122 passed**
sonucunu verdi. Beklenen iki başarısızlıktan (xfail) biri, B1-F faz kriterinin başarısızlığını bilerek kayıt altında tutar.

::: {.not}
**Bu sonuçların anlamı.**

- **C4** yalnızca **doğrulayıcı olmayan (non-validating) bir sayısal çapraz kontroldür**. Bağımsız doğrulama veya PASS
  kanıtı değildir.
- T1–T6, mutasyon paketi ve pytest **yazılım ve protokol kontrolleridir**. Sertifika kodunun şartnameye uyduğunu ve bilinen
  hata türlerini yakaladığını gösterirler; fiziksel modelin veya W/Z'nin doğruluğunu tek başlarına ispatlamazlar.
- Testler ile test edilen kod **aynı yapay zekâ kodlama ajanı** (Claude Code) tarafından yazılmıştır (R0 paketi, Bölüm A).
:::

# Açık sorunlar ve sınırlamalar

**G-1: "duplicate set" tanımsız.**

- Şartname §4.7 madde 4(b), dps-100 kontrolünü "k ≡ 0 (mod 4), s = 1, 2⁻⁴⁰ **ve duplicate set** üzerinde" ister. Bu küme
  şartnamenin hiçbir yerinde tanımlanmamıştır.
- Tanım **uydurulmadı**. Driver bu yüzden çalışmayı reddediyor (`DUPLICATE_SET = None`).
- Konu bağımsız R0 incelemesinde ya da açık bir şartname değişikliğiyle (*spec amendment*) çözülmelidir.

**A-1: `previous_files.sha256` zamanlaması.**

- §8.3, önceki dosyaların hash'lerinin "herhangi bir koşudan önce" kaydedilmesini ve "her aşamadan sonra" doğrulanmasını
  ister.
- Bu kayıt manifest aşamasında oluşturulmadı. İlk çözücü koşusundan hemen önce çalıştırılacak ayrı bir komut olarak
  bırakıldı.
- Manifest aşamasının da bu kapsama girip girmediği henüz bağımsız olarak karara bağlanmamıştır.

**Diğer şartname boşlukları (G-2–G-5) ve uygulayıcı seçimleri (C-1–C-16).** Bunlar keyfi biçimde çözülmemiş, R0 paketi
Bölüm I.5'te tek tek listelenmiştir.

| kod | konu | etkilediği |
|---|---|---|
| G-2 / C-1 | ızgara koşullama: binary64 SVD, tam kafesin bir kez yuvarlanmasıyla | **manifest (dondurulmuş)** |
| G-4 / C-2 | $n_{\mathrm{open}}$ = $N_{y}$ için $E_{\mathrm{upper}}$'ın kesin rasyonel hesaplanması | **manifest** |
| C-3 | filtre adımlarının sırası (yalnızca kayıt nedeni etkilenir) | **manifest** (kayıt) |
| C-4 | M10 matrisinde Λ olarak enclosure orta noktası | **manifest** |
| G-3 / C-6 | `lu_solve` toleransı 1e-55 (devralınan değer) | driver |
| G-5 | oracle yalnızca belirlenmiş bloğu verir; "dört blok için D*" | driver / çıktı şeması |
| C-5, C-7–C-16 | k ≡ 0 (mod 4) okuması; M10 için Mag; regime-model uyumu; bootstrap anahtarı; metrik çevirileri; first-deviation kapalı formları; L'de tam lead kökleri (C-15) vb. | yalnızca driver |

: Tablo 7. Açık konular. Manifesti etkileyen seçimler (C-1–C-4) dondurulmuştur. R0 bunlardan birini farklı okursa yeni bir
tasarım revizyonu gerekir.

**R0 ve bağımsızlık.**

- Uygulama düzeyi R0 **tamamlanmamıştır**. `REVIEW_coupling_scaling_v2.md` yoktur, R0 kontrol listesi boştur ve bağımsız
  hakem imzası bulunmamaktadır.
- Şartnamede ayrıca bir **tasarım düzeyi** R0 incelemesinden söz edilir (`f101224`, sonuç CONDITIONAL PASS). Bu inceleme
  yalnızca belgeye karşı (a)–(b) maddelerini kapsar ve uygulama düzeyi onayın yerini tutmaz.
- Clean-room transcription, oracle ve sertifika kodunu yazan ajanla **aynı modelden**, ayrı başlatılmış bir yapay zekâ
  alt-ajanı tarafından yazılmıştır. Alt-ajana tek bir malzeme dosyası verilmiştir. Yalıtım bir sandbox ile değil,
  talimatla sağlanmıştır.
- Bu nedenle transcription **tam bağımsız akademik doğrulama sayılamaz** (`wz_mutants_mp_PROVENANCE.md`).

**Belge içi tutarsızlık.** Dondurulmuş şartnamenin başlığında hâlâ tasarım dönemine ait şu ifade yer alıyor: "nothing in
this document has been implemented or run … no matrix of it has been generated." Bu artık doğru değildir; altyapı yazılmış
ve manifestler üretilmiştir. Şartname dondurulmuş olduğu için değiştirilmemiştir. Bu çelişkinin R0'da not edilmesi önerilir.

# R0 inceleme kapsamı

Şartname §4.3'e göre, herhangi bir çözücü koşusundan önce bağımsız bir hakem `REVIEW_coupling_scaling_v2.md` dosyasında
kimlik, tarih ve commit belirterek şu maddeleri onaylamalıdır:

- **(a)** §4.2'deki mutant değişiklik sütununun gerçek `wz_mutants.py` koduyla karşılaştırılması;
- **(b)** her oracle formülünün bağımsız matematiksel türetimi;
- **(c)** oracle uygulamasının doğruluğu (`wz_mutants_oracle.py`, ilk sapma kapalı formları dahil);
- **(d)** clean-room transcription'ın doğruluğu ve bağımsızlık sınırları;
- **(e)** sertifika, validated-residual, noktasal sertifika (§5.6) kodu ile T1–T6 altyapısı.

Bunlara ek olarak:

- **G-1** ve **A-1**'in karara bağlanması;
- C-1–C-16 ve G-2–G-5'in değerlendirilmesi;
- incelenen SHA'da resmi koşunun hakem tarafından **yeniden çalıştırılması**.

Hakemin yalnızca testleri çalıştırması yeterli değildir. Türetimleri ve kodu bağımsız olarak değerlendirmesi gerekir;
çünkü testler, oracle ve transcription tek bir yazar ekosisteminden gelmektedir. Ortak bir şartname yanlış okuması üç
uygulamada da aynı şekilde bulunabilir. §4.3'ün de belirttiği gibi, bu yalnızca madde (a) ile yakalanabilir. Gerekli
belgeler R0 paketinde (Bölüm B–I) ve kontrol listesinde (Bölüm H) hazırdır. **Bu belge onay değil, inceleme talebidir.**

# Sonuç ve akademik değerlendirme talebi

Bugüne kadar varılan sonuçlar:

- **Kaynak tespiti (tek noktada).** Üretim RGF çözücüsünün uzun cihazdaki faz hatası, tek bir noktada tam ayrıştırmayla
  saklanan köşegenin yuvarlamasına bağlanmıştır.
- **B1-F.** Önceden belirlenmiş kriterle **FAIL** olarak kalmaktadır.
- **W/Z.** Cebirsel olarak eşdeğerdir. Temiz tek noktada fazı yaklaşık 10⁵ kat iyileştirirken genliği kötüleştirmiştir.
  Birbiriyle değişmeyen bloklarda ε düzeyinde doğru davranmış, ancak önceden kayıtlı tasarım karar üretememiştir
  (DESIGN INSUFFICIENT).
- **Genel durum.** W/Z **doğrulanmamıştır**. Bunu sınayacak coupling-scaling deneyinin altyapısı ve manifestleri hazırdır,
  ancak deney R0 tamamlanmadan çalıştırılmayacaktır.

Değerlendirmesi istenen konular:

1. W/Z türetiminin ve mutant oracle formüllerinin matematiksel doğruluğu;
2. doğrulanmış sertifika yaklaşımının (kesin rasyonel ters norm sınırları, Neumann örtüsü) yeterliliği;
3. G-1, A-1 ve C-1–C-16 için karar;
4. önceden kayıtlı karar kurallarının bilimsel uygunluğu;
5. bağımsızlık sınırlarının kabul edilebilirliği.

Görüşler `REVIEW_coupling_scaling_v2.md` dosyasına, yalnızca hakemin kendisi tarafından kaydedilmelidir.

# Kaynaklar ve depo bağlantıları

**Depo kaynakları** (dal `claude/quantum-schrodinger-audit-iy6kn4`, HEAD `1668540`):

- `experiments/wz_rgf_diagnostic/PREDICTIONS_coupling_scaling_v2.md`: dondurulmuş şartname, Rev. 3.7 (`0f502c9`).
- `experiments/wz_rgf_diagnostic/R0_REVIEW_PACKET.md`: R0 inceleme paketi, Bölüm A–I (`1668540`).
- `docs/RGF_B1_VALIDATION.md` ve `results/rgf_b1_validation/{summary.json, phase_audit_hp_summary.json}`.
- `experiments/wz_rgf_diagnostic/REPORT.md`, `ERRATUM_worst_point.md` ve `results/` (scalar deney).
- `experiments/wz_rgf_diagnostic/amplitude_roundoff_report.md` ve `results_amplitude_roundoff/` (63 kanal).
- `experiments/wz_rgf_diagnostic/PREDICTIONS_coupled_blocks.md` ve `results_coupled_blocks/` (coupled-block).
- `experiments/wz_rgf_diagnostic/results_coupling_scaling_v2/` (manifestler; `fac5f5e`).
- Kod: `stationary.py`, `wz_rgf.py`, `wz_mutants.py`, `wz_mutants_oracle.py`, `wz_mutants_mp.py` (+ `_MATERIALS.md`,
  `_PROVENANCE.md`), `coupling_scaling_v2_{certify,covering,refs,manifest,analysis}.py`, `coupling_scaling_v2.py`,
  `cs2_t6_reference/`, `cs2_t6_harness/`, `cs2_t1t5/`.
- Grafikler: `docs/turkce_ozet_figurler/make_figures.py` (yalnızca commit edilmiş sonuç dosyalarını okur).

**Bilimsel kaynaklar** ([1], [2] ve [6] depodaki `docs/LITERATURE.md` listesinden alınmıştır; [3]–[5] bu belge için web
aramasıyla yayın kaydı üzerinden doğrulanmıştır):

1. P. W. Anderson, "Absence of diffusion in certain random lattices," *Phys. Rev.* **109**, 1492 (1958).
   doi:10.1103/PhysRev.109.1492
2. P. A. Lee, T. V. Ramakrishnan, "Disordered electronic systems," *Rev. Mod. Phys.* **57**, 287 (1985).
   doi:10.1103/RevModPhys.57.287
3. D. S. Fisher, P. A. Lee, "Relation between conductivity and transmission matrix," *Phys. Rev. B* **23**, 6851 (1981).
   doi:10.1103/PhysRevB.23.6851
4. P. A. Lee, D. S. Fisher, "Anderson localization in two dimensions," *Phys. Rev. Lett.* **47**, 882 (1981).
   doi:10.1103/PhysRevLett.47.882
5. A. MacKinnon, "The calculation of transport properties and density of states of disordered solids," *Z. Phys. B* **59**,
   385 (1985). doi:10.1007/BF01328846
6. B. Kramer, A. MacKinnon, "Localization: theory and experiment," *Rep. Prog. Phys.* **56**, 1469 (1993).
   doi:10.1088/0034-4885/56/12/001

::: {.dipnot}
Bu belge yalnızca depodaki dosyalardan derlenmiştir. Hazırlanması sırasında hiçbir çözücü, ana referans hesabı veya
driver çalıştırılmamış; manifestler, dondurulmuş şartname ve korunan dosyalar değiştirilmemiştir. Metin, bir yapay zekâ
kodlama asistanı (Claude Code) yardımıyla hazırlanmıştır.
:::
