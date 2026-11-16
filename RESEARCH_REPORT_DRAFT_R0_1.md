# Araştırma Raporu — Taslak R0.1

> **Taslak durumu.** Bu bir **iç revizyon taslağıdır**, yayıma hazır bir makale değildir. R0.1, R0 taslağının
> (`RESEARCH_REPORT_DRAFT_R0.md`) revizyonudur; R0 dosyası değiştirilmemiştir.
>
> - **Taslak notları.** Her ana bölümün sonunda "Taslak notu" başlığıyla kaynak dosyalar listelenmiştir. Bu notlar
>   revizyon sürecine yöneliktir ve nihai metinden çıkarılabilir.
> - **İşaretler.**
>   - **[TO BE COMPLETED]:** eksik bilgi;
>   - **[REFERENCE DETAILS TO VERIFY]:** bibliyografik ayrıntısı depodan doğrulanamayan kaynak;
>   - **[REVIEW]:** çözülmemiş bir tutarsızlık.
> - **Kapsam.** Bu taslak yazılırken hiçbir çözücü, eğitim, benchmark veya deney çalıştırılmamıştır. Tüm sayılar commit
>   edilmiş dosyalardan alınmıştır.
> - **Depo durumu** (R0.1 yazılırken):
>
>   | öğe | commit |
>   |---|---|
>   | dal | `claude/quantum-schrodinger-audit-iy6kn4` |
>   | rapor hazırlık denetimi (`REPORT_READINESS_AUDIT.md`) | `e7aa728` |
>   | R0 taslağı (`RESEARCH_REPORT_DRAFT_R0.md`) | `11b342c` (R0.1 yazılırken HEAD) |
>   | `origin/main` | `5bee7eb` (PR #26 birleştirmesi; dosya ağacı `11b342c` ile aynı) |
>   | R0.1 taslağı (bu dosya) | bu dosyayı ekleyen commit; kendi SHA'sını içeremez, commit mesajında ve teslim notunda verilir |
>
> - **R0 → R0.1 değişiklikleri** (ayrıntı: Ek B):
>   - TÜBİTAK 2204-A başlık kuralına uygun 11 kelimelik başlık;
>   - V6 ve PINN V4 değerleri kaynak dosyalarına bağlandı;
>   - kaynağı bulunamayan "kalıntı 0.22" iddiası çıkarıldı;
>   - özgün motivasyon, belgelenmiş depo geçmişi ile öğrencinin anlatımı ayrılarak yeniden yazıldı;
>   - CAP tanımı koda göre düzeltildi (s boyutsuz derinliktir);
>   - doğrulama türleri (iç doğrulama / yakınsama / bağımsız durağan doğrulama) ayrıldı;
>   - RGF/W-Z ayrıntıları Ek A'ya taşındı;
>   - Adam, Holm, Welch ve bootstrap p-değeri kaynakları eklendi.

---

## 1. Başlık

**Korelasyon Uzunluğunun İki Boyutlu Kuantum Dalga Paketi Taşınımına Etkisinin Sayısal İncelenmesi**

*Gayriresmî İngilizce karşılığı (çeviri; resmî başlık değildir):* Numerical Investigation of the Effect of Correlation
Length on Two-Dimensional Quantum Wave-Packet Transport

| | |
|---|---|
| Hazırlayan | [TO BE COMPLETED: Ad Soyad — depoda yok] |
| Okul / kurum | [TO BE COMPLETED — depoda yok] |
| Danışman | [TO BE COMPLETED — depoda yok] |
| Yarışma | TÜBİTAK 2204-A Lise Öğrencileri Araştırma Projeleri Yarışması, 2026–2027 [TO BE COMPLETED: alan/kategori bilgisi] |

> **İç not: TÜBİTAK başlık uyumu.**
>
> - **Kural.** Başlık en fazla 15 kelime olmalıdır.
> - **Kelime sayımı.** Kelimeler boşlukla ayrılarak sayılmıştır:
>   1. Korelasyon
>   2. Uzunluğunun
>   3. İki
>   4. Boyutlu
>   5. Kuantum
>   6. Dalga
>   7. Paketi
>   8. Taşınımına
>   9. Etkisinin
>   10. Sayısal
>   11. İncelenmesi
>
>   **Toplam: 11 kelime ≤ 15.**
> - **Uyum kontrolü.** Başlık, uygulanan çalışmayı tanımlar: ζ'nin iletime etkisi, 2B geometri ve sayısal yöntem.
>   - 3B bir model iddia etmez.
>   - PINN'in klasik çözücüden üstün olduğunu ima etmez; PINN başlıkta yer almaz.
>   - Anderson lokalizasyonunun kanıtlandığını ima etmez.
>   - Depo kanıtlarıyla ciddi bir bilimsel uyumsuzluk bulunmamıştır.
> - **Önceki başlık.** R0'daki uzun başlık (15 kelimeyi aşıyordu) resmî başlık olarak **kullanılmamaktadır.** Bu dosyada
>   15 kelimeyi aşan başka bir başlık resmî başlık olarak yer almaz.

---

## 2. Özet

Bu çalışma, iki boyutlu (2B) sert duvarlı bir dalga kılavuzunun sonlu bir bölgesine yerleştirilmiş Gauss-korelasyonlu
statik rastgele potansiyelin, bir kuantum dalga paketinin iletimini (*transmission*) ve yansımasını (*reflection*) nasıl
etkilediğini sayısal olarak inceler. Odak, düzensizliğin **korelasyon uzunluğu ζ**'dir; düzensizlik gücü W, gelen
momentum k ve düzensiz bölge uzunluğu L_d da taranmıştır.

**Yöntem.**

- Zamana bağlı Schrödinger denklemi (TDSE), tip-I ayrık sinüs dönüşümüne (DST-I) dayalı Strang split-step çözücüyle
  çözülmüştür. Uçlarda karmaşık soğurucu potansiyeller (CAP) kullanılmıştır.
- Çözücü analitik çözümlere, zaman ve uzay yakınsamasına ve olasılık korunumuna karşı sınanmıştır. Üretim ayarlarında
  taşınım gözlemlenebilirlerinin ayrıklaştırma hatası yaklaşık 10⁻⁵, istatistiksel belirsizlik ise yaklaşık 10⁻²
  mertebesindedir.

**Bulgular.**

- **W:** Ansambl ortalaması ⟨T⟩, W ile monoton azalır. 48 realizasyonla W = 25 için 0.938, W = 300 için 0.122 ölçülmüştür.
- **ζ:** Sabit W = 100'de ⟨T⟩, ζ'nin monoton olmayan bir fonksiyonudur ve kζ ≈ 1 civarında bir minimum gösterir. Bu
  minimum sonuçlar görüldükten sonra seçilen (*post-hoc*) eşleştirilmiş bootstrap kontrastlarıyla desteklenmiştir; ancak
  yalnızca 7 noktalı ızgara çözünürlüğünde belirlenmiştir.
- **Born karşılaştırması:** Parametresiz Born + difüzyon tahmini, kℓ_tr ≳ 10 için ölçülen ⟨T⟩'yi 0.95–1.12 oran
  aralığında yakalar. kℓ_tr ≈ 1–3'te oran 1.21–1.31'e çıkar; bu sapmanın mekanizması belirlenmemiştir.
- **L_d (W = 200):** ⟨T⟩'nin L_d ile azalması üstel bir biçime daha iyi uyar. Bu, lokalizasyonun başlangıcıyla
  **tutarlıdır ancak kanıt değildir.**
- **Kapalı kutu:** 2B kapalı kutularda erişilebilen boyutlarda lokalizasyon saptanmamıştır.

**PINN.** Aynı düzensiz potansiyel üzerinde değerlendirilen PINN, mevcut CPU bütçesinde taşınım fiziği için yeterli
doğruluğa ulaşamamıştır. En kolay ayarda bile son anda bağıl L2 hatası yaklaşık %49'dur; klasik çözücünün hatası ise
yaklaşık 10⁻⁶'dır.

**Doğrulama.** Çözücü iç testlerle (analitik çözümler, korunum) ve yakınsama testleriyle doğrulanmıştır. Bağımsız
sabit-enerji doğrulaması için başlatılan özyinelemeli Green fonksiyonu (RGF) hattı ise **tamamlanmamıştır:**

- B1-F: FAIL; B2: PASS; B3 ve B3X: dx kriterinde FAIL; B3Y: PASS.
- **B4 (düzensizlik içeren sabit realizasyonda TDSE ile durağan çözümün karşılaştırması) yapılmamıştır.**

Bu nedenle düzensizlik içeren taşınım sonuçları **bağımsız durağan bir yöntemle henüz doğrulanmamıştır** (Bölüm 9, Ek A).

---

## 3. Anahtar kelimeler

Kuantum taşınımı; korelasyonlu düzensizlik; zamana bağlı Schrödinger denklemi; dalga kılavuzu; split-step spektral
yöntem; Born yaklaşımı; Anderson lokalizasyonu (yalnızca bağlam olarak); fizik-bilgili sinir ağları; sayısal doğrulama;
özyinelemeli Green fonksiyonu

---

## 4. Giriş

### 4.1 Arka plan ve motivasyon

Düzensizlik birçok fiziksel sistemde bulunur: iletkenlerdeki safsızlıklar, dalga kılavuzlarındaki pürüzlülük, soğuk
atomlar için kullanılan rastgele optik (*speckle*) potansiyeller.

- **Girişim.** Rastgele bir potansiyelden defalarca saçılan dalgalar girişim yapar.
- **Rejimler.** Düzensizliğin gücüne, korelasyon uzunluğuna ve boyuta bağlı olarak taşınım üç biçimde olabilir:
  - balistik;
  - difüzif;
  - Anderson lokalizasyonuyla bastırılmış [1–3].
- **Korelasyon uzunluğunun rolü.** Korelasyonlu düzensizlikte tanecik boyutunun dalga boyuna oranı (kζ), saçılmanın ne
  kadar güçlü olduğunu ve hangi yönlere yöneldiğini belirler [9].

Zamana bağlı dinamiğin doğrudan simülasyonu iki şeye olanak verir: bu süreçleri izlemek ve iletim istatistiklerini çok
sayıda düzensizlik realizasyonu üzerinden ölçmek.

Karşılaştırma için kullanılan kuram yerleşiktir. Bu çalışmanın katkısı yeni bir kuram değildir. Sayısal deney,
pertürbatif öngörülerin belirli ve sonlu bir geometride nerede geçerli kaldığını sınar.

### 4.2 Kuantum taşınımı ve düzensiz sistemler

Aşağıdaki sonuçlar literatürden alınmıştır. Bunlar bu projenin bulguları değildir.

- **Boyut ve lokalizasyon.** Ölçekleme kuramına göre, ortogonal simetri sınıfında 1B ve 2B'de keyfi zayıf korelasyonsuz
  düzensizlik için bütün durumlar lokalizedir [2].
- **2B lokalizasyon uzunluğu.** Zayıf düzensizlikte 2B lokalizasyon uzunluğu ξ ≈ ℓ exp(πkℓ/2) kadar büyüktür [3]. Bu
  nedenle kℓ ~ 1 olmadıkça sonlu örneklerde gözlenemez.
- **Yarı-1B teller.** N açık kanallı bir telde lokalizasyon uzunluğu Nℓ mertebesindedir [6]. Rejimler:

  | Rejim | Koşul | Davranış |
  |---|---|---|
  | balistik | L ≪ ℓ | — |
  | difüzif (Ohm) | ℓ ≪ L ≪ Nℓ | ⟨T⟩ ∝ ℓ/L |
  | lokalize | L ≫ Nℓ | ln T ∝ −L/ξ; T geniş, log-normal dağılır [6] |

### 4.3 Araştırma sorusu

**Ana soru.**

> Gauss-korelasyonlu düzensizliğin uzaysal korelasyon uzunluğu ζ, sonlu iki boyutlu sert duvarlı bir dalga
> kılavuzundaki düzensiz bölgeden geçen kuantum dalga paketinin iletimini nasıl etkiler?

**Kapsam üzerine not.** Bu odak `docs/PHASE2_METHOD_SELECTION.md` girişindeki Phase 2 sorusudur. `README.md` §1'deki soru
daha geniştir:

- ζ'ye ek olarak düzensizlik gücü W'yu da kapsar;
- Born kuramıyla tutarlılığı sorar;
- lokalizasyonun sonlu-boyut etkilerinden ayırt edilip edilemeyeceğini sorar.

Bu taslak, proje başlığıyla tutarlı olarak ζ'yi ana soru alır; W, k ve L_d taramalarını destekleyici sonuçlar olarak
sunar. README'deki daha geniş formülasyon (W, Born tutarlılığı, lokalizasyon) bu destekleyici sonuçlar içinde ele
alınmıştır. Bu bir kapsam seçimidir; README'deki soruyla çelişmez.

**İkincil (yöntemsel) soru.**

> Aynı düzensiz 2B TDSE'nin PINN çözümü, doğrulanmış klasik çözücüye göre ne kadar doğrudur ve bu doğruluğu ne sınırlar?

### 4.4 Amaçlar ve kapsam

**Amaçlar:**

1. Doğrulanmış bir klasik çözücüyle ⟨T⟩ ve ⟨R⟩'nin ζ, W, k ve L_d'ye bağımlılığını güven aralıklarıyla ölçmek.
2. Gözlenen eğilimleri parametresiz bir Born + difüzyon tahminiyle karşılaştırmak.
3. Erişilebilir boyutlarda lokalizasyon belirtilerini sınamak.
4. PINN'in aynı problemdeki doğruluğunu ve hata kaynaklarını nicel olarak ortaya koymak.

**Kapsam dışı:**

- etkileşimler (Gross–Pitaevskii doğrusal olmayan terimi);
- zamana bağlı düzensizlik;
- faz bozulması (*dephasing*) veya açık sistem etkileri;
- üç boyut;
- deneysel doğrulama (README §10).

### 4.5 Özgün motivasyon ile uygulanan 2B çalışma arasındaki ilişki

Bu bölüm üç tür bilgiyi ayrı tutar: (a) depo geçmişinde **belgelenmiş** olanlar, (b) öğrencinin **anlatımı**
(belgelenmemiş), (c) **güncel** uygulanan çalışma.

**(a) Belgelenmiş depo geçmişi** (git geçmişinden; tarihler commit tarihleridir):

| tarih | commit | belgelenen içerik |
|---|---|---|
| 2026-07-18 | `00a52e9` | README: "Solving 2D Time-Dependent Schrödinger Equation using Physics-Informed Neural Networks (PINNs) with interactive 3D Plotly visualizations." |
| 2026-07-18 | `7b9f41a` | `wave_function.py`: \|Ψ(x,y,t)\|²'nin etkileşimli **3D yüzey görselleştirmesi**. Bu bir görselleştirmedir; çözülen denklem 2B'dir |
| 2026-07-18 | `5f56539` | README motivasyonu: "…to explore how machine learning can accelerate quantum chemistry and molecular biology simulations." Problem: 2B sonsuz kuyu içinde dikdörtgen bariyerden 2B Gauss paketi tünellemesi |
| 2026-08-08 | `4fe3eb2` | DST tabanlı klasik referans çözücü ve "ENAQT-esinli" düzensizlik deneyi eklenir |
| 2026-09-01 | `42cea06` | denetlenen sürüm: README başlığı "2D Tunneling & ENAQT-Inspired Disorder Experiment"; metin **1B** TDSE ve periyodik sınırlar anlatır, kod **2B** ve Dirichlet duvarlıdır (`docs/AUDIT.md` 1.1) |

**Bu sürümün denetim bulguları** (`docs/AUDIT.md`):

- ENAQT motivasyonu modelle uyuşmaz: modelde çevre, faz bozulması, Lindblad veya HEOM dinamiği yoktur (2.4).
- "Düzensizlik iletimi artırır" sonucu, kapalı kutuda yaklaşık 40 duvar yansımasından sonra alınmış bir anlık görüntünün
  yapaylığıdır (2.2–2.3; README §9.9).
- Özgün yapılandırmadaki referans çözüm yakınsamamıştı (README §7).

**(b) Öğrencinin anlatımı (belgelenmemiş).** Öğrenciye göre çalışmanın çıkış noktası, Schrödinger denkleminin daha
yüksek boyutlu uzaylarda makine öğrenmesi tabanlı yaklaşımlarla modellenmesine duyulan ilgidir.

- Depoda bu ifade **açıkça yer almaz.**
- En yakın belgelenmiş ifade, makine öğrenmesinin "kuantum kimyası ve moleküler biyoloji simülasyonlarını hızlandırması"
  hedefidir (`5f56539`).
- "3D" sözcüğü ilk commit'lerde yalnızca **görselleştirme** için kullanılmıştır.
- [TO BE COMPLETED: Öğrencinin anlatımı tarihli bir kaynakla, örneğin proje başvurusu, defter veya e-posta ile
  desteklenmedikçe, nihai metinde yalnızca "kişisel motivasyon" olarak sunulmalıdır.]

**(c) Güncel çalışma.**

- Statik düzensizlikte **üniter** dinamik üzerinedir. ENAQT çerçevesi kaldırılmıştır (README "Project history";
  `docs/LITERATURE.md` "Removed from the earlier bibliography").
- Uygulanan fiziksel çalışma **iki boyutludur.** Üç boyutlu Schrödinger denklemi çözülmemiştir.
- PINN bileşeni, makine öğrenmesi yöntemlerinin bu 2B problem üzerindeki **değerlendirmesidir.**
- Daha yüksek boyutlara genelleme, bu çalışmanın verisiyle desteklenmeyen bir gelecek çalışma konusudur.

**Evrimin özeti.** Proje, bir PINN gösterimi olarak başlamıştır. Daha sonra bağımsız bir klasik referans çözücü ve bir
düzensizlik deneyi eklenmiştir. Bir denetim, ilk sonuçların ve motivasyonun desteklenmediğini göstermiştir. Bunun
ardından odak, doğrulanmış klasik çözücüyle yapılan 2B taşınım çalışmasına kaymış, PINN ise ikincil bir değerlendirme
konumuna gelmiştir. Bu sıralama git geçmişiyle tutarlıdır. Ara adımların gerekçeleri ise yalnızca belgelerde yazıldığı
kadarıyla bilinmektedir.

> **Taslak notu (Bölüm 1–4).**
> - **Proje belgeleri:** `README.md` §1–§2, §10, "Project history"; `docs/PHASE2_METHOD_SELECTION.md` (giriş);
>   `docs/AUDIT.md` 1.1, 2.2–2.4; `docs/LITERATURE.md` §1 ve "Removed from the earlier bibliography".
> - **Depo geçmişi:** `git log --reverse`; `git show 00a52e9:README.md`, `git show 5f56539:README.md`,
>   `git show 4fe3eb2:README.md`, `git show 42cea06:README.md`; `git show 7b9f41a:wave_function.py`.
> - **Özet sayıları:** `results/transport/strength/summary.csv`; `results/transport/analysis.json`;
>   `results/pinn/summary.csv`.

---

## 5. Kuramsal arka plan

### 5.1 Zamana bağlı Schrödinger denklemi

Boyutsuz birimler kullanılmıştır (ħ = m = 1). Uzunluklar dalga kılavuzu genişliği birimindedir. Dalga fonksiyonu
ψ(x, y, t) şu denklemi sağlar:

$$ i\,\partial_t\psi = -\tfrac{1}{2}\left(\partial_x^2+\partial_y^2\right)\psi + V(x,y)\,\psi \qquad (1)$$

Burada V(x, y) statik potansiyeldir: düzensizlik, isteğe bağlı bir bariyer ve açık geometride soğurucu terimler.

### 5.2 Dalga fonksiyonu ve olasılık yoğunluğu

ρ = |ψ|² olasılık yoğunluğudur. Soğurucu katman yoksa dinamik üniterdir ve ∫|ψ|² dx dy = 1 korunur.

Soğurucu katmanlar varsa, her uçta soğurulan olasılık ayrıca kaydedilir. Böylece T + R + P_in = 1 yuvarlama hatası
düzeyinde sağlanır (README §3). Burada P_in, t_end anında hâlâ düzensiz bölgede bulunan olasılıktır.

### 5.3 Dalga paketi yayılımı

Gelen durum, enine 1. modda bir Gauss paketidir:

$$ \psi_0 = (\pi\sigma_x^2)^{-1/4} e^{-(x-x_0)^2/2\sigma_x^2} e^{i k_{x0} x}\,\sqrt{2/L_y}\,\sin(\pi y/L_y) \qquad (2)$$

- Paketin momentum yayılımı σ_k = 1/(√2 σ_x)'tir.
- Paket tek bir enerjiyi değil, bir **enerji bandını** örnekler. Ölçülen T ve R bu nedenle enerji üzerinden
  ortalanmıştır; mezoskopik kuramın sabit-enerji nicelikleri değildir (README §10).

### 5.4 Sert duvarlı dalga kılavuzu

- **Kılavuz.** y ∈ [0, L_y], L_y = 1 aralığında ψ = 0 koşulu (sert duvar) uygulanır.
- **Enine modlar ve enerjiler:** χ_n(y) = √(2/L_y) sin(nπy/L_y), E_⊥,n = (nπ/L_y)²/2.
- **Açık kanallar.** E₀ > E_⊥,n koşulunu sağlayan kanallardır. Üretim ayarında (k_x0 = 20) toplam dalga sayısı
  k = (k_x0² + π²)^{1/2} ≈ 20.25 ve E₀ = k²/2 ≈ 205'tir; bu da **N = 6 açık enine kanal** verir (README §3).

### 5.5 Gauss-korelasyonlu düzensizlik

Düzensizlik, bir tohum (*seed*), RMS gücü W ve korelasyon uzunluğu ζ ile tanımlanan **sürekli** bir rastgele alandır.
Shinozuka–Deodatis spektral temsil yöntemiyle üretilir [14]:

$$ V_{\mathrm{dis}}(\mathbf r) = W\sqrt{\frac{2}{M}}\sum_{j=1}^{M}\cos(\mathbf q_j\cdot\mathbf r+\varphi_j)\,S(x),\qquad \mathbf q_j\sim\mathcal N(0,\zeta^{-2}\mathbb{1}_2),\ \varphi_j\sim U[0,2\pi) \qquad (3)$$

**Parametreler:**

- M = 512 mod kullanılır.
- S(x) = ½[tanh((x−a)/δ) − tanh((x−b)/δ)] (δ = 0.02) düzensizliği a < x < b bölgesiyle sınırlayan yumuşak bir
  penceredir.

**İstatistiksel özellikler.** Her M değeri için ansambl ortalamasında:

$$ \langle V_{\mathrm{dis}}\rangle = 0,\qquad \langle V_{\mathrm{dis}}(\mathbf r)V_{\mathrm{dis}}(\mathbf r')\rangle = W^2 e^{-|\mathbf r-\mathbf r'|^2/2\zeta^2},\qquad S(q) = 2\pi W^2\zeta^2 e^{-q^2\zeta^2/2}. \qquad (4)$$

**Realizasyon düzeyinde:**

- Tek tek realizasyonlar yeniden normalize edilmez.
- Tek-nokta istatistiği yalnızca yaklaşık olarak Gauss'tur (README §10).
- Aynı tohum, her ızgara çözünürlüğünde aynı fiziksel alanı verir.
- NumPy (referans çözücü) ve PyTorch (PINN) hesapları noktasal olarak 10⁻⁹ bağıl düzeyde uyuşur (README §4;
  `tests/test_potential.py`).

### 5.6 İletim ve yansıma

Açık dalga kılavuzunda gözlemlenebilirler t_end anında tanımlanır (`transport.py`; README §8). a ve b düzensiz bölgenin
sınırlarıdır.

| Nicelik | Tanım |
|---|---|
| iletim T | sağ uçta soğurulan olasılık + P(x > b) |
| yansıma R | sol uçta soğurulan olasılık + P(x < a) |
| P_in | 1 − T − R |

P_in küçük olduğunda T ve R, paketin enerji yayılımı üzerinden ortalanmış asimptotik saçılma olasılıklarına yaklaşır.

### 5.7 İlgili fiziksel ölçekler ve boyutsuz parametreler

**Born ortalama serbest yolları.** Sınırsız 2B Gauss-korelasyonlu düzensizlikte k dalga sayılı bir parçacık için Fermi'nin
altın kuralı aşağıdaki oranları verir (`disorder.born_rates`; README §4). Burada a = k²ζ² ve I₀, I₁ değiştirilmiş Bessel
fonksiyonlarıdır.

$$ \frac{1}{\tau_s}=2\pi W^2\zeta^2 e^{-a}I_0(a),\qquad \frac{1}{\tau_{tr}}=2\pi W^2\zeta^2e^{-a}\left[I_0(a)-I_1(a)\right],\qquad \ell_{s,tr}=k\,\tau_{s,tr} \qquad (5)$$

**Difüzif tahmin.** Bölge iletimi için basit bir difüzif tahmin kullanılır:

$$ T_{\mathrm{est}} \approx \frac{\ell}{L_d+\ell}, \qquad \ell = \frac{\pi}{2}\,\ell_{tr}. \qquad (6)$$

**Geçerlilik koşulları:**

- Bu ifadeler yalnızca zayıf düzensizlik (kℓ ≫ 1) için geçerlidir.
- Dalga kılavuzunun ayrık kanal yapısını içermezler.
- π/2 öneki konvansiyona bağlıdır.
- Bu yüzden **mertebe kıyası** olarak kullanılmalı ve yalnızca eğilimleri yorumlanmalıdır (README §4, §9.2).

**Boyutsuz parametreler:**

| Parametre | Anlamı |
|---|---|
| kζ | tanecik boyutu / dalga boyu |
| kℓ_tr | saçılma gücü; kℓ ≫ 1 zayıf saçılmadır |
| W/E₀ | potansiyel genliğinin kinetik enerjiye oranı |
| L_d/ℓ ve L_d/ξ | sistem uzunluğunun taşınım ölçeklerine oranı |

**Sayısal örnekler** (`results/transport/analysis.json`, alan `k_l_tr`):

- W = 100, ζ = 0.05'te kℓ_tr = 10.4;
- W = 300'de kℓ_tr = 1.16.

> **Taslak notu (Bölüm 5).**
> - `README.md` §3, §4, §8, §10;
> - `config.py` (`waveguide_config`);
> - `disorder.py` (`born_rates`);
> - `transport.py`;
> - `results/transport/analysis.json` (`k_l_tr`).

---

## 6. Fiziksel model ve sayısal yöntemler

### 6.1 Geometri ve sınır koşulları

| Geometri | Kullanım | Alan | Başlangıç durumu | Sınırlar |
|---|---|---|---|---|
| Açık dalga kılavuzu (`waveguide_config`) | taşınım deneyleri | L_y = 1, L_x = 2.9 + L_d | mod-1 Gauss, x₀ = 1.2, σ_x = 0.2, k_x0 = 20 | y'de sert duvar; x uçlarında 0.6 genişlikte CAP |
| Kapalı kutu (`closed_box_config`) | yayılma / sonlu-boyut testi | L × L, L ∈ {1, 2} | 2B Gauss, merkezde, σ = 0.08, k_x0 = 20 | sert duvarlar; düzensizlik her yerde |
| PINN problemi (`pinn_config`) | PINN ile referans karşılaştırması | 1 × 1, t ∈ [0, 0.02] | 2B Gauss, (0.35, 0.5), σ = 0.07, k_x0 = 15 (E₀ = 112.5) | sert duvarlar |

**Dalga kılavuzu yerleşimi** (x ekseni boyunca):

1. [0, 0.6] sol CAP;
2. x₀ = 1.2 paket merkezi;
3. [1.9, 1.9 + L_d] düzensiz bölge;
4. 0.4 boşluk;
5. 0.6 sağ CAP.

**Sayısal ayarlar:** ızgara dx = dy = 1/64, ny = 63, Δt = 2×10⁻⁴.

**Bitiş zamanı t_end** (`summary.csv`, sütun `t_end`):

| Sweep | t_end |
|---|---|
| strength, ζ, convergence | 1.0 |
| momentum | k = 12 için 1.667, k = 16 için 1.25, diğerleri için 1.0 |
| length | 1 + 0.5 L_d, yani 1.125–2.5 |
| barrier | 1.5 |

### 6.2 Düzensizlik üretimi

Denklem (3) kullanılır (`disorder.py`). Realizasyon i, her parametre noktasında tohum i'yi kullanır (README §8).

**Eşleştirilmiş tasarımın sonuçları:**

- **Avantaj.** Farklı parametre noktaları aynı rastgele çekilişleri paylaşır; W ve ζ ile yalnızca ölçeklenirler. Bu,
  komşu noktalar arasındaki farkların gürültüsünü azaltır.
- **Bedel.** Bir tarama içindeki noktalar **istatistiksel olarak bağımsız değildir.**

### 6.3 Başlangıç dalga paketi

Dalga kılavuzu için denklem (2), kapalı kutu ve PINN problemi için 2B Gauss paketi kullanılır:

$$ \psi_0 = (\pi\sigma_x\sigma_y)^{-1/2}\exp\!\left[-\frac{(x-x_0)^2}{2\sigma_x^2}-\frac{(y-y_0)^2}{2\sigma_y^2}\right]e^{i(k_{x0}x+k_{y0}y)}. $$

Her iki paket de 1'e normalizedir; testler kutu üzerinde < 10⁻⁶ düzeyinde doğrular (README §3; `initial_condition.py`).

### 6.4 Hamiltonyen ve TDSE

H = −½∇² + V'dir. Potansiyel üç bileşenden oluşur:

- V_dis: düzensizlik;
- V_b: isteğe bağlı bariyer; tanh kenarlı, genişliği 0.1 ve kenar genişliği 0.01 olan bir şerittir;
- −iW_CAP(x): açık geometride soğurucu terim.

Düzensizlik ve bariyer tanımı iki çözücü tarafından ortak kullanılır (`potential.py`).

### 6.5 DST tabanlı Strang split-step çözücü

**Şema** (`reference_solver.py`; README §5) [15–17]:

- **Uzaysal temsil.** 2B tip-I ayrık sinüs dönüşümü kullanılır.
- **Kinetik adım.** sin(nπx/L) modları Dirichlet Laplasyeninin özfonksiyonlarıdır. Bu yüzden kinetik adım sinüs serisi
  temsili için tamdır ve duvarlarda ψ = 0 kendiliğinden sağlanır.
- **Zaman adımı.** Strang ayrıştırması kullanılır: Δt'de ikinci dereceden, uzayda spektral doğruluktadır.

**Soğurucu katman.** Denklem (1)'e −iW_CAP(x)ψ terimi eklenir. W_CAP(x) = W₀ s(x)² biçimindedir; W₀ = 300'dür. s(x)
boyutsuz derinliktir: katman içindeki derinliğin katman genişliğine (0.6) oranıdır, yani 0 ≤ s ≤ 1
(`reference_solver.py::cap_profile`; `cfg.cap_strength = 300`, `cap_power = 2`, `cap_width = 0.6`) [19]. Soğurulan
olasılık her uçta ayrı ayrı biriktirilir.

**Diğer özellikler.** Birden fazla düzensizlik realizasyonu tek bir yığın (*batch*) olarak yayılır. Referans çözücü
**kesin varsayılmaz**; doğruluğu ölçülmüştür (§6.6).

### 6.6 Sayısal doğrulama: üç ayrı tür

Bu raporda üç doğrulama türü birbirinden ayrı tutulur:

1. **İç çözücü doğrulaması** (*internal verification*): çözücünün analitik çözümleri ve korunum yasalarını doğru
   uygulayıp uygulamadığı. Aynı çözücü kendi ürettiği sonuçlarla değil, kesin bilinen çözümlerle sınanır.
2. **Yakınsama testleri:** ayrıklaştırma parametreleri (Δt, dx, M) inceltildiğinde gözlemlenebilirlerin değişimi. Bu
   testler ayrıklaştırma hatasının büyüklüğünü verir. Aynı denklemin farklı bir yöntemle çözümüne karşı bir sınama
   değildir.
3. **Bağımsız durağan doğrulama:** kod paylaşmayan, sabit enerjide çalışan bir yöntemle (RGF) karşılaştırma
   (Bölüm 9). Düzensizlik içeren taşınım için gerekli olan karşılaştırma B4'tür ve **yapılmamıştır.**

Tablo 6.1, 6.2 ve 6.3'ün kaynakları: `results/validation/validation_summary.json`, `v1_eigenmode.csv` …
`v7_n_modes.csv`; `figures/validation_solver.png`; üretici betik `validate_solver.py`.

**(1) İç çözücü doğrulaması**

| Kontrol | Ayar | Sonuç |
|---|---|---|
| V1 tam özmod (sabit V = 50) | 63², Δt = 10⁻⁴, t ≤ 0.1 | max \|hata\| = 4.58×10⁻¹⁴ (`v1_eigenmode_max_abs_err`) |
| V2 analitik serbest Gauss | 127² ve 255² | bağıl L2 ≤ 9.39×10⁻¹⁰ (`v2_free_gaussian_max_rel_l2`) |
| V6 korunum, kapalı kutu | PINN problemi, W = 112.5 (= E₀), CAP yok | max \|N(t) − 1\| = 1.51×10⁻¹² (`v6.closed_norm_max_dev`) |
| V6 korunum, açık kılavuz | W = 200, t ≤ 1 | max \|N(t) + soğurulan_sol + soğurulan_sağ − N(0)\| = 1.60×10⁻¹³ (`v6.open_bookkeeping_max_dev`) |

: Tablo 6.1. V6'nın tanımı `validate_solver.py::v6_conservation` fonksiyonundadır. Aynı değerler
`results/validation_log.txt` ve `results/reproduction/reproduction_log.txt` dosyalarında da kayıtlıdır. V6 için ayrı bir
CSV dosyası yoktur.

**(2) Yakınsama testleri**

| Kontrol | Ayar | Sonuç |
|---|---|---|
| V3 Δt yakınsaması, PINN problemi (W = E₀) | 127², Δt = 1.25×10⁻⁶'ya göre | 2. derece; üretim Δt = 10⁻⁵'te bağıl L2 5.29×10⁻⁷ (`v3`, `case = pinn_W112`) |
| V3 Δt yakınsaması, dalga kılavuzu (W = 200) | dx = 1/64, t = 1 | 2. derece; üretim Δt = 2×10⁻⁴'te \|ΔT\| = 1.96×10⁻⁵, \|ΔR\| = 2.25×10⁻⁵ (`v3`, `case = waveguide_W200`) |
| V4 uzaysal yakınsama, dalga kılavuzu (W = 200; ζ = 0.02 ve 0.05) | dx = 1/32, 1/64, 1/128 | dx = 1/64'te \|ΔT\| ≤ 1.4×10⁻⁵; dx = 1/32'de ≤ 2×10⁻³ (README §7; ham T değerleri `v4`, `case = waveguide_W200_zeta*`) |
| V4 uzaysal yakınsama, PINN problemi (W = E₀) | dx = 1/32, 1/64, 1/128 | dx = 1/128'de (127² ızgara) bağıl L2 1.90×10⁻⁶ (`v4_space_convergence.csv`, `case = pinn_W112_psi`; bu satırda hata değeri `T` sütununda saklanır) |
| V5 CAP sahte geri dönüşü (temiz kılavuz) | k_x0 = 10 … 30 | k = 20: 1.0×10⁻⁴; k = 16: 1.3×10⁻³; **k = 12: 1.33×10⁻²; k = 10: 3.07×10⁻²** (`v5_cap_reflection.csv`, `spurious_return`) |
| V7 düzensizlik mod sayısı M | W = 100, 24 realizasyon | M = 64: 0.488 ± 0.016; 256: 0.515 ± 0.018; 512: 0.518 ± 0.015; 2048: 0.501 ± 0.014; anlamlı bağımlılık yok (`v7_n_modes.csv`, `T_mean`, `T_se`) |

: Tablo 6.2. V5 bir yakınsama testi değil, soğurucu sınırın sistematik hatasının ölçümüdür; burada sayısal parametrelerin
sınanması olarak listelenmiştir.

**(3) Bağımsız durağan doğrulama** (ayrıntı Bölüm 9 ve Ek A)

| Karşılaştırma | Kapsam | Durum |
|---|---|---|
| B2: TDSE paket iletimi ile enerji-ortalamalı bağımsız 1B ODE | düzensizlik **yok**; y'den bağımsız yumuşak bariyer (V_b = 150, 250) | uzaklık 3.2×10⁻⁴ / 7.0×10⁻⁴; kriter 10⁻³; geçti (`docs/RGF_B2_VALIDATION.md`) |
| B4: TDSE ile durağan RGF, sabit **düzensizlik** realizasyonunda | Bölüm 7'deki taşınım sonuçlarının bağımsız doğrulaması | **YAPILMADI** |

: Tablo 6.3.

**Değerlendirme:**

- Üretim ayarlarında taşınım gözlemlenebilirlerinin ayrıklaştırma hatası yaklaşık 10⁻⁵'tir. Bu, sonlu ansambllardan
  gelen istatistiksel belirsizliğin (SE ~10⁻²) yaklaşık üç mertebe altındadır.
- İhmal edilemeyen tek sistematik hata, **yavaş bileşenler için CAP'tır.** k_x0 = 12'de olasılığın yaklaşık %1.3'ü sahte
  olarak geri döner; bu, o noktada R'yi yaklaşık 0.01 yukarı kaydırır.
- **Bağımsız doğrulamanın sınırı.** Bağımsız bir yöntemle yapılan tek TDSE karşılaştırması (B2) düzensizlik içermez.
  Düzensiz taşınım sonuçları yalnızca iç doğrulama ve yakınsama testlerine dayanır; **B4 tamamlanmış değildir.**

### 6.7 PINN formülasyonu ve eğitim yöntemi

**PDE kalıntıları.** ψ = u + iv yazıldığında denklem (1) iki gerçel PDE'ye ayrılır:

$$ f_u \equiv \partial_t u + \tfrac{1}{2}\nabla^2 v - V v = 0, \qquad f_v \equiv \partial_t v - \tfrac{1}{2}\nabla^2 u + V u = 0. \qquad (7)$$

**Ağ:**

- Mimari: 4 gizli katman × 64 tanh birim.
- Girdiler [0, 1]'e normalize edilir: x/L_x, y/L_y, t/T.
- Ardından isteğe bağlı bir Fourier öznitelik gömmesi gelir: birim başına 1'den 6'ya 6 doğrusal aralıklı frekans [23].

**Sert kısıtlar:**

- u = u₀ + A·B(x,y)·τ·N_u ve v = v₀ + A·B·τ·N_v; burada B = 16x̂(1−x̂)ŷ(1−ŷ) ve τ = t/T.
- Başlangıç koşulu tam olarak sağlanır.
- Dirichlet koşulu duvarlarda |ψ₀| kadar sağlanır (< 4×10⁻⁵).

**Kayıp fonksiyonu.** Ölçeklenmiş kalıntının ortalama karesi, (T/A)·(f_u, f_v), artı rastgele zamanlarda bir normalizasyon
cezası, (∫|ψ|² − 1)².

**Kolokasyon:**

- Her Adam adımında 4000 noktalık taze, düzgün dağılımlı bir yığın çekilir;
- alternatif olarak kalıntı-uyarlamalı örnekleme (RAD) kullanılır [28];
- isteğe bağlı nedensel zaman ağırlıklandırması uygulanır [27].

**Optimizasyon:**

- 2000 Adam adımı. Öğrenme oranı 10⁻³'tür ve her 1000 adımda yarıya iner. PDE ağırlığı 300 adımda doğrusal olarak
  ısıtılır; gradyan kırpma uygulanır.
- Adam algoritmasının kaynağı [37]'dir.
- Ardından sabit 16000 noktada 200 L-BFGS iterasyonu yapılır [30].

**Potansiyel.** Referans çözücüyle aynı `Potential` nesnesidir, yani aynı düzensizlik realizasyonu kullanılır.

**Eğitim protokolü bulgusu** (bir ayar seçimi değil, bir gözlem):

- İlk reçetede sabit 4000 kolokasyon noktası vardı ve eğitim sabit noktalarda L-BFGS ile bitiyordu.
- Bu reçetede son eğitim kaybı 9.07×10⁻⁴'tür (`results/pinn/diag_fixed_points_lbfgs/train_info.json`, `final_loss`).
  Buna karşın çözüm hatası 1.218'dir (`eval_summary.json`, `rel_l2_final`); bu, ψ = 0 tahmininden daha kötüdür.
- Ağ, kolokasyon noktalarına **aşırı uyum** sağlamıştır: düşük eğitim kaybı, düşük çözüm doğruluğuyla birlikte
  görülmüştür.
- Bu koşu için ayrılmış noktalarda (*held-out*) bir kalıntı değeri saklanmamıştır (`train_info.json`'da `holdout_pde`
  alanı yoktur). Bu nedenle bu koşu için ayrılmış nokta kalıntısı raporlanmamaktadır (Ek B, çıkarılan iddialar).
- Her adımda taze yığın kullanıldığında eğitim ve ayrılmış (*held-out*) kalıntılar benzer çıkar: 5.67×10⁻³ ve 9.28×10⁻³
  (`results/pinn/diag_stochastic/train_info.json`, `final_train_pde`, `holdout_pde`). Tüm çalışma koşularında bu protokol
  kullanılmıştır.
- Nedensel ağırlıklandırma (ε = 1) bu bütçede yardımcı olmamıştır: son hata 0.565, nedensel ağırlıklandırma olmadan 0.4885
  (`diag_stochastic_causal`, `diag_stochastic`; `eval_summary.json`, `rel_l2_final`).

> **Taslak notu (Bölüm 6).**
> - **Proje belgeleri:** `README.md` §3, §5–§8.
> - **Kod:** `config.py`, `reference_solver.py`, `potential.py`, `initial_condition.py`, `transport.py`, `model.py`,
>   `pinn_core.py`, `train.py`.
> - **Doğrulama:** `results/validation/*`.
> - **PINN teşhis koşuları:** `results/pinn/diag_*/{train_info,eval_summary}.json`.
> - **Transport t_end:** `results/transport/*/summary.csv` (`t_end`).

---

## 7. Sonuçlar

**Ortak bilgiler.**

- **Ham veriler.** Tüm taşınım sonuçları DST referans çözücüsüyle üretilmiştir (`experiment_disorder.py`). Ham veriler
  `results/transport/<sweep>/<point>.npz`, istatistikler `summary.csv` dosyalarındadır.
- **Güven aralıkları.** Aksi belirtilmedikçe ansambl ortalamasının %95 percentile-bootstrap aralıklarıdır (sütunlar
  `T_ci_boot_lo`, `T_ci_boot_hi`).
- **Bağımsızlık.** Aralıklar her nokta için ayrı ayrı geçerlidir. Bir sweep içindeki noktalar eşleştirilmiş tasarım
  nedeniyle bağımsız değildir.
- **Analizin niteliği.** Bölüm 7'deki taramaların tümü **keşifseldir** (*exploratory*). Önceden kayıtlı doğrulayıcı bir
  kampanya önerilmiş, ancak dondurulmamış ve çalıştırılmamıştır (`docs/PHYSICAL_CAMPAIGN_PLAN.md`, revizyon 3).

### 7.1 Düzensizlik gücüne bağımlılık

Ayarlar: ζ = 0.05, L_d = 1, k_x0 = 20, t_end = 1, n = 48 realizasyon (W = 0 için tek deterministik koşu).

| W | W/E₀ | kℓ_tr (Born) | ⟨T⟩ [%95 CI] | ⟨R⟩ | sd(T) | P_in(t_end) |
|---|---|---|---|---|---|---|
| 0 | 0 | ∞ | 1.000 | 0.000 | – | 7.9×10⁻⁷ |
| 25 | 0.12 | 167 | 0.938 [0.935, 0.941] | 0.060 | 0.011 | 2.0×10⁻³ |
| 50 | 0.24 | 42 | 0.792 [0.783, 0.801] | 0.207 | 0.032 | 1.1×10⁻³ |
| 75 | 0.37 | 18.5 | 0.629 [0.614, 0.645] | 0.370 | 0.054 | 7.9×10⁻⁴ |
| 100 | 0.49 | 10.4 | 0.494 [0.475, 0.513] | 0.505 | 0.067 | 7.3×10⁻⁴ |
| 150 | 0.73 | 4.6 | 0.316 [0.295, 0.337] | 0.683 | 0.074 | 1.2×10⁻³ |
| 200 | 0.98 | 2.6 | 0.222 [0.202, 0.242] | 0.776 | 0.071 | 1.6×10⁻³ |
| 300 | 1.46 | 1.2 | 0.122 [0.104, 0.141] | 0.875 | 0.066 | 3.2×10⁻³ |

: Tablo 7.1. Kaynak: `results/transport/strength/summary.csv` (sütunlar `T_mean`, `T_ci_boot_lo`, `T_ci_boot_hi`,
`R_mean`, `T_std`, `P_in_mean`); kℓ_tr için `analysis.json["strength"][*].k_l_tr`.

**Gözlem:**

- ⟨T⟩, incelenen aralıkta W ile monoton ve düzgün biçimde azalır.
- t_end anında olasılığın neredeyse tamamı bölgeyi terk etmiştir (P_in ≤ %0.32). Bu nedenle T ve R asimptotik değerlerine
  yakındır.

**Önceki sürümle ilişki.** Bu sonuç, ilk sürümün "düzensizlik iletimi artırır (0.424 vs 0.277)" iddiasını geçersiz kılar.
O iddia, kapalı kutuda alınmış tek bir anlık görüntünün yapaylığıydı (README §9.9; `docs/AUDIT.md` 2.3).

**Born karşılaştırması** (`results/transport/analysis.json`, alanlar `ratio`, `ratio_avg`, `k_l_tr`):

- Paketin momentum dağılımı üzerinden ortalanmış parametresiz Born + difüzyon tahmini (6), kℓ_tr ≥ 10 olan noktalarda
  ⟨T⟩_sim/T_est oranını **0.950–1.121** aralığında verir (strength, ζ ve momentum taramalarının tümü).
- Oran, kℓ_tr azaldıkça **sistematik olarak artar.** kℓ_tr ≈ 1–3'te 1.21–1.31'e ulaşır:

  | Nokta | kℓ_tr | oran |
  |---|---|---|
  | W = 200 | 2.6 | 1.213 |
  | W = 300 | 1.16 | 1.309 |
  | k_x0 = 12 | 1.71 | 1.303 |

- Enerji ortalaması sapmanın bir kısmını açıklar (W = 300'de 1.48'den 1.31'e), tamamını açıklamaz.
- **Yorum.** Bu sapmanın mekanizması belirlenmemiştir. Aday nedenler:
  - Born yaklaşımının çökmesi;
  - dalga kılavuzunun 6 kanallı ayrık yapısı;
  - ℓ/(L+ℓ) biçiminin kabalığı.

  Bunların hiçbiri test edilmemiştir.

### 7.2 Korelasyon uzunluğuna bağımlılık (ana soru)

Ayarlar: W = 100, L_d = 1, k_x0 = 20 (k ≈ 20.25), t_end = 1, n = 48.

| ζ | kζ | kℓ_tr (Born) | ⟨T⟩ [%95 CI] | ⟨R⟩ | sd(T) |
|---|---|---|---|---|---|
| 0.02 | 0.40 | 20.8 | 0.670 [0.656, 0.684] | 0.329 | 0.050 |
| 0.03 | 0.61 | 12.4 | 0.545 [0.526, 0.565] | 0.454 | 0.068 |
| 0.05 | 1.01 | 10.4 | 0.494 [0.475, 0.513] | 0.505 | 0.067 |
| 0.07 | 1.42 | 14.4 | 0.535 [0.506, 0.562] | 0.465 | 0.100 |
| 0.10 | 2.02 | 24.0 | 0.590 [0.560, 0.619] | 0.409 | 0.106 |
| 0.15 | 3.04 | 38.9 | 0.694 [0.650, 0.735] | 0.306 | 0.153 |
| 0.20 | 4.05 | 53.0 | 0.734 [0.687, 0.777] | 0.265 | 0.159 |

: Tablo 7.2. Kaynak: `results/transport/corr_length/summary.csv`; kℓ_tr için `analysis.json["corr_length"]`. kζ
değerleri k = 20.245 ile hesaplanmıştır (`paired_corr_length.json`, alan `k`).

**Gözlemler:**

1. ⟨T⟩ ζ'nin **monoton olmayan** bir fonksiyonudur ve kζ ≈ 1'de en küçük değerini alır. Uçlardaki ve ortadaki noktaların
   %95 aralıkları örtüşmez.
2. **Örnekten örneğe dağılım ζ ile büyür:** sd(T) 0.050'den 0.159'a çıkar.

**Yorum.**

- Gözlem 1, Born öngörüsüyle niteliksel olarak uyumludur. Born öngörüsüne göre ℓ_tr, kζ ≈ 0.8 civarında en küçüktür:
  - çok kısa menzilli düzensizliğin toplam gücü zayıftır (beyaz gürültü sınırı, ∝ ζ²);
  - çok düzgün düzensizlik ise ağırlıklı olarak ileri yönde saçar.

  Bu açıklama yerleşik kurama aittir (README §9.3).
- Yumuşak düzensizlikte (kζ ≳ 2) tahmin T'yi yaklaşık %5 fazla öngörür (ratio_avg ≈ 0.95–0.97).
- Gözlem 2 için önerilen açıklama: bölgedeki bağımsız korelasyon alanlarının sayısı ∝ L_d L_y/ζ² azalır ve iletim daha az
  öz-ortalanır. Bu bir **yorumdur**, ayrıca test edilmemiştir.

**Eşleştirilmiş yeniden analiz (POST-HOC).** ζ noktaları tohumları paylaştığı için karşılaştırma, tohum bazında
eşleştirilmiş bootstrap farklarıyla yinelenmiştir.

- Ayarlar: `paired_analysis.py`, 48 tohum, B = 20000 bootstrap, tohum 20261105; üç birincil kontrast üzerinde Holm
  düzeltmesi [38]; bootstrap p-değeri (b+1)/(B+1) biçimindedir [40, 41].
- Kontrastlar, eşleştirilmemiş sonuçlar **görüldükten sonra** seçilmiştir. Bu bağımsız bir tekrar değildir.

| Kontrast | Tanım | Ortalama fark | Eşleştirilmiş %95 bootstrap CI | Holm-düzeltilmiş bootstrap p |
|---|---|---|---|---|
| C1 | ⟨T⟩(0.05) − ⟨T⟩(0.02) | −0.176 | [−0.198, −0.154] | 1.5×10⁻⁴ |
| C2 | ⟨T⟩(0.05) − ⟨T⟩(0.20) | −0.240 | [−0.284, −0.193] | 1.5×10⁻⁴ |
| C3 | eğrilik K = [T(0.03) + T(0.07)]/2 − T(0.05) | +0.046 | [0.023, 0.068] | 1.5×10⁻⁴ |
| C3b | ⟨T⟩(0.07) − ⟨T⟩(0.05) (destekleyici) | +0.041 | [0.006, 0.074]; %98.33 CI [−0.002, 0.081] | — |

: Tablo 7.3. Kaynak: `results/analysis/paired_corr_length.json` (`contrasts.C1`…`C3b_T0.07_minus_T0.05`; alanlar `mean`,
`boot_ci_0.9500`, `boot_ci_0.9833`, `p_boot_holm`). Veri girdisi: `results/reproduction/transport/corr_length` (alan
`data`). Bu dizinin `results/transport/corr_length` ile bit düzeyinde aynı olduğu `docs/REPRODUCIBILITY_AUDIT.md` §5'te
belgelenmiştir.

**Analizin çıktıları:**

- **p değerleri.** Bootstrap kuyruk tahmini (b+1)/(B+1), çözünürlük sınırındadır (≈ 5×10⁻⁵; Holm sonrası 1.5×10⁻⁴). Bu,
  tam bir sonlu-örneklem p değeri değildir (README §9.3).
- **Minimumun konumu.** Bootstrap'ta minimum 98.9% oranında ζ = 0.05'te, 1.1% oranında ζ = 0.07'de bulunur
  (`argmin_distribution`).
- **Parabolik tepe kestirimi.** kζ = 0.95 [0.88, 1.12] (`vertex_k_zeta`).
- **Çoklu karşılaştırma.** ζ = 0.07 ile fark, düzeltilmiş aralıkta sıfırı içerir. Minimum bu nedenle yalnızca **kaba
  biçimde** (kζ ≈ 0.6–1.4) belirlenmiştir.

**Sonuç (sınırlamalarıyla).**

- Bu geometri ve paket için, W = 100'de ⟨T⟩(ζ)'nin kζ ≈ 1 civarında bir minimumu olduğu verilerle desteklenmektedir.
- Kanıt keşifsel ve post-hoc düzeydedir.
- Doğrulayıcı bir test, önceden kaydedilmiş fiziksel kampanya ile yapılabilirdi. Bu kampanya planlanmış, ancak
  yapılmamıştır (§10, §12).

### 7.3 Momentuma bağımlılık

Ayarlar: W = 100, ζ = 0.05, L_d = 1, n = 32. Her k için temiz kılavuz referansı (W = 0) da vardır.

| k_x0 | ⟨T⟩ [%95 CI] | ⟨R⟩ | kℓ_tr | temiz kılavuz T / R | t_end |
|---|---|---|---|---|---|
| 12 | 0.251 [0.222, 0.280] | 0.743 | 1.71 | 0.984 / 0.014 | 1.667 |
| 16 | 0.360 [0.333, 0.389] | 0.638 | 4.31 | 0.998 / 0.001 | 1.25 |
| 20 | 0.505 [0.483, 0.529] | 0.494 | 10.4 | 1.000 / 0.000 | 1.0 |
| 25 | 0.691 [0.676, 0.707] | 0.309 | 29.7 | 1.000 / 0.000 | 1.0 |
| 30 | 0.835 [0.826, 0.843] | 0.165 | 77.4 | 0.999 / 0.001 | 1.0 |

: Tablo 7.4. Kaynak: `results/transport/momentum/summary.csv` (noktalar `k*`, `k*_clean`); kℓ_tr için
`analysis.json["momentum"]`.

**Gözlem.** Daha hızlı paketler daha az saçılır: ⟨T⟩, k ile monoton artar.

**Sistematik hata.** k_x0 = 12'de CAP, temiz kılavuzda bile olasılığın yaklaşık %1.3–1.4'ünü sahte olarak geri döndürür
(V5; temiz R = 0.014). Bu, o noktada R'yi yaklaşık 0.01 yukarı kaydırır.

**Ayrıntı notu.** k_x0 = 20 noktası (n = 32) strength taramasındaki W = 100 noktasından (n = 48) farklı bir realizasyon
sayısına sahiptir: 0.505 ile 0.494. İki değer güven aralıkları içinde uyumludur.

### 7.4 Sistem uzunluğuna bağımlılık

Ayarlar: ζ = 0.05, k_x0 = 20, n = 24, t_end = 1 + 0.5 L_d.

| L_d | W = 100: ⟨T⟩ [%95 CI] | W = 200: ⟨T⟩ [%95 CI] |
|---|---|---|
| 0.25 | 0.827 [0.801, 0.852] | 0.602 [0.544, 0.656] |
| 0.5 | 0.694 [0.665, 0.722] | 0.414 [0.368, 0.460] |
| 1 | 0.518 [0.490, 0.548] | 0.249 [0.221, 0.278] |
| 2 | 0.320 [0.302, 0.338] | 0.104 [0.092, 0.116] |
| 3 | 0.232 [0.214, 0.250] | 0.056 [0.050, 0.064] |

: Tablo 7.5. Kaynak: `results/transport/length/summary.csv`.

| | W = 100 | W = 200 |
|---|---|---|
| Difüzif fit ℓ/(L+ℓ): ℓ, χ²_red | 1.00 ± 0.03, **3.6** | 0.24 ± 0.01, 16.1 |
| Üstel fit A e^{−L/ξ}: ξ, χ²_red | 2.01 ± 0.06, 6.8 | **1.15 ± 0.04, 5.7** |
| Born ℓ = (π/2)ℓ_tr | 0.81 | 0.20 |
| var(ln T) / ⟨ln T⟩ (L_d = 3) | 0.042 / −1.48 | 0.133 / −2.93 |

: Tablo 7.6. Kaynak: `results/transport/analysis.json["length"]` (alanlar `ohmic_l`, `ohmic_l_err`, `ohmic_chi2_red`,
`exp_xi`, `exp_xi_err`, `exp_chi2_red`, `born_l`, `lnT_var`, `lnT_mean`). Fitler ağırlıklı en küçük karelerdir
(`figures/CAPTIONS.md`).

**Gözlemler.**

- **W = 100.** Difüzif (Ohm) biçim daha iyi uyar. Uydurulan ℓ değeri Born değerine yakındır. Ölçülen aralıkta taşınım
  difüzif davranışla uyumludur.
- **W = 200.** ⟨T⟩ difüzif biçimden daha hızlı azalır ve üstel bir biçim (ξ ≈ 1.15) daha iyi uyar. ξ değeri, 6 kanal için
  beklenen yarı-1B lokalizasyon uzunluğuyla (ξ ~ Nℓ ~ 1.2, mertebe kestirimi [6]) aynı mertebededir.

**Yorum.** W = 200'deki davranış **lokalizasyonun başlangıcıyla tutarlıdır, ancak kanıt değildir.** Gerekçeler:

1. Yalnızca L_d/ξ ≤ 2.6'ya ulaşılmıştır.
2. var(ln T), |⟨ln T⟩|'nin çok altında kalır. Oysa derin lokalize rejimde ln T geniş dağılır. Paketin enerji bandı
   üzerinden ortalama alması bu dalgalanmaları bastırır.
3. İki fit de istatistiksel olarak tam değildir (χ²_red > 1). Ayrıca fitler, eşleştirilmiş tasarımın getirdiği
   korelasyonları hesaba katmaz.
4. W ≈ E₀'da bölgenin bir kısmı klasik olarak yasaktır (V > E₀); klasik engellenme de katkıda bulunur.
   - README §9.5 bu oranı "yaklaşık %15" olarak verir, ancak hesabın kaynağını belirtmez.
   - Bu değer, tek-nokta dağılımı Gauss varsayılırsa P(V > E₀) = ½ erfc(E₀/(√2 W)) kestirimiyle örtüşür: E₀ = 204.93 ve
     W = 200 için 0.153. Bu taslak için yalnızca bu formülle aritmetik olarak hesaplanmıştır.
   - Bu bir **analitik kestirimdir**, realizasyonlar üzerinde ölçülmüş bir değer değildir. Ayrıca tek-nokta istatistiği
     yalnızca yaklaşık olarak Gauss'tur (§5.5) ve bölge kenarlarındaki S(x) penceresi hesaba katılmamıştır. [REVIEW:
     README'nin bu değeri nasıl elde ettiği belgelenmelidir.]

**Kapsam.** Bu tarama n = 24 ile diğer taramalardan küçüktür ve tekrar üretilmemiştir (§13). Taslak "Anderson
lokalizasyonu gözlendi" ifadesini **kullanmaz**; yalnızca "üstel benzeri azalma" der.

### 7.5 Kapalı kutu sonuçları

Ayarlar: L × L kapalı kutu, L ∈ {1, 2}, düzensizlik her yerde, ζ = 0.05, k_x0 = 20, t_end = 0.5, n = 16 (W = 0 için tek
koşu). Ölçülen nicelik, katılım oranıdır (*participation ratio*): PR = (∫ρ)²/∫ρ².

| W | geç zaman PR/alan, L = 1 | geç zaman PR/alan, L = 2 | PR(L = 2)/PR(L = 1) |
|---|---|---|---|
| 0 | 0.375 | 0.366 | 3.9 |
| 100 | 0.415 ± 0.007 | 0.425 ± 0.022 | 4.1 |
| 200 | 0.334 ± 0.020 | 0.312 ± 0.031 | 3.7 |
| 300 | 0.267 ± 0.023 | 0.228 ± 0.032 | 3.4 |

: Tablo 7.7. ± = SE. Kaynak: `results/localization/summary.csv` (sütunlar `PR_late_over_area_mean`,
`PR_late_over_area_se`, `PR_late_mean`).

**Gözlem.** Alan dört katına çıktığında PR de yaklaşık dört katına çıkar (oran 3.4–4.1). Yani PR alanla yaklaşık
orantılıdır. Bu boyutlarda **2B lokalizasyon saptanmamıştır.**

**Yorum.**

- Bu beklenen bir **negatif sonuçtur.** Yerleşik kurama göre 2B lokalizasyon uzunluğu ξ ≈ ℓ e^{πkℓ/2}, W ≤ 150 için kutu
  boyutunu mertebelerce aşar [3]. Kutu boyutuna ancak Born yaklaşımının çöktüğü bölgede yaklaşır.
- W = 300'deki PR/alan düşüşünün klasik bir katkısı vardır; bu düşüş kuantum lokalizasyonuna atfedilemez.
  - README §9.6, V > E₀ olan alan oranını "yaklaşık %25" olarak verir.
  - Bu değer, aynı Gauss kestirimiyle (W = 300 için 0.247) örtüşür.
  - Ölçülmüş bir değer değildir. [REVIEW: kaynak belgelenmelidir.]
- **PR/alan'ın mutlak değeri** (≈ 0.23–0.43), tam yayılmış bir durum için bile 1 değildir. `observables.py` belgesine
  göre ergodik (rastgele dalga) bir durumun uzun-zaman ortalamalı yoğunluğu için PR ≈ A/2 ile A/3 arasındadır. Bu nedenle
  lokalizasyon testi mutlak değere değil, **L ile ölçeklemeye** dayanır.
- W = 100'de kutu, temiz kutudan **daha düzgün** doldurulur. Akla yatkın bir açıklama: düzensizlik momentum yönünü
  rastgeleleştirir, oysa temiz paket ağırlıklı olarak x boyunca hareket eder. Bu bir yorumdur (README §9.6).

### 7.6 Bariyer ve tünelleme incelemesi

Ayarlar: genişliği 0.1 olan, bölgenin ortasında yumuşak dikdörtgen bariyer; ζ = 0.05, L_d = 1, n = 24, t_end = 1.5.
Bariyer yüksekliği V_b, boylamasına enerji k_x0²/2 = 200'ün altında (150) veya üstündedir (250).

| V_b | W = 0 | W = 50 | W = 100 | W = 200 |
|---|---|---|---|---|
| 150 | 0.708 | 0.587 [0.561, 0.613] | 0.402 [0.365, 0.440] | 0.203 [0.174, 0.234] |
| 250 | 0.355 | 0.315 [0.296, 0.333] | 0.247 [0.218, 0.277] | 0.160 [0.132, 0.190] |

: Tablo 7.8. ⟨T⟩ [%95 CI]. Kaynak: `results/transport/barrier/summary.csv`.

**Gözlem.** Her iki bariyer yüksekliğinde de düzensizlik iletimi **yalnızca azaltmıştır.** Bu parametrelerde
**düzensizlik destekli tünelleme gözlenmemiştir.**

**Kapsam.** Bu, iki bariyer yüksekliği ve üç düzensizlik gücüyle sınırlı bir negatif sonuçtur.

### 7.7 Realizasyon yakınsaması ve istatistiksel analiz

Ayarlar: W = 100, ζ = 0.05, L_d = 1; 192 realizasyon.

**Sonuçlar** (`results/transport/convergence/summary.csv`):

- ⟨T⟩ = 0.509 [0.500, 0.519], SE = 0.0048.
- SE, 1/√n ile azalır: n = 48'de 0.010, n = 192'de 0.005.
- 48 realizasyonluk ortalama (0.494, strength taraması) tam ortalamayla 1.5 SE içinde uyumludur.
- 24–48 realizasyonlu taramalar, ⟨T⟩'de ≳ 0.03'lük farkları ayırt edebilir (README §9.8).

**İstatistik araçları** (`stats_utils.py`):

- ortalama, örneklem standart sapması, standart hata;
- Student-t ve percentile-bootstrap %95 güven aralıkları;
- medyan;
- realizasyonlar üzerindeki dağılım.

**İstatistiksel nitelendirmeler:**

1. Güven aralıkları nokta başınadır; eğri boyunca eşzamanlı değildir.
2. Eşleştirilmiş tasarım nedeniyle bir sweep'teki noktalar koreledir.
3. Bölüm 7.2'deki eşleştirilmiş kontrastlar dışında çoklu karşılaştırma düzeltmesi yapılmamıştır.
4. Önceden kayıtlı (doğrulayıcı) bir analiz yoktur.

### 7.8 Önerilen şekiller (Bölüm 7)

| Şekil | Dosya yolu | Gösterdiği | Veri kaynağı | Önerilen başlık |
|---|---|---|---|---|
| Ş7.1 | `figures/transport_vs_strength.png` | T, R ve P_in'in W'ya bağımlılığı; Born + difüzyon tahmini | `results/transport/strength/summary.csv` | "Düzensiz bölge (L_d = 1, ζ = 0.05) boyunca dalga paketi taşınımı. Ansambl ortalamaları ve %95 bootstrap aralıkları (n = 48). Kesikli çizgi: parametresiz Born + difüzyon tahmini, uydurma değildir." |
| Ş7.2 | `figures/born_ratio_vs_kl.png` | ⟨T⟩_sim/T_est'in kℓ_tr'ye bağımlılığı (üç tarama) | `results/transport/analysis.json` | "Simüle edilen ⟨T⟩'nin Born + difüzyon tahminine oranı. kℓ_tr ≳ 10 için oran 0.95–1.12; kℓ_tr ≈ 1–3'te 1.2–1.3. Mekanizma belirlenmemiştir." [REVIEW: mevcut `CAPTIONS.md` metni "≈ 1.05'e yaklaşır" der; bkz. §10.4 K-2.] |
| Ş7.3 | `figures/transport_vs_correlation_length.png` | T ve R'nin kζ'ye bağımlılığı; Born ℓ_tr | `results/transport/corr_length/summary.csv` | "Sabit W = 100'de korelasyon uzunluğunun etkisi (n = 48). ⟨T⟩ kζ ≈ 1 civarında minimumdur; minimum yalnızca 7 noktalı ızgara çözünürlüğünde belirlenmiştir." |
| Ş7.4 | `figures/paired_corr_length.png` | eşleştirilmiş bootstrap kontrastları; minimum konumu | `results/analysis/paired_corr_length.json` | "ζ taramasının tohum bazında eşleştirilmiş yeniden analizi (POST-HOC; 48 tohum; B = 20000)." |
| Ş7.5 | `figures/transport_vs_momentum.png` | T'nin k'ya bağımlılığı; temiz kılavuz değerleri | `results/transport/momentum/summary.csv` | "Gelen dalga sayısının etkisi (W = 100, ζ = 0.05, n = 32). k_x0 = 12'de CAP sahte geri dönüşü R'yi yaklaşık 0.01 yukarı kaydırır." |
| Ş7.6 | `figures/transport_vs_length.png` | ⟨T⟩'nin L_d'ye bağımlılığı; difüzif ve üstel fitler | `results/transport/length/summary.csv`, `analysis.json` | "Düzensiz bölge uzunluğuna bağımlılık (n = 24). W = 200'de üstel biçim daha iyi uyar; bu lokalizasyon kanıtı değildir (L_d/ξ ≤ 2.6, χ²_red > 1)." |
| Ş7.7 | `figures/localization_closed_box.png` | PR/alan'ın zamana bağımlılığı; geç zaman ortalamaları | `results/localization/` | "Kapalı kutularda katılım oranı (L = 1, 2; n = 16). PR alanla orantılıdır; bu boyutlarda 2B lokalizasyon saptanmamıştır." |
| Ş7.8 | `figures/transport_barrier.png` | bariyer + düzensizlik iletimi | `results/transport/barrier/summary.csv` | "Bariyer (V_b = 150, 250) ve düzensizlik birlikte (n = 24). Düzensizlik iletimi yalnızca azaltmıştır." |
| Ş7.9 | `figures/realization_convergence.png` | koşan ortalama; SE ∝ 1/√n | `results/transport/convergence/` | "Ansambl ortalamalarının realizasyon sayısıyla yakınsaması (W = 100, n = 192)." |
| Ş7.10 | `figures/disorder_examples.png` | üç ζ için V_dis örnekleri (tohum 0, W = 100) | `disorder.py` ile deterministik üretim; saklı veri dosyası yok, ancak tohumla tam olarak yeniden üretilebilir | "Aynı tohumla üç korelasyon uzunluğu için düzensizlik realizasyonları (illüstrasyon)." |
| Ş7.11 | `figures/waveguide_snapshots.png` | tek bir realizasyonda \|ψ\|²'nin dört anlık görüntüsü | **saklı veri yok**; çözücü `make_figures.py::fig_snapshots` içinde çalıştırılır | "**İllüstrasyon.** Tek bir realizasyon için (W = 150, ζ = 0.05, tohum 0) olasılık yoğunluğu. Nicel analizde kullanılmamıştır." |
| Ş6.1 | `figures/validation_solver.png` | V2–V4 doğrulama ve yakınsama | `results/validation/*.csv` | "DST split-step çözücünün doğrulanması: serbest Gauss'a karşı hata, Δt'de ikinci derece yakınsama, uzaysal yakınsama." |

> **Taslak notu (Bölüm 7).**
> - **Taramalar:** `results/transport/{strength,corr_length,momentum,length,barrier,convergence}/summary.csv` (sütunlar
>   `T_mean`, `T_ci_boot_lo`, `T_ci_boot_hi`, `R_mean`, `T_std`, `T_se`, `P_in_mean`, `lnT_mean`, `lnT_std`, `T_n`,
>   `t_end`).
> - **Analizler:** `results/transport/analysis.json`; `results/analysis/paired_corr_length.json`;
>   `results/localization/summary.csv`.
> - **Belgeler:** `README.md` §8–§9; `figures/CAPTIONS.md`; `make_figures.py`.

---

## 8. PINN değerlendirmesi

### 8.1 Eğitim düzeni

- **Problem:** 1 × 1 kutu, t ∈ [0, 0.02], Gauss paket (σ = 0.07, k_x0 = 15, E₀ = 112.5).
- **Referans:** 127² ızgara, Δt = 10⁻⁵ ile DST referansı. Hatası yaklaşık 10⁻⁶'dır (V3, V4).
- **Koşular:** 18 çalışma koşusu ve 3 teşhis koşusu.
- **Değerlendirme** (`evaluate.py`): bağıl L2(ψ) ve bağıl L2(|ψ|²), norm sapması, kütle merkezi hatası, P(x > 0.5)
  hatası, eğitim süresi ve referans çözücü süresi.
- **Donanım:** yalnızca CPU (GPU yok).

**Önemli kapsam notu.** PINN, taşınım deneylerinin dalga kılavuzu geometrisinde **çalıştırılmamıştır.** Yalnızca daha
kısa ve daha küçük bir kutu problemi değerlendirilmiştir. Bu nedenle PINN sonuçları ana taşınım sorusuna doğrudan katkı
sağlamaz; yöntemsel bir değerlendirmedir.

### 8.2 Klasik çözücüyle karşılaştırma

Tablo 8.1, t = T anındaki bağıl L2 hatalarını, aynı yapılandırmanın aynı potansiyelli referans çözümüne göre verir.

| Yapılandırma | koşu(lar) | bağıl L2 ψ(T) | bağıl L2 \|ψ\|²(T) | maks \|norm − 1\| |
|---|---|---|---|---|
| temiz, 3 ağ tohumu | `clean_fourier_uniform_s{0,1,2}` | 0.489, 0.487, 0.567 (ort. 0.514, sd 0.046) | 0.57–0.76 | ≤ 0.061 |
| W = 28 (E₀/4) | `W28_fourier_uniform_s0` | 0.530 | 0.607 | 0.038 |
| W = 56 (E₀/2), 3 ağ tohumu | `W56_fourier_uniform_s{0,1,2}` | 0.586, 0.574, 0.610 (ort. 0.590, sd 0.018) | 0.67–0.75 | ≤ 0.072 |
| W = 56, farklı düzensizlik realizasyonu | `W56_d1_fourier_uniform_s0` | 0.582 | 0.760 | 0.038 |
| W = 112 (E₀) | `W112_fourier_uniform_s0` | 0.713 | 0.957 | 0.053 |
| W = 56, ζ = 0.04 / 0.08 / 0.16 | `W56_zeta0.04…`, `W56_fourier_uniform_s0`, `W56_zeta0.16…` | 0.685 / 0.586 / 0.497 | 0.90 / 0.69 / 0.58 | ≤ 0.04 |
| Fourier özniteliği yok: temiz / W = 56 | `*_nofourier_uniform_s0` | **1.152** / 0.724 | 1.42 / 0.78 | ≤ 0.05 |
| RAD örnekleme: temiz / W = 56 | `*_fourier_rad_s0` | 0.469 / 0.525 | 0.45 / 0.56 | ≤ 0.025 |
| zaman penceresi T = 0.01 / 0.02 / 0.04 (temiz) | `clean_T0.01…`, `clean_fourier_uniform_s0`, `clean_T0.04…` | **0.184** / 0.489 / **0.986** | 0.15 / 0.57 / 1.02 | 0.009 / 0.023 / 0.123 |
| 3 kat eğitim bütçesi (temiz) | `clean_long_fourier_uniform_s0` | 0.281 | 0.322 | 0.007 |

: Tablo 8.1. Kaynak: `results/pinn/summary.csv` (sütunlar `rel_l2_final`, `rel_l2_rho_final`, `norm_max_dev`,
`train_time_s`, `reference_time_s`).

**Doğruluk:**

- **PINN, bu bütçede taşınım fiziği için yeterince doğru değildir.** İncelenen en kolay ayarda (temiz, T = 0.02) bile son
  hata yaklaşık %49'dur. Referans çözücünün hatası yaklaşık 10⁻⁶'dır.
- Tek bir tohumla yapılan T = 0.04 koşusunda son hata 0.986'dır. Bu, ψ = 0 tahmininden daha iyi değildir.

**Maliyet:**

- Referans çözücü 0.7–2.9 s sürmüştür (`reference_time_s`). PINN eğitimi aynı CPU'da 470–3100 s sürmüştür
  (`train_time_s`); standart bütçeli koşular yaklaşık 1000 s'dir.
- Bu, yaklaşık 10⁵ kat daha büyük bir hata için kabaca 500 kat daha fazla hesaplama anlamına gelir (README §9.10).

**Sonuç.** Bu ileri problem için PINN'in klasik çözücüye göre **bir üstünlüğü bulunmamıştır.**

### 8.3 Gözlenen başarısızlık modları

1. **Zaman penceresi uzunluğu.** Hata T = 0.01'de %18'dir, T = 0.04'te %99'a çıkar ve zaman içinde düzenli olarak büyür
   (`figures/pinn_error_vs_time.png`; `results/pinn/<run>/eval_timeseries.csv`). PINN, paket duvarlara ulaştığında ortaya
   çıkan girişim saçaklarını yakalayamaz (README §9.10).
2. **Düzensizlik.**
   - Hata W ile artar: 0.53 → 0.59 → 0.71.
   - Korelasyon uzunluğu kısaldıkça da artar: ζ küçüldükçe 0.50 → 0.59 → 0.69.
   - Üçer tohumluk temiz ile W = 56 farkı (0.514 ile 0.590) yalnızca bir eğilimdir ve %5 düzeyinde **anlamlı değildir**:
     Welch t = 2.70, p = 0.086, n = 3 + 3 [39]. Bu değer bu taslak için `results/pinn/summary.csv`'den yeniden hesaplanmış ve
     README ile aynı bulunmuştur.
3. **Spektral yanlılık.** Fourier öznitelikleri temiz problem için zorunludur: öznitelik olmadan hata 1.15'tir [22, 23].
4. **Kolokasyon noktalarına aşırı uyum.** Sabit noktalarda eğitim düşük eğitim kaybına rağmen 1.22 hata vermiştir (§6.7).
5. **Kısmi iyileştirmeler.**
   - RAD örnekleme ılımlı bir iyileşme sağlamıştır; bu sonuç tek tohumludur.
   - 3 kat bütçe hatayı 0.49'dan 0.28'e düşürmüştür. Yani PINN bu bütçede hesaplamayla hâlâ iyileşmektedir.
6. **Olasılık korunumu.** Norm yalnızca yumuşak bir kısıt olmasına karşın %2–7 içinde korunmuştur (T = 0.04 için %12).

**Önerilen şekiller:**

| Şekil | Dosya yolu | Veri kaynağı | Önerilen başlık |
|---|---|---|---|
| Ş8.1 | `figures/pinn_errors.png` | `results/pinn/summary.csv` | "Her PINN koşusunun t = T'deki bağıl L2 hatası, aynı yapılandırma ve aynı düzensizlik realizasyonuna sahip DST referansına göre." |
| Ş8.2 | `figures/pinn_error_vs_time.png` | `results/pinn/<run>/eval_timeseries.csv` | "PINN dalga fonksiyonunun bağıl L2 hatasının zamana bağımlılığı (temiz, düzensiz ζ = 0.08 ve iki katına çıkarılmış zaman penceresi)." |
| Ş8.3 | `figures/pinn_vs_reference_W56.png` | `results/pinn/W56_fourier_uniform_s0/` | "Düzensiz PINN problemi (W = 56, ζ = 0.08, düzensizlik tohumu 0) için t = T'de PINN ve DST referansı. Tek bir koşudur, temsili değildir." |

### 8.4 Hesaplama kısıtları

- **Bütçe ve donanım.** Yalnızca CPU kullanılmıştır; koşu başına eğitim yaklaşık 15–50 dakikadır.
- **Kapsam.** Tek bir mimari incelenmiştir ve yapılandırma başına en fazla 3 ağ tohumu vardır.
- **İş parçacığı bağımlılığı.** PINN eğitimi PyTorch iş parçacığı sayısına bağlıdır:
  - 2 iş parçacığıyla ağırlıklar bit düzeyinde yeniden üretilmiştir;
  - 4 iş parçacığıyla son hata 0.48850'den 0.48987'ye değişmiştir (`docs/REPRODUCIBILITY_AUDIT.md` §5).

  Bu fark tohumdan tohuma dağılımın çok altındadır, ancak iş parçacığı sayısı önbellek belirtimine eklenmiştir
  (`docs/PHASE2_METHOD_SELECTION.md` §1).
- **Değerlendirilmeyenler.** Daha büyük ağlar, daha uzun eğitim, GPU, zaman adımlı (*time-marching*) eğitim ve dalga
  kılavuzu geometrisi denenmemiştir. Sonuçlar **yalnızca bu bütçe ve bu standart teknikler** içindir; PINN'ler hakkında
  genel bir hüküm değildir [26, 29].

> **Taslak notu (Bölüm 8).**
> - **Sonuçlar:** `results/pinn/summary.csv`; `results/pinn/*/train_info.json`, `eval_summary.json`,
>   `eval_timeseries.csv`.
> - **Belgeler:** `README.md` §6, §9.10, §10; `docs/REPRODUCIBILITY_AUDIT.md` §5; `figures/CAPTIONS.md`.
> - **Yeniden hesaplama:** Welch testi, `scipy.stats.ttest_ind(equal_var=False)` ile CSV değerlerinden.

---

## 9. Bağımsız sayısal doğrulama hattı (RGF): kısa durum

> **Bu bölüm ana fiziksel araştırmadan ayrıdır.**
>
> - RGF çalışması, Bölüm 7'deki zamana bağlı sonuçları bağımsız bir sabit-enerji yöntemiyle doğrulamak için
>   başlatılmıştır ve **tamamlanmamıştır.**
> - Bölüm 7'deki hiçbir sayıyı değiştirmez: taşınım verileri DST çözücüsünden (dx = 1/64) gelir, RGF testleri ayrı bir
>   durağan çözücüye aittir.
> - Ayrıntılı ölçümler, W/Z teşhis çalışmaları ve uygulama ayrıntıları **Ek A**'dadır.

**Varlık nedeni.**

- §6.6'daki iç doğrulama ve yakınsama testleri, çözücünün kendi ayrıklaştırma ve uygulama hatalarını sınar; kod
  paylaşmayan bağımsız bir yöntemle düzensiz taşınım karşılaştırması içermez.
- Bu boşluk için mod uzayında özyinelemeli Green fonksiyonu (RGF) ve Fisher–Lee bağıntısı seçilmiş, önceden
  belirlenmiş B1–B4 ölçütleri tanımlanmıştır (`docs/PHASE2_METHOD_SELECTION.md` §2–§3) [32–34, 36].
- Belgenin kuralına göre taşınım sonuçları **yalnızca B1–B4'ün hepsi geçerse** bağımsız olarak doğrulanmış sayılır
  (`docs/RGF_B2_VALIDATION.md`).

| Test | Kısa tanım | Durum |
|---|---|---|
| B1 | temiz kılavuz, 6 konfigürasyon | A–E geçti; **B1-F FAIL** (faz kriteri 10⁻¹⁰, değiştirilmedi) |
| B2 | y'den bağımsız bariyer; RGF ve TDSE, bağımsız 1B ODE'ye karşı | **PASS** |
| B3 | sabit düzensizlik realizasyonları; RGF ile seyrek Green referansının karşılaştırması | **FAIL** (dx kriteri) |
| B3X | 6 konfigürasyon; bağımsız dalga fonksiyonu referansı | **FAIL** (dx kriteri, C11) |
| B3Y | başarısız B3X noktalarında ince ızgara, 4 yöntem | **PASS** |
| B4 | sabit düzensizlik realizasyonunda TDSE ile durağan çözümün karşılaştırması | **YAPILMADI** |

: Tablo 9.1. Kaynaklar ve ölçülen değerler: Ek A.1.

**Sonuç.** Bağımsız durağan doğrulama **tamamlanmamıştır.** Düzensizlik içeren taşınım sonuçları bağımsız bir
sabit-enerji yöntemiyle henüz doğrulanmamıştır; bunun için B4 gereklidir.

**W/Z teşhisi (bir cümleyle).** B1-F başarısızlığını anlamak için RGF'nin cebirsel olarak eşdeğer bir yeniden yazımı
(W/Z) sayısal bir teşhis çalışması olarak incelenmiştir. Sonuçlar:

| Konu | Durum |
|---|---|
| W/Z | **NOT VALIDATED** |
| coupled-block deneyi | **DESIGN INSUFFICIENT** (W/Z FAIL anlamına gelmez) |
| ana coupling-scaling deneyi | **NOT RUN** (yalnızca manifestler üretilmiştir) |
| uygulama düzeyi R0 | **PENDING** |

Bu çalışma fiziksel bir sonuç üretmez (Ek A.2).

> **Taslak notu (Bölüm 9).** `docs/PHASE2_METHOD_SELECTION.md` §2–§3; `docs/RGF_B1/B2/B3/B3X/B3Y_VALIDATION.md`;
> `results/rgf_b{1,2,3,3x,3y}_validation/summary.json`; `experiments/wz_rgf_diagnostic/`; `REPORT_READINESS_AUDIT.md` §4.

---

## 10. Tartışma

### 10.1 Ana bulguların yorumu

**Korelasyon uzunluğu.** ⟨T⟩'nin ζ'ye göre monoton olmaması ve kζ ≈ 1 civarındaki minimumu, Born kuramının taşınım ortalama
serbest yolu ℓ_tr için öngördüğü monoton olmayan davranışla niteliksel olarak uyumludur [9].

- **Belirlenen:** minimumun varlığı. Keşifsel düzeyde ve post-hoc eşleştirilmiş kontrastlarla desteklenmiştir.
- **Belirlenmeyen:** minimumun tam konumu (kζ ≈ 0.6–1.4) ve bu davranışın başka W, L_d veya kanal sayılarında sürüp
  sürmediği.

**Born tahminiyle karşılaştırma.**

- Zayıf saçılmada (kℓ_tr ≳ 10) tahmin, ölçülen ⟨T⟩'yi yaklaşık %12 içinde yakalar.
- Güçlü saçılmada simülasyon tahminden sistematik olarak yüksek iletim verir. Bu sapma üç taramada ortak bir kℓ_tr eğilimi
  izler.
- Bu gözlem **bu geometriye özgü bir sonuçtur.** Nedeni belirlenmemiştir. Yalnızca eğilim yorumlanmalıdır, çünkü π/2 öneki
  konvansiyona bağlıdır.

**Lokalizasyon.**

- Kapalı 2B kutularda lokalizasyonun saptanmaması, yerleşik kuramın öngördüğü negatif sonuçtur.
- Yarı-1B dalga kılavuzunda W = 200'deki üstel benzeri azalma, lokalizasyonun başlangıcıyla **tutarlıdır.** Ancak erişilen
  L_d/ξ, ln T dağılımı ve klasik yasak bölge katkısı nedeniyle **kanıt değildir.**

**PINN.** Mevcut bütçede klasik çözücüye karşı bir üstünlük bulunmamıştır. Başarısızlık modları literatürde bildirilen PINN
güçlükleriyle niteliksel olarak benzerdir: spektral yanlılık, uzun zaman pencereleri ve kolokasyona aşırı uyum [22, 26].
Bu taslak, bu benzerliğin nedenselliğini ayrıca test etmemiştir.

### 10.2 Gözlem ve yorum ayrımı

| Gözlem (ölçülmüş) | Yorum (test edilmemiş veya kuramdan alınmış) |
|---|---|
| ⟨T⟩ W ile azalır | daha güçlü saçılma; Born ℓ_tr küçülür (kuram) |
| ⟨T⟩(ζ) minimumu kζ ≈ 1'de | ℓ_tr'nin minimumu (kuram) |
| sd(T) ζ ile artar | daha az bağımsız korelasyon alanı ve daha zayıf öz-ortalama (yorum) |
| Born oranı kℓ küçüldükçe artar | mekanizma belirsiz (aday nedenler listelendi) |
| W = 200'de üstel benzeri L_d bağımlılığı | lokalizasyonun başlangıcıyla tutarlı; kanıt değil |
| kapalı kutuda PR ∝ alan | bu boyutlarda 2B lokalizasyon yok (kuramla uyumlu) |
| PINN hatası T ve W ile artar | girişim saçakları; spektral yanlılık (yorum) |

### 10.3 Keşifsel ve doğrulayıcı analiz

**Bu çalışmadaki bütün fiziksel analizler keşifseldir.**

- ζ için eşleştirilmiş kontrastlar sonuçlar görüldükten sonra seçilmiştir (post-hoc).
- Önceden kaydedilmiş doğrulayıcı bir fiziksel kampanya tasarlanmıştır (`docs/PHYSICAL_CAMPAIGN_PLAN.md`, revizyon 3):
  - 68 testlik Holm ailesi;
  - iki aşamalı bootstrap;
  - kapalı-başarısızlık (*fail-closed*) mantığı.
- Ancak bu kampanya **dondurulmamış, hash'lenmemiş ve çalıştırılmamıştır.** Gerekli kök entropi ve C0 tohumu henüz
  kaydedilmemiştir (`docs/CAMPAIGN_STATUS_CONTRACT.md`).

### 10.4 Taslak sırasında bulunan tutarsızlıklar [REVIEW]

Bunlar uzlaştırılmadan raporlanmıştır.

| # | Tutarsızlık | Kaynaklar |
|---|---|---|
| K-1 | README'nin depo yapısı satırı (satır 596) "B3 and B3X: dx criterion fails; B4 open" der ve B3Y'den (PASS) söz etmez; `docs/` listesinde B3Y yoktur | `README.md` §12; `docs/RGF_B3Y_VALIDATION.md`; `results/rgf_b3y_validation/summary.json` |
| K-2 | Born oranı grafiğinin açıklaması "≈ 1.05'e yaklaşır" der; README "kℓ_tr ≳ 10 için 0.95–1.12" der | `figures/CAPTIONS.md`; `README.md` §9.2; `analysis.json` |
| K-3 | Dondurulmuş v2 spec'inin başlığı "nothing implemented or run" der; oysa altyapı ve manifestler mevcuttur | `PREDICTIONS_coupling_scaling_v2.md`; `fac5f5e` |
| K-4 | `docs/rgf_b3y_protocol.json` içindeki tarihsel `status` metni, kampanyanın çalıştırılmadığını söyler | `docs/RGF_B3Y_VALIDATION.md` (çelişkiyi açıklar) |
| K-6 | Paired analizin girdisi `results/reproduction/...` dizinidir | `paired_corr_length.json`; `REPRODUCIBILITY_AUDIT.md` §5 (bit düzeyinde aynı) |
| K-7 | README §6, sabit kolokasyonlu PINN koşusu için "20000 taze noktada kalıntı 0.22" der. R0.1'de depo, git geçmişi (değer `8f013cf` ile README'ye girmiştir) ve PINN kayıtları yeniden tarandı; destekleyen bir dosya bulunamadı (`train_info.json`'da `holdout_pde` alanı yok). **İddia bu taslaktan çıkarıldı.** README'de kalmaya devam eder | `README.md` §6; `results/pinn/diag_fixed_points_lbfgs/` |
| K-8 | Araştırma sorusunun kapsamı: README §1 W ve ζ'yi birlikte, Phase 2 belgesi yalnızca ζ'yi sorar. R0.1'de başlıkla uyumlu olarak ζ ana soru alındı (§4.3); çelişki değil, kapsam seçimi | `README.md` §1; `docs/PHASE2_METHOD_SELECTION.md` |
| K-9 | "Yüksek boyutlu uzaylarda ML" motivasyonu depoda belgelenmemiştir. Belgelenen motivasyonlar: "ML ile kuantum kimyası / moleküler biyoloji simülasyonlarını hızlandırmak" (`5f56539`) ve ENAQT (`4fe3eb2`, `42cea06`). R0.1'de öğrencinin anlatımı olarak ayrıldı (§4.5) | `git log`; `docs/AUDIT.md` 2.4 |
| K-10 | README'deki "alanın ≈ %15 / %25'i V > E₀" değerleri bir sonuç dosyasına bağlı değildir. Değerler Gauss tek-nokta kestirimiyle (0.153 / 0.247) örtüşür; ölçülmüş değildir (§7.4, §7.5) | `README.md` §9.5, §9.6 |

> **Taslak notu (Bölüm 10).** `README.md` §9.2–§9.11; `REPORT_READINESS_AUDIT.md` §7; `docs/PHYSICAL_CAMPAIGN_PLAN.md`;
> `docs/CAMPAIGN_STATUS_CONTRACT.md`; `results/pinn/diag_fixed_points_lbfgs/train_info.json`.

---

## 11. Sınırlamalar

1. **Geometri ve boyut.**
   - Dalga kılavuzunun yalnızca 6 açık kanalı vardır ve düzensiz bölge en fazla L_d = 3'tür.
   - Kapalı kutular en fazla 2 × 2'dir.
   - Bu boyutlar savunulabilir bir 2B Anderson lokalizasyonu çalışmasına izin vermez; yarı-1B lokalizasyon için yalnızca
     bir ipucu verir.
2. **Enerji ortalaması.** Gauss paketi k ≈ 20 etrafında σ_k ≈ 3.5 genişliğinde bir bant örnekler. T, R ve
   dalgalanmaları enerji üzerinden ortalanmıştır.
3. **Sonlu gözlem süresi.** t_end'de P_in en fazla %1.15'tir (en kötü durum W = 200, L_d = 3:
   `length/summary.csv`, `P_in_mean`). Hapsolmuş yavaş bileşenler "hâlâ içeride" sayılır.
4. **Soğurucu sınırlar.** CAP'ın geri döndürdüğü olasılık k_x0 ≥ 20 için ≤ 5.2×10⁻⁴ (k_x0 = 30'da 5.19×10⁻⁴;
   README "≤ 5×10⁻⁴" der), k_x0 = 16 için 1.3×10⁻³, k_x0 = 12 için 1.3×10⁻²'dir (`v5_cap_reflection.csv`). En düşük momentumda sistematik bir yanlılık vardır.
5. **Kuram karşılaştırması.**
   - Born oranları sınırsız 2B içindir. ℓ/(L+ℓ) bir mertebe tahminidir.
   - İkisi de dalga kılavuzu kanal yapısını, zayıf lokalizasyonu veya ekstrapolasyon uzunluklarını içermez.
   - Uyum ve sapmalar ampirik gözlemlerdir.
6. **İstatistik.**
   - Nokta başına 24–48 realizasyon vardır; yakınsama testi için 192.
   - Eşleştirilmiş tasarım bir tarama içindeki noktaları korele kılar.
   - Güven aralıkları nokta başınadır.
   - Analizlerin tümü keşifseldir; ζ kontrastları post-hoc'tur.
7. **Düzensizlik modeli.** Yalnızca Gauss-korelasyonlu düzensizlik incelenmiştir; speckle veya korelasyonsuz kafes
   düzensizliği yoktur. Tek-nokta istatistiği yaklaşık olarak Gauss'tur.
8. **Bağımsız doğrulama eksikliği.** Düzensiz taşınım için bağımsız durağan çapraz kontrol (B4) yapılmamıştır. Bağımsız
   bir yöntemle yapılan tek TDSE karşılaştırması (B2) düzensizlik içermez. RGF hattında B1-F ile B3 ve B3X başarısızdır
   (§9, Ek A).
9. **Tekrar üretilebilirlik.**
   - Bit düzeyinde yeniden üretilenler: strength, ζ ve kapalı kutu taramaları ile çözücü doğrulaması.
   - Yeniden çalıştırılmayanlar: momentum, length, barrier ve convergence taramaları ile 17 PINN koşusu.
   - 44 transport `.npz` dosyasının 22'sinde eski meta veri bulunur (`config_json` T = 1.0; gerçek t_end 1.25–2.5). Sayılar
     etkilenmez (`docs/REPRODUCIBILITY_AUDIT.md` §2).
10. **PINN.**
    - Yalnızca CPU bütçesiyle, tek mimariyle ve yapılandırma başına en fazla 3 tohumla çalışılmıştır.
    - Yalnızca kutu problemi incelenmiştir; dalga kılavuzu incelenmemiştir.
    - Sonuç iş parçacığı sayısına bağlıdır.
11. **Modellenmeyenler:**
    - etkileşimler;
    - zamana bağlı düzensizlik;
    - faz bozulması / açık sistem etkileri;
    - üç boyut.

> **Taslak notu (Bölüm 11).** `README.md` §10; `results/transport/length/summary.csv` (`P_in_mean`);
> `results/validation/v5_cap_reflection.csv`; `docs/REPRODUCIBILITY_AUDIT.md` §2, §5; `REPORT_READINESS_AUDIT.md` §7.1.

---

## 12. Sonuç ve gelecek çalışma

### 12.1 Sonuç

Bu çalışma, 6 kanallı 2B sert duvarlı bir dalga kılavuzunda Gauss-korelasyonlu statik düzensizlik boyunca dalga paketi
iletimini, doğrulanmış bir DST split-step çözücüyle keşifsel olarak haritalamıştır.

**Ana bulgular (bu geometri ve bu parametreler için):**

1. ⟨T⟩, düzensizlik gücüyle monoton azalır.
2. ⟨T⟩, korelasyon uzunluğuna göre monoton değildir ve kζ ≈ 1 civarında bir minimum gösterir. Kanıt keşifsel ve post-hoc
   düzeydedir; minimumun konumu kaba çözünürlüktedir.
3. Born + difüzyon tahmini zayıf saçılmada yaklaşık %12 içinde uyar. Güçlü saçılmada ise kℓ'ye bağlı, mekanizması
   belirlenmemiş sistematik bir sapma vardır.
4. W = 200'de L_d'ye bağımlılık üstel benzeridir. Bu lokalizasyonla tutarlıdır, ancak kanıt değildir.
5. Kapalı kutularda bu boyutlarda 2B lokalizasyon saptanmamıştır.
6. Bariyer içeren düzensizlikte tünelleme desteği gözlenmemiştir.
7. PINN, mevcut bütçede taşınım fiziği için yeterli doğruluğa ulaşmamıştır.

**Doğrulama durumu.** Taşınım sonuçlarının bağımsız durağan doğrulaması **tamamlanmamıştır:** B4 yapılmamıştır. RGF/W-Z
teşhis çalışmaları ayrı bir sayısal güvenilirlik hattıdır:

| Konu | Durum |
|---|---|
| B1-F | FAIL |
| coupled-block | DESIGN INSUFFICIENT |
| W/Z | NOT VALIDATED |
| coupling-scaling v2 | NOT RUN |
| R0 | PENDING |

**Bu çalışma yeni bir kuram, yeni bir PINN mimarisi, Anderson lokalizasyonunun kanıtı veya deneysel doğrulama iddia
etmez.**

### 12.2 Gelecek çalışma

Bu çalışmalar henüz yapılmamıştır (README §9.11; `docs/PHYSICAL_CAMPAIGN_PLAN.md`).

1. **B4.** Seçilmiş noktalarda zamana bağlı çözüm ile durağan çözümün (RGF) enerji genişliği hesaba katılarak
   karşılaştırılması. Taşınım sonuçlarının bağımsız doğrulaması için gereklidir.
2. **Önceden kayıtlı ζ kampanyası.** Fiziksel doğrulayıcı kampanyanın dondurulması (kök entropi ve C0 tohumunun
   kaydedilmesi) ve çalıştırılması.
3. **Daha uzun bölgeler ve dar paketler.** W = 150–200'de L_d ≈ 6–8'e kadar ve ≥ 64 realizasyonla uzunluk taraması;
   sabit enerjiye yaklaşmak ve log-normal ln T istatistiğini sınamak için daha dar paketler (σ_x ≥ 0.5).
4. **Kanal sayısı.** ξ ∝ N'yi sınamak için dalga kılavuzu genişliği taraması (N = 3, 6, 12 kanal).
5. **Daha büyük kutular.** W ≥ 200'de L = 4 kapalı kutular ve klasik perkolasyon karşılaştırması.
6. **PINN.** Yapılandırma başına ≥ 5 tohum, GPU'da daha uzun bütçeler, zaman adımlı eğitim ve dalga kılavuzu
   geometrisinde PINN.
7. **W/Z.** Bağımsız uygulama düzeyi R0 incelemesi ile G-1 ve A-1'in çözümü; ancak bunlardan sonra coupling-scaling v2
   deneyi.
8. **Yüksek boyut.** Özgün motivasyonla ilgili yüksek boyutlu genelleme, bu çalışmanın kapsamı dışındadır. [TO BE
   COMPLETED: yazarlar ilgili hedefi ve yöntemi tanımlamalıdır.]

> **Taslak notu (Bölüm 12).** `README.md` §9.11; `docs/PHYSICAL_CAMPAIGN_PLAN.md`; `docs/CAMPAIGN_STATUS_CONTRACT.md`;
> `REPORT_READINESS_AUDIT.md` §8, §10.

---

## 13. Tekrar üretilebilirlik ve veri erişimi

**Depo.** `https://github.com/mrthnkzlrmk865-cmd/Quantum-Schrodinger-PINN`

- Dal: `claude/quantum-schrodinger-audit-iy6kn4`.
- Revizyon bilgisi:

  | belge | commit |
  |---|---|
  | rapor hazırlık denetimi | `e7aa728` |
  | R0 taslağı | `11b342c` |
  | R0.1 taslağı | bu dosyayı ekleyen commit (teslim notunda verilir) |
  | `origin/main` | `5bee7eb` (dosya ağacı `11b342c` ile aynı) |
- Fiziksel şekiller `make_figures.py` ile `1c72630` commit'inde üretilmiştir.

**Ortam.** Kesin sürümler `requirements-lock.txt`'tedir: Python 3.13, NumPy 2.5.3, SciPy 1.18.1, PyTorch 2.15.0 (CPU).
Zamanlamalar GPU'suz 4 çekirdekli bir CPU içindir (README §11).

**Komutlar** (README §11; bu taslak için **çalıştırılmamıştır**):

```bash
python -m pytest -q tests
python validate_solver.py
python experiment_disorder.py --sweep all
python experiment_localization.py --n-real 16
python experiment_pinn.py
python make_figures.py
python paired_analysis.py
```

**Köken bilgisi ve önbellek** (`provenance.py`):

- Her sonuç şunları kaydeder: tam yapılandırma, tohum listesi, realizasyon sayısı, yapılandırma parmak izi (SHA-256),
  kaynak kod parmak izi, git commit'i, UTC zamanı ve kütüphane sürümleri.
- Önbellek yalnızca parmak izi tam eşleştiğinde yeniden kullanılır.
- **Rapora girecek her sonuç için `--strict` önerilir** (README §11).

**Tekrar üretim durumu** (`docs/REPRODUCIBILITY_AUDIT.md` §5):

| Sonuç | Durum |
|---|---|
| V1–V7, strength, ζ, kapalı kutu | bit düzeyinde aynı |
| 2 iş parçacığıyla PINN `clean_fourier_uniform_s0` | ağırlıklar bit düzeyinde aynı |
| momentum, length, barrier, convergence; diğer PINN koşuları | yeniden çalıştırılmadı (**UNKNOWN**) [TO BE COMPLETED: `--strict` ile yeniden üretim] |

**Veri konumları:**

| Veri | Konum |
|---|---|
| taşınım (ham `.npz` ve `summary.csv`) | `results/transport/` |
| analiz | `results/analysis/` |
| kapalı kutu | `results/localization/` |
| PINN | `results/pinn/` |
| doğrulama | `results/validation/` |
| yeniden üretim | `results/reproduction/` |
| RGF | `results/rgf_b*_validation/` |
| W/Z | `experiments/wz_rgf_diagnostic/` |

> **Taslak notu (Bölüm 13).** `README.md` §11–§12; `docs/REPRODUCIBILITY_AUDIT.md`; `requirements-lock.txt`;
> `provenance.py`.

---

## 14. Kaynaklar

[1]–[31] `docs/LITERATURE.md` listesinden, [32]–[36] `docs/PHASE2_METHOD_SELECTION.md` §8'den aynen alınmıştır. Bu
belgeler bibliyografik verilerin kontrol edildiğini belirtir.

1. P. W. Anderson, "Absence of diffusion in certain random lattices," *Phys. Rev.* **109**, 1492 (1958). doi:10.1103/PhysRev.109.1492
2. E. Abrahams, P. W. Anderson, D. C. Licciardello, T. V. Ramakrishnan, "Scaling theory of localization: absence of quantum diffusion in two dimensions," *Phys. Rev. Lett.* **42**, 673 (1979). doi:10.1103/PhysRevLett.42.673
3. P. A. Lee, T. V. Ramakrishnan, "Disordered electronic systems," *Rev. Mod. Phys.* **57**, 287 (1985). doi:10.1103/RevModPhys.57.287
4. B. Kramer, A. MacKinnon, "Localization: theory and experiment," *Rep. Prog. Phys.* **56**, 1469 (1993). doi:10.1088/0034-4885/56/12/001
5. F. Evers, A. D. Mirlin, "Anderson transitions," *Rev. Mod. Phys.* **80**, 1355 (2008). doi:10.1103/RevModPhys.80.1355
6. C. W. J. Beenakker, "Random-matrix theory of quantum transport," *Rev. Mod. Phys.* **69**, 731 (1997). doi:10.1103/RevModPhys.69.731
7. E. Akkermans, G. Montambaux, *Mesoscopic Physics of Electrons and Photons*, Cambridge University Press (2007).
8. A. Lagendijk, B. van Tiggelen, D. S. Wiersma, "Fifty years of Anderson localization," *Physics Today* **62**(8), 24 (2009). doi:10.1063/1.3206091
9. R. C. Kuhn, O. Sigwarth, C. Miniatura, D. Delande, C. A. Müller, "Coherent matter wave transport in speckle potentials," *New J. Phys.* **9**, 161 (2007). doi:10.1088/1367-2630/9/6/161
10. J. Billy et al., "Direct observation of Anderson localization of matter waves in a controlled disorder," *Nature* **453**, 891 (2008). doi:10.1038/nature07000
11. G. Roati et al., "Anderson localization of a non-interacting Bose–Einstein condensate," *Nature* **453**, 895 (2008). doi:10.1038/nature07071
12. T. Schwartz, G. Bartal, S. Fishman, M. Segev, "Transport and Anderson localization in disordered two-dimensional photonic lattices," *Nature* **446**, 52 (2007). doi:10.1038/nature05623
13. L. Sanchez-Palencia, M. Lewenstein, "Disordered quantum gases under control," *Nature Physics* **6**, 87 (2010). doi:10.1038/nphys1507
14. M. Shinozuka, G. Deodatis, "Simulation of stochastic processes by spectral representation," *Appl. Mech. Rev.* **44**(4), 191 (1991). doi:10.1115/1.3119501
15. M. D. Feit, J. A. Fleck Jr., A. Steiger, "Solution of the Schrödinger equation by a spectral method," *J. Comput. Phys.* **47**, 412 (1982). doi:10.1016/0021-9991(82)90091-2
16. D. Kosloff, R. Kosloff, "A Fourier method solution for the time dependent Schrödinger equation as a tool in molecular dynamics," *J. Comput. Phys.* **52**, 35 (1983). doi:10.1016/0021-9991(83)90015-3
17. G. Strang, "On the construction and comparison of difference schemes," *SIAM J. Numer. Anal.* **5**, 506 (1968). doi:10.1137/0705041
18. T. Jahnke, C. Lubich, "Error bounds for exponential operator splittings," *BIT Numer. Math.* **40**, 735 (2000). doi:10.1023/A:1022396519656
19. J. G. Muga, J. P. Palao, B. Navarro, I. L. Egusquiza, "Complex absorbing potentials," *Phys. Rep.* **395**, 357 (2004). doi:10.1016/j.physrep.2004.03.002
20. M. Raissi, P. Perdikaris, G. E. Karniadakis, "Physics-informed neural networks: A deep learning framework for solving forward and inverse problems involving nonlinear partial differential equations," *J. Comput. Phys.* **378**, 686 (2019). doi:10.1016/j.jcp.2018.10.045
21. G. E. Karniadakis, I. G. Kevrekidis, L. Lu, P. Perdikaris, S. Wang, L. Yang, "Physics-informed machine learning," *Nat. Rev. Phys.* **3**, 422 (2021). doi:10.1038/s42254-021-00314-5
22. N. Rahaman et al., "On the spectral bias of neural networks," *ICML* (2019). arXiv:1806.08734
23. M. Tancik et al., "Fourier features let networks learn high frequency functions in low dimensional domains," *NeurIPS* 33 (2020). arXiv:2006.10739
24. S. Wang, H. Wang, P. Perdikaris, "On the eigenvector bias of Fourier feature networks: From regression to solving multi-scale PDEs with physics-informed neural networks," *Comput. Methods Appl. Mech. Eng.* **384**, 113938 (2021). doi:10.1016/j.cma.2021.113938
25. S. Wang, Y. Teng, P. Perdikaris, "Understanding and mitigating gradient flow pathologies in physics-informed neural networks," *SIAM J. Sci. Comput.* **43**, A3055 (2021). doi:10.1137/20M1318043
26. A. S. Krishnapriyan, A. Gholami, S. Zhe, R. M. Kirby, M. W. Mahoney, "Characterizing possible failure modes in physics-informed neural networks," *NeurIPS* 34 (2021). arXiv:2109.01050
27. S. Wang, S. Sankaran, P. Perdikaris, "Respecting causality for training physics-informed neural networks," *Comput. Methods Appl. Mech. Eng.* **421**, 116813 (2024). doi:10.1016/j.cma.2024.116813
28. C. Wu, M. Zhu, Q. Tan, Y. Kartha, L. Lu, "A comprehensive study of non-adaptive and residual-based adaptive sampling for physics-informed neural networks," *Comput. Methods Appl. Mech. Eng.* **403**, 115671 (2023). doi:10.1016/j.cma.2022.115671
29. K. Shah, P. Stiller, N. Hoffmann, A. Cangi, "Physics-informed neural networks as solvers for the time-dependent Schrödinger equation," NeurIPS 2022 Workshop on Machine Learning and the Physical Sciences. arXiv:2210.12522
30. D. C. Liu, J. Nocedal, "On the limited memory BFGS method for large scale optimization," *Math. Program.* **45**, 503 (1989). doi:10.1007/BF01589116
31. B. Efron, R. J. Tibshirani, *An Introduction to the Bootstrap*, Chapman & Hall (1993).
32. P. A. Lee, D. S. Fisher, "Anderson localization in two dimensions," *Phys. Rev. Lett.* **47**, 882 (1981). doi:10.1103/PhysRevLett.47.882
33. D. S. Fisher, P. A. Lee, "Relation between conductivity and transmission matrix," *Phys. Rev. B* **23**, 6851 (1981). doi:10.1103/PhysRevB.23.6851
34. A. MacKinnon, "The calculation of transport properties and density of states of disordered solids," *Z. Phys. B* **59**, 385 (1985). doi:10.1007/BF01328846
35. A. C. Davison, D. V. Hinkley, *Bootstrap Methods and their Application*, Cambridge University Press (1997).
36. S. Datta, *Electronic Transport in Mesoscopic Systems*, Cambridge University Press (1995).
37. D. P. Kingma, J. Ba, "Adam: A method for stochastic optimization," 3rd International Conference on Learning Representations (ICLR), San Diego (2015). arXiv:1412.6980
38. S. Holm, "A simple sequentially rejective multiple test procedure," *Scand. J. Statist.* **6**(2), 65–70 (1979). JSTOR 4615733
39. B. L. Welch, "The generalization of 'Student's' problem when several different population variances are involved," *Biometrika* **34**(1/2), 28–35 (1947). [REFERENCE DETAILS TO VERIFY: iki farklı DOI adayı bulundu (10.1093/biomet/34.1-2.28 ve 10.2307/2332510); doğru olan doi.org üzerinden doğrulanmadan DOI verilmemiştir.]
40. B. V. North, D. Curtis, P. C. Sham, "A note on the calculation of empirical P values from Monte Carlo procedures," *Am. J. Hum. Genet.* **71**, 439 (2002). PMID 12111669
41. B. Phipson, G. K. Smyth, "Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn," *Stat. Appl. Genet. Mol. Biol.* **9**(1), Article 39 (2010). doi:10.2202/1544-6115.1585

**Kaynakların doğrulanma yolları:**

- [1]–[31]: `docs/LITERATURE.md` "References" listesinden aynen alınmıştır.
- [32]–[36], [40], [41]: `docs/PHASE2_METHOD_SELECTION.md` §8 ve §8.1'den alınmıştır. [40] ve [41], kullanılan (b+1)/(B+1)
  bootstrap p-değeri düzeltmesinin kaynaklarıdır.
- [37]–[39]: R0.1 sırasında web aramasıyla yayın kayıtları üzerinden doğrulanmıştır.
  - [37]: arXiv kaydı.
  - [38]: JSTOR 4615733.
  - [39]: cilt ve sayfa doğrulanmış, DOI doğrulanmamıştır.
- Kullanım yerleri: Adam [37] §6.7'de, Holm düzeltmesi [38] §7.2'de, Welch testi [39] §8.3'te.

> **Taslak notu (Bölüm 14).** `docs/LITERATURE.md` "References"; `docs/PHASE2_METHOD_SELECTION.md` §8, §8.1. Yeni dış
> kaynaklar [37]–[39], doğrulama yollarıyla birlikte yukarıda belirtilmiştir.

---

## Ek A. RGF doğrulama hattı ve W/Z teşhis çalışmalarının ayrıntıları

> Bu ek, ana metinden ayrılan teknik ayrıntıları içerir. Buradaki sonuçların hiçbiri Bölüm 7'deki fiziksel bulguları
> değiştirmez.

### A.1 B1–B4 ölçümleri

| Test | Amaç | Durum | Ölçülen | Kaynak |
|---|---|---|---|---|
| B1 | temiz kılavuz, 6 konfigürasyon, 50 enerji; analitik kafes cevabı | A–E geçti; **B1-F FAIL** (faz kriteri 10⁻¹⁰, değiştirilmedi) | F (N_y = 63, dx = 1/1024, L_d = 3, N = 3483): faz büyüklüğü 1.76851694417934×10⁻¹⁰; 100 basamaklı kafes fazına göre 1.5424928001550615×10⁻¹⁰ | `results/rgf_b1_validation/summary.json`, `phase_audit_hp_summary.json` (alanlar `b1_phase_quantity_max`, `hp_max_abs_phase_err`); `docs/RGF_B1_VALIDATION.md` |
| B2 | y'den bağımsız bariyer; RGF ile bağımsız 1B ODE'nin karşılaştırması | **PASS** (tüm bağlayıcı kriterler) | dx = 1/512'de maks \|ΔT\| = 3.59×10⁻⁵ (V_b = 150), 5.95×10⁻⁵ (V_b = 250); kriter 2×10⁻⁴; TD ile enerji-ortalamalı ODE uzaklığı 3.2×10⁻⁴ / 7.0×10⁻⁴, kriter 10⁻³ | `results/rgf_b2_validation/summary.json`; `docs/RGF_B2_VALIDATION.md` |
| B3 | sabit düzensizlik realizasyonları; RGF ile doğrudan seyrek Green referansının karşılaştırması | **FAIL: dx_1_512_vs_1_1024** (10 kriter grubundan 9'u geçti) | S uyumu ≤ 4.9×10⁻¹²; dx = 1/512 → 1/1024 iyileştirmesinde maks \|ΔT₁\| = 5.97×10⁻⁴ > 2×10⁻⁴ (W = 200) | `results/rgf_b3_validation/summary.json` (`overall_B3_status`); `docs/RGF_B3_VALIDATION.md` |
| B3X | 6 konfigürasyon (W ∈ {100, 200} × ζ ∈ {0.02, 0.05, 0.2}), bağımsız dalga fonksiyonu referansı | **FAIL: C11_dx_refinement** (14 kriterden 13'ü geçti) | S uyumu ≤ 5.3×10⁻¹²; C11: 9.17×10⁻⁴ > 2×10⁻⁴ | `results/rgf_b3x_validation/summary.json` (`overall_B3X_status`); `docs/RGF_B3X_VALIDATION.md` |
| B3Y | başarısız B3X konfigürasyonlarında ince ızgara; 4 yöntem, biri genişletilmiş hassasiyette | **PASS** (`overall_status`; `stopped_at = null`) | dört yöntem S'de ≤ 9.1×10⁻¹¹ uyumlu; B3X'teki dx değişimleri dört yöntemde aynıdır; yani bunlar tek bir çözücünün değil, **sonlu-fark kafes probleminin** özelliğidir | `results/rgf_b3y_validation/summary.json`; `docs/RGF_B3Y_VALIDATION.md` Part II |
| B4 | sabit realizasyonda durağan çözüm ile zamana bağlı çözümün karşılaştırması | **YAPILMADI** | — | `docs/PHASE2_METHOD_SELECTION.md` §3 (B4 tanımı); README §12 ("B4 open") |

: Tablo A.1. RGF doğrulama hattının ayrıntılı durumu.

**Testlerin kapsamı:**

- **B1** bir iç tutarlılık testidir; bağımsız fiziksel doğrulama değildir.
- **B2** kanal bağlaşımını sınamaz.
- **B3Y**, kafes probleminin sürekli Schrödinger denklemine yakınsamasını göstermez. Tek bir realizasyonla sınırlıdır.

**Sonuç: Bağımsız durağan doğrulama tamamlanmamıştır.** Bölüm 7'deki taşınım sonuçları **bağımsız bir sabit-enerji
yöntemiyle henüz doğrulanmamıştır.** Bunun için B4 gereklidir.

### A.2 W/Z teşhis çalışmaları

B1-F'nin faz kriterini geçememesi üzerine, üretim RGF çözücüsünün uzun cihazlardaki yuvarlama davranışını anlamak için
ayrı bir teşhis çalışması başlatılmıştır (`experiments/wz_rgf_diagnostic/`).

**W/Z formülasyonu.** Aynı ayrık problemin, küçük büyüklüklü değişkenlerle yazılmış cebirsel bir eşdeğeridir. Δ_j ve W_j
değişkenleri şöyle tanımlanır:

$$ \Delta_j = (E\,I - \mathrm{diag}(E_\perp) - V_j)\,dx^2, \qquad W_j = D_j/t + I. $$

Bu çalışma **fiziksel bir sonuç üretmez.** Yalnızca bir sayısal yöntemin yuvarlama özelliklerini inceler.

| Çalışma | Durum (aynen) | Başlıca ölçüm | Sınırlama | Kaynak |
|---|---|---|---|---|
| B1-F en kötü nokta ayrıştırması | — | 110 basamaklı tam teleskopik ayrıştırma: hatanın +1.6108×10⁻¹⁰'u (≈ %104) saklanan köşegenin binary64 temsilinden gelir | **tek nokta**; tüm satırlar için tek nedene bağlanmamıştır | `experiments/wz_rgf_diagnostic/REPORT.md` §4 |
| Scalar W/Z (tek nokta, 1×1 kanal) | **INCONCLUSIVE** | faz −1.19×10⁻¹⁵ (orijinal +1.5425×10⁻¹⁰); genlik +4.82×10⁻¹⁴ (orijinal −1.08×10⁻¹⁵, 44 kat kötü) | yalnızca kanal 5; kapalı kanallar çalışmadı | `REPORT.md`; düzeltmeler `ERRATUM_worst_point.md` (E1–E3) |
| 63 kanallı genlik yuvarlaması | önceden kayıtlı H_A: SUPPORTED (kanal 5'te bağımsız olarak ayırt edilemedi) | W/Z fazı ≤ 3.24×10⁻¹⁵ (orijinal 1.0×10⁻¹¹–1.54×10⁻¹⁰); W/Z genliği maks 1.93×10⁻¹³ (üretimin ≈ 3.4 katı); \|P\| − 1 = 1.37×10⁻¹⁷, 3482 adımda 4.757×10⁻¹⁴ | temiz cihaz (köşegen bloklar); bağlı bloklar, düzensizlik ve N-ölçeklenmesi test edilmedi | `amplitude_roundoff_report.md`; `results_amplitude_roundoff/amplitude_roundoff_results.json` |
| Coupled-block (birbiriyle değişmeyen bloklar) | **DESIGN INSUFFICIENT** | 76 execution / 72 konfigürasyon; referans kontrolleri ≤ 2.0×10⁻⁶¹ / 3.2×10⁻⁶¹; W/Z uç blokları ≤ 2.5×10⁻¹⁶; 480/480 algılanabilir mutant algılandı; ölçekleme eğimi 0.78–0.79 < 0.8 (önceden kayıtlı alt sınır) | N ≤ 4, N_y ≤ 3; **W/Z FAIL anlamına gelmez**: tasarım karar üretememiştir | `results_coupled_blocks/coupled_blocks_report.md` |
| Coupling-scaling v2 | manifestler **GENERATED AND FROZEN** (`fac5f5e`); ana deney **NOT RUN** | XP 60, HO 120 konfigürasyon; çözücü veya referans çalıştırılmadı | driver kapılıdır: §4.7 4(b) "duplicate set" tanımsız (G-1), R0 imzası yok | `PREDICTIONS_coupling_scaling_v2.md` (Rev 3.7); `results_coupling_scaling_v2/`; `R0_REVIEW_PACKET.md` |
| Uygulama düzeyi R0 | **PENDING** | — | `REVIEW_coupling_scaling_v2.md` yok; kontrol listesi boş | `R0_REVIEW_PACKET.md` |
| T1–T6 yazılım testleri | PASS (yalnızca yazılım ve protokol) | 54/507/28/468/3528/36 kontrol, 0 başarısız; mutasyon testleri 36/36; pytest 691 passed + 2 xfailed | fiziksel doğrulama değildir; test edilen kodla aynı yapay zekâ ekosistemi tarafından yazılmıştır | `R0_REVIEW_PACKET.md` I.7 (resmi koşu commit `1668540`) |

: Tablo A.2. W/Z teşhis çalışmalarının durumu.

**Özet durum:** W/Z **NOT VALIDATED**; coupled-block **DESIGN INSUFFICIENT**; ana coupling-scaling deneyi **NOT RUN**;
R0 **PENDING**.

**Kaldırılan iddia.** Önceki bir belge talebinde geçen "B1 için korelasyon 0.91–0.93, eğim 0.99–1.05" değerleri commit
edilmiş hiçbir kaynakta bulunamamıştır. Bu taslakta **kullanılmamıştır** (`REPORT_READINESS_AUDIT.md` §5, K-5; Ek B).

**Önerilen şekiller (ek için):**

| Şekil | Dosya yolu | Veri kaynağı | Önerilen başlık |
|---|---|---|---|
| ŞEk.1 | `docs/turkce_ozet_figurler/sekil1_b1_faz.png` | `results/rgf_b1_validation/phase_audit_hp_summary.json` | "B1 konfigürasyonlarında maksimum faz hatası; F, önceden belirlenmiş 10⁻¹⁰ kriterini aşar." |
| ŞEk.2 | `docs/turkce_ozet_figurler/sekil2_63kanal.png` | `amplitude_roundoff_results.json` | "63 kanallı tek noktada açık kanal başına faz ve genlik hataları (W/Z ile orijinal RGF)." |
| ŞEk.3 | `docs/turkce_ozet_figurler/sekil3_coupled.png` | `coupled_blocks_comparisons.csv` | "Coupled-block deneyinde W/Z uç blok hatası, sezgisel sınır τ_end ve 10⁻¹⁰ kriteri." |
| ŞEk.4 | `results/rgf_b2_validation/{error_vs_dx,error_vs_energy,order_vs_energy}.png` | `results/rgf_b2_validation/` | "B2: RGF ile 1B ODE arasındaki iletim farkı ve yakınsama derecesi." |

> **Taslak notu (Ek A).**
> - **Planlama:** `docs/PHASE2_METHOD_SELECTION.md` §2–§3, §8.
> - **Doğrulama raporları:** `docs/RGF_B1_VALIDATION.md`, `RGF_B2_VALIDATION.md`, `RGF_B3_VALIDATION.md`,
>   `RGF_B3X_VALIDATION.md`, `RGF_B3Y_VALIDATION.md`.
> - **Sonuç dosyaları:** `results/rgf_b{1,2,3,3x,3y}_validation/summary.json`.
> - **W/Z teşhisi:** `experiments/wz_rgf_diagnostic/{REPORT.md, ERRATUM_worst_point.md, amplitude_roundoff_report.md,
>   results_coupled_blocks/coupled_blocks_report.md, PREDICTIONS_coupling_scaling_v2.md, R0_REVIEW_PACKET.md}`.
> - **Denetim:** `REPORT_READINESS_AUDIT.md` §4.

---

## Ek B. R0 → R0.1 revizyon kaydı (iç not)

**Değişiklikler:**

1. **Başlık.** TÜBİTAK 2204-A kuralına uygun 11 kelimelik başlık kullanıldı (§1). R0'daki uzun başlık kaldırıldı.
2. **Künye.** Yazar, okul ve danışman için yer tutucular eklendi; hiçbir kişisel bilgi uydurulmadı.
3. **V6.** Kaynağı bulundu: `results/validation/validation_summary.json` (`v6`), `validate_solver.py::v6_conservation`.
   Değerler 1.51×10⁻¹² ve 1.60×10⁻¹³'tür. R0'daki "kaynak yok" işareti kaldırıldı.
4. **PINN V4.** Kaynağı bulundu: `v4_space_convergence.csv`, `pinn_W112_psi`, dx = 1/128 satırında 1.90×10⁻⁶.
5. **PINN "kalıntı 0.22".** Kaynağı bulunamadı ve çıkarıldı. Yerine başka bir değer konmadı.
6. **Özgün motivasyon.** Belgelenmiş git geçmişi ile öğrencinin anlatımı ayrıldı. Projenin evrimi temkinli biçimde
   anlatıldı (§4.5).
7. **CAP tanımı.** Koda göre düzeltildi: s, katman genişliğine bölünmüş boyutsuz derinliktir (§6.5).
8. **Doğrulama türleri.** İç doğrulama, yakınsama ve bağımsız durağan doğrulama ayrıldı (§6.6). B2'deki TDSE ile 1B ODE
   karşılaştırmasının düzensizlik içermediği belirtildi.
9. **RGF/W-Z.** Ana metinde kısaltıldı (§9). Ayrıntılı tablolar ve şekiller Ek A'ya taşındı. Tüm durumlar korundu.
10. **Klasik yasak alan oranları.** "%15 / %25" değerleri Gauss tek-nokta kestirimi olarak açıklandı ve ölçülmüş değer
    olmadığı belirtildi (§7.4, §7.5).
11. **PR'nin mutlak değeri.** Neden 1 olmadığı `observables.py` belgesine dayanarak açıklandı (§7.5).
12. **Kaynaklar.** Adam [37], Holm [38], Welch [39] eklendi (Welch DOI'si doğrulanmadı). Bootstrap p-değeri kaynakları
    [40, 41] depo belgesinden eklendi.
13. **Commit bilgileri.** Denetim, R0 ve R0.1 commit'leri ayrıldı (başlık bloğu, §13).

**Denklem kontrolleri** (koda karşı):

| Denklem | Kontrol edilen kaynak | Sonuç |
|---|---|---|
| (1) TDSE | `reference_solver.py` başlık belgesi; CAP'lı biçim −iW_cap | uyumlu |
| (2) ve §6.3 başlangıç paketleri | `initial_condition.py` (A = 1/√(πσ_xσ_y); (πσ_x²)^{−1/4}√(2/L_y)) | uyumlu |
| (3) düzensizlik | `disorder.py` (√(2/M), q = g/ζ, φ ~ U[0, 2π), tanh penceresi) | uyumlu |
| (5) Born oranları | `disorder.born_rates` (`ive` = e^{−a}I_n(a); ℓ = k/oran) | uyumlu |
| (6) difüzif tahmin | `analyze_transport.py` (ℓ = (π/2)ℓ_tr, T = ℓ/(L+ℓ); enerji ortalaması; fitler 1/SE² ağırlıklı) | uyumlu |
| (7) PINN kalıntıları ve sert kısıt | `pinn_core.py` ((T/A) ölçekleme, normalizasyon kaybı); `model.py` (B = 16x̂(1−x̂)ŷ(1−ŷ)) | uyumlu |
| CAP | `reference_solver.py::cap_profile` | **R0'da yanlış ifade edilmişti; düzeltildi** |
| T, R, PR | `transport.py`, `observables.py` | uyumlu |

**Kanıt yetersizliği nedeniyle çıkarılan iddialar:**

- "20000 taze noktada kalıntı 0.22" (PINN sabit kolokasyon koşusu);
- "B1 için korelasyon 0.91–0.93, eğim 0.99–1.05" (R0'da da kullanılmamıştı; kayıt için).

> **Taslak notu (Ek B).** `git log`; `results/validation/validation_summary.json`; `results/validation/v4_space_convergence.csv`;
> `results/validation_log.txt`; `results/pinn/diag_fixed_points_lbfgs/train_info.json`; kod dosyaları tablo içinde.
