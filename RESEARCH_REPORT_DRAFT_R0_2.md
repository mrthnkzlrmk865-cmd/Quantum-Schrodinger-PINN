# Araştırma Raporu — Taslak R0.2

> **Taslak durumu.** Bu bir **iç revizyon taslağıdır**, yayıma hazır bir metin değildir.
>
> - **Ne yapıldı.** R0.2, R0.1 taslağının (`RESEARCH_REPORT_DRAFT_R0_1.md`) bilimsel anlatı ve yapı revizyonudur. R0 ve
>   R0.1 değiştirilmemiştir.
> - **Kapsam.** Bu taslak yazılırken hiçbir simülasyon, eğitim, benchmark veya doğrulama betiği çalıştırılmamıştır. Tüm
>   sayılar commit edilmiş dosyalardan alınmıştır ve R0.1'deki doğrulanmış değerlerle aynıdır.
> - **İşaretler.**
>   - **[TO BE COMPLETED]:** eksik bilgi;
>   - **[REFERENCE DETAILS TO VERIFY]:** bibliyografik ayrıntısı doğrulanamayan kaynak;
>   - **[REVIEW]:** çözülmemiş bir tutarsızlık.
> - **Taslak notları.** Her ana bölümün sonundaki "Taslak notu" kutuları kaynak dosyaları listeler. Bu kutular revizyon
>   içindir ve nihai metinden çıkarılabilir.
> - **Depo durumu** (R0.2 yazılırken):
>
>   | belge | commit |
>   |---|---|
>   | rapor hazırlık denetimi (`REPORT_READINESS_AUDIT.md`) | `e7aa728` |
>   | R0 taslağı | `11b342c` |
>   | R0.1 taslağı | `509011c` (R0.2 yazılırken HEAD) |
>   | `origin/main` | `7fec08a` (PR #27 birleştirmesi; dosya ağacı `509011c` ile aynı) |
>   | R0.2 taslağı (bu dosya) | bu dosyayı ekleyen commit (teslim notunda verilir) |
>
>   Dal: `claude/quantum-schrodinger-audit-iy6kn4`. R0.1 → R0.2 değişiklikleri Ek F'dedir.

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
| Yarışma | TÜBİTAK 2204-A Lise Öğrencileri Araştırma Projeleri Yarışması, 2026–2027 |
| Alan / kategori | [TO BE COMPLETED — depoda yok] |

> **İç not: TÜBİTAK başlık uyumu.**
>
> - **Kural.** Başlık en fazla 15 kelime olmalıdır.
> - **Sayım.** Kelimeler boşlukla ayrılarak sayılmıştır:
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
> - **Uyum.** Başlık 2B modeli ve ana soruyu (ζ'nin taşınıma etkisi) tanımlar. 3B model, PINN üstünlüğü veya Anderson
>   lokalizasyonunun kanıtı ima edilmez.
> - **Tek resmî başlık.** Bu dosyada 15 kelimeyi aşan başka hiçbir başlık resmî başlık olarak yer almaz.

---

## 2. Özet

**Problem.** Bu çalışma, iki boyutlu (2B) bir dalga kılavuzunun sonlu bir bölgesine yerleştirilmiş rastgele bir
potansiyelin, içinden geçen bir kuantum dalga paketinin **iletimini** nasıl etkilediğini sayısal olarak inceler. İletim,
paketin düzensiz bölgeyi geçme olasılığıdır.

- **Potansiyel.** Gauss-korelasyonlu, statik bir rastgele alandır. Ana değişken, potansiyelin **korelasyon uzunluğu ζ**'dir;
  ζ kabaca rastgele tepe ve çukurların tipik genişliğidir.
- **Destekleyici taramalar.** Düzensizlik gücü W, gelen dalga sayısı k ve düzensiz bölge uzunluğu L_d da taranmıştır.

**Yöntem.**

- Zamana bağlı Schrödinger denklemi, ayrık sinüs dönüşümüne (DST) dayalı bir Strang adım-bölme (*split-step*) çözücüyle
  çözülmüştür.
- Çözücü üç tür sınamadan geçmiştir:
  - **İç doğrulama:** analitik çözümler ve korunum yasaları;
  - **yakınsama testleri:** zaman adımı, ızgara aralığı ve düzensizlik mod sayısı;
  - **bağımsız durağan karşılaştırma:** yalnızca düzensizlik içermeyen bir bariyer için yapılmıştır.
- **Düzensiz taşınım için bağımsız durağan çapraz doğrulama (B4) yapılmamıştır.**

**Ana sonuç.** Sabit düzensizlik gücünde (W = 100), 48 rastgele realizasyonun ortalaması olan iletim ⟨T⟩, ζ'nin monoton
olmayan bir fonksiyonudur:

| ζ | kζ | ⟨T⟩ |
|---|---|---|
| 0.02 | 0.40 | 0.670 |
| 0.05 | 1.01 | 0.494 (en küçük değer) |
| 0.20 | 4.05 | 0.734 |

**Ana sonucun yorumu ve sınırları.**

- **Kuramla ilişkisi.** Bu davranış, Born saçılma kuramının taşınım ortalama serbest yolu için öngördüğü minimumla
  niteliksel olarak uyumludur.
- **Sınırları.**
  - Yalnızca yedi ζ noktası örneklenmiştir.
  - Eşleştirilmiş istatistiksel karşılaştırmalar sonuçlar görüldükten sonra seçilmiştir (*post-hoc*).
  - Önceden kayıtlı doğrulayıcı kampanya çalıştırılmamıştır.
  - Bu nedenle minimumun varlığı **keşifsel** düzeyde desteklenmiş, konumu ise yalnızca **kaba** biçimde (kζ ≈ 0.6–1.4)
    belirlenmiştir.

**Destekleyici sonuçlar.**

- **W:** ⟨T⟩, W ile monoton azalır.
- **Born karşılaştırması:** basit bir Born + difüzyon tahmini, zayıf saçılmada ölçümü yaklaşık %12 içinde yakalar. Güçlü
  saçılmada nedeni belirlenmemiş sistematik bir sapma vardır.
- **k:** daha hızlı paketler daha az saçılır.
- **L_d:** güçlü düzensizlikte iletimin bölge uzunluğuyla azalması üstel biçime daha iyi uyar. Bu, lokalizasyonla
  **tutarlıdır ancak kanıt değildir.**
- **Kapalı kutu:** 2B kapalı kutularda bu boyutlarda lokalizasyon saptanmamıştır.

**İkincil yöntemsel sonuç.** Fizik-bilgili sinir ağı (PINN), mevcut hesaplama bütçesinde klasik çözücüye yakın bir doğruluğa
ulaşamamıştır. En kolay ayarda bile hata yaklaşık %49'dur.

---

## 3. Anahtar kelimeler

Kuantum taşınımı; korelasyonlu düzensizlik; korelasyon uzunluğu; zamana bağlı Schrödinger denklemi; dalga kılavuzu;
split-step spektral yöntem; Born yaklaşımı; fizik-bilgili sinir ağları

---

## 4. Giriş

### 4.1 Arka plan ve motivasyon

Düzensizlik doğada ve teknolojide yaygındır: iletkenlerdeki safsızlıklar, dalga kılavuzlarındaki pürüzlülük, soğuk atom
deneylerinde kullanılan rastgele ışık desenleri (*speckle*) [9, 13].

**Girişim ve taşınım rejimleri.** Bir dalga rastgele bir ortamda defalarca saçıldığında saçılan dalgalar birbirleriyle
girişim yapar. Düzensizliğin gücüne, boyutuna ve uzaysal yapısına bağlı olarak taşınım üç biçimde olabilir [1–5, 8]:

| Rejim | Dalga ne yapar? |
|---|---|
| balistik | neredeyse saçılmadan geçer |
| difüzif | çok kez saçılarak, rastgele yürüyüşe benzer biçimde ilerler |
| Anderson lokalizasyonu | girişim nedeniyle bir bölgeye hapsolur [1] |

Bu olgular elektronlarda, ışıkta ve madde dalgalarında gözlenmiştir [8, 10–12].

**Korelasyon uzunluğu.** Korelasyon uzunluğu ζ, rastgele potansiyelin tepe ve çukurlarının tipik genişliğidir. Önemli olan
ζ'nin dalga boyuna göre büyüklüğüdür; bu, kζ boyutsuz sayısıyla ölçülür [9]:

- kζ ≪ 1: potansiyel dalga boyuna göre çok ince taneli bir "gürültü" gibidir;
- kζ ≫ 1: potansiyel dalga boyuna göre düzgün, yavaş değişen bir alandır.

Bu iki uç durumda dalga farklı biçimde saçılır.

**Bu çalışmanın yaklaşımı.**

- Zamana bağlı dinamiğin doğrudan simülasyonu, bir dalga paketinin bu ortamdan geçişini izlemeye ve iletim istatistiğini
  çok sayıda rastgele potansiyel örneği üzerinden ölçmeye olanak verir.
- Karşılaştırma için kullanılan kuram yerleşiktir; bu çalışma yeni bir kuram önermez.
- Sayısal deney, bu kuramın öngörülerinin belirli ve sonlu bir geometride nasıl göründüğünü inceler.

### 4.2 Literatürden alınan bağlam

Aşağıdakiler yerleşik bilgilerdir; bu projenin bulguları değildir.

- **Ölçekleme kuramı.** Zamanı tersine çevirmeye göre simetrik sistemlerde, 1B ve 2B'de keyfi zayıf korelasyonsuz
  düzensizlik için bütün durumlar lokalizedir [2].
- **2B'de lokalizasyon uzunluğu.** Zayıf düzensizlikte, dalganın hapsolduğu bölgenin tipik boyu ξ ≈ ℓ exp(πkℓ/2) kadar
  büyüktür; burada ℓ ortalama serbest yoldur [3]. Bu nedenle kℓ ~ 1 olmadıkça sonlu örneklerde gözlenmez.
- **Az kanallı dalga kılavuzları** (yarı-1B teller). N açık kanallı bir kılavuzda lokalizasyon uzunluğu Nℓ mertebesindedir
  [6]. Taşınım şöyle ayrılır:

  | Rejim | Koşul | Ortalama iletim |
  |---|---|---|
  | balistik | L ≪ ℓ | — |
  | difüzif | ℓ ≪ L ≪ Nℓ | ⟨T⟩ ∝ ℓ/L |
  | lokalize | L ≫ Nℓ | ln T ∝ −L/ξ; T geniş, log-normal dağılır [6] |

### 4.3 Araştırma sorusu

**Ana soru.**

> Gauss-korelasyonlu statik düzensizliğin korelasyon uzunluğu ζ, sonlu iki boyutlu sert duvarlı bir dalga kılavuzundaki
> düzensiz bölgeden geçen kuantum dalga paketinin iletimini nasıl etkiler?

**Destekleyici sorular.**

1. Gözlenen ζ bağımlılığı, Born saçılma kuramının öngörüsüyle niteliksel olarak uyumlu mu?
2. Düzensizlik gücü W, gelen dalga sayısı k ve bölge uzunluğu L_d iletimi nasıl etkiler?
3. Erişilebilir sistem boyutlarında lokalizasyon belirtileri sonlu-boyut etkilerinden ayırt edilebilir mi?

**İkincil yöntemsel soru.**

> Aynı 2B problemde fizik-bilgili sinir ağı (PINN) çözümü, klasik çözücüye göre ne kadar doğrudur ve bu doğruluğu ne
> sınırlar?

`README.md` §1'deki soru W ve ζ'yi birlikte sorar. Bu rapor, başlıkla tutarlı olarak ζ'yi ana soru alır ve diğerlerini
destekleyici sonuç olarak sunar. Bu bir kapsam seçimidir; README ile çelişmez.

### 4.4 Amaçlar ve kapsam

**Amaçlar:**

1. ⟨T⟩'nin ζ'ye bağımlılığını güven aralıklarıyla ölçmek (ana amaç).
2. Bu bağımlılığı Born + difüzyon tahminiyle karşılaştırmak.
3. W, k ve L_d bağımlılıklarını destekleyici bağlam olarak ölçmek.
4. Erişilebilir boyutlarda lokalizasyon belirtilerini sınamak.
5. PINN'in aynı tür problemdeki doğruluğunu ve hata kaynaklarını nicel olarak belirlemek.

**Kapsam dışı:**

- parçacıklar arası etkileşimler;
- zamanla değişen düzensizlik;
- faz bozulması (*dephasing*) veya çevreyle etkileşim;
- üç boyut;
- deneysel doğrulama (README §10).

### 4.5 Projenin gelişimi

Bu bölüm belgelenmiş depo geçmişini, öğrencinin anlatımını ve güncel çalışmayı ayrı tutar.

**(a) Belgelenmiş depo geçmişi** (git geçmişinden):

| tarih | commit | belgelenen içerik |
|---|---|---|
| 2026-07-18 | `00a52e9` | "2B zamana bağlı Schrödinger denkleminin PINN ile çözümü, etkileşimli 3B Plotly görselleştirmesiyle" |
| 2026-07-18 | `7b9f41a` | \|Ψ\|²'nin 3B **yüzey görselleştirmesi**; çözülen denklem 2B'dir |
| 2026-07-18 | `5f56539` | motivasyon: "makine öğrenmesinin kuantum kimyası ve moleküler biyoloji simülasyonlarını nasıl hızlandırabileceğini keşfetmek"; problem: 2B bariyerden tünelleme |
| 2026-08-08 | `4fe3eb2` | DST tabanlı klasik referans çözücü ve "ENAQT-esinli" düzensizlik deneyi |
| 2026-09-01 | `42cea06` | denetlenen sürüm: metin 1B ve periyodik sınırlar anlatır, kod ise 2B ve sert duvarlıdır (`docs/AUDIT.md` 1.1) |

**Denetimin bulguları** (`docs/AUDIT.md`):

- ENAQT (çevre-destekli kuantum taşınımı) motivasyonu modelle uyuşmuyordu; modelde çevre veya faz bozulması yoktur.
- "Düzensizlik iletimi artırır" sonucu, kapalı kutuda alınmış tek bir anlık görüntünün yapaylığıydı.
- Özgün referans çözüm yakınsamamıştı (README §7, §9.9).

**(b) Öğrencinin anlatımı (belgelenmemiş).** Öğrenciye göre çıkış noktası, Schrödinger denkleminin daha yüksek boyutlu
uzaylarda makine öğrenmesiyle modellenmesine duyulan ilgidir.

- Depoda bu ifade açıkça yer almaz.
- "3B" sözcüğü ilk commit'lerde yalnızca görselleştirme için geçer.
- [TO BE COMPLETED: Tarihli bir kaynakla desteklenmedikçe nihai metinde "kişisel motivasyon" olarak sunulmalıdır.]

**(c) Güncel çalışma.**

- Statik düzensizlikte üniter (olasılığı koruyan) dinamik üzerinedir.
- **İki boyutludur**; üç boyutlu Schrödinger denklemi çözülmemiştir.
- Ana sonuçlar klasik çözücüyle elde edilmiştir. PINN ikincil bir yöntemsel değerlendirmedir.
- Yüksek boyutlara genelleme, bu çalışmanın verisiyle desteklenmeyen bir gelecek çalışma konusudur.

> **Taslak notu (Bölüm 1–4).** `README.md` §1–§2, §7, §9.9, §10, "Project history"; `docs/PHASE2_METHOD_SELECTION.md`
> (giriş); `docs/AUDIT.md` 1.1, 2.2–2.4; `docs/LITERATURE.md`; `git log --reverse`; `git show 00a52e9:README.md`,
> `5f56539:README.md`, `4fe3eb2:README.md`, `42cea06:README.md`, `7b9f41a:wave_function.py`.

---

## 5. Kuramsal arka plan

### 5.1 Birimler ve semboller

Boyutsuz birimler kullanılmıştır (ħ = m = 1). Uzunluklar dalga kılavuzu genişliği L_y = 1 birimindedir.

| Sembol | Anlamı |
|---|---|
| ψ(x, y, t) | dalga fonksiyonu; \|ψ\|² olasılık yoğunluğudur |
| V(x, y) | potansiyel |
| W | düzensizliğin RMS (karekök ortalama kare) gücü: potansiyelin tipik genliği |
| ζ | korelasyon uzunluğu: rastgele tepe ve çukurların tipik genişliği |
| k_x0 | paketin x yönündeki ortalama dalga sayısı |
| k | toplam dalga sayısı, k = (k_x0² + π²)^{1/2}; π terimi 1. enine moddan gelir |
| E₀ = k²/2 | paketin ortalama enerjisi |
| L_d | düzensiz bölgenin uzunluğu |
| T, R | iletim ve yansıma olasılıkları |
| ⟨·⟩ | rastgele realizasyonlar üzerinden ortalama |
| ℓ_s, ℓ_tr | saçılma ve taşınım ortalama serbest yolları |

### 5.2 Zamana bağlı Schrödinger denklemi

Dalga fonksiyonu zamana bağlı Schrödinger denklemini (TDSE) sağlar:

$$ i\,\partial_t\psi = -\tfrac{1}{2}\left(\partial_x^2+\partial_y^2\right)\psi + V(x,y)\,\psi \qquad (1)$$

- Sağ taraftaki ilk terim kinetik enerji, ikinci terim potansiyel enerjidir.
- Soğurucu sınırlar yoksa ∫|ψ|² dx dy = 1 korunur: toplam olasılık sabittir.

### 5.3 Dalga paketi ve sert duvarlı dalga kılavuzu

**Kılavuz.** y ∈ [0, 1] aralığında ψ = 0 koşulu uygulanır; bu bir "sert duvar"dır, dalga duvardan geçemez. Bu koşul
y yönünde ayrık **enine modlar** oluşturur:

- modlar: χ_n(y) = √2 sin(nπy);
- enine enerjiler: E_⊥,n = (nπ)²/2.

Paketin enerjisi E_⊥,n'den büyük olan modlara **açık kanal** denir; dalga yalnızca bu kanallarda ilerleyebilir. Üretim
ayarında k_x0 = 20 için k ≈ 20.25 ve E₀ ≈ 205'tir; bu **6 açık kanal** verir (README §3).

**Gelen paket.** Enine 1. modda bir Gauss paketidir:

$$ \psi_0 = (\pi\sigma_x^2)^{-1/4} e^{-(x-x_0)^2/2\sigma_x^2} e^{i k_{x0} x}\,\sqrt{2}\,\sin(\pi y) \qquad (2)$$

- σ_x paketin x yönündeki genişliğidir; momentum yayılımı σ_k = 1/(√2 σ_x)'tir.
- Paket tek bir enerjiyi değil bir **enerji bandını** taşır. Bu yüzden ölçülen T ve R enerji üzerinden ortalanmış
  niceliklerdir (README §10).

### 5.4 Gauss-korelasyonlu düzensizlik

Düzensizlik, bir tohum (*seed*: rastgele sayı üretecinin başlangıç değeri), W ve ζ ile tanımlanan **sürekli** bir rastgele
alandır. Shinozuka–Deodatis yöntemiyle, rastgele yönlü M dalganın toplamı olarak üretilir [14]:

$$ V_{\mathrm{dis}}(\mathbf r) = W\sqrt{\frac{2}{M}}\sum_{j=1}^{M}\cos(\mathbf q_j\cdot\mathbf r+\varphi_j)\,S(x),\qquad \mathbf q_j\sim\mathcal N(0,\zeta^{-2}\mathbb{1}_2),\ \varphi_j\sim U[0,2\pi) \qquad (3)$$

**Bileşenler:**

- q_j rastgele dalga vektörleridir. Ortalaması 0 ve kovaryansı ζ⁻² olan Gauss dağılımından çekilir; küçük ζ daha büyük
  dalga vektörleri, yani daha ince tanecikler verir.
- φ_j rastgele fazlardır.
- M = 512 kullanılmıştır.
- S(x) = ½[tanh((x−a)/δ) − tanh((x−b)/δ)] (δ = 0.02), düzensizliği a < x < b bölgesiyle sınırlayan yumuşak bir penceredir.

**İstatistiksel özellikler.** Çok sayıda realizasyonun ortalamasında:

$$ \langle V_{\mathrm{dis}}\rangle = 0,\qquad \langle V_{\mathrm{dis}}(\mathbf r)V_{\mathrm{dis}}(\mathbf r')\rangle = W^2 e^{-|\mathbf r-\mathbf r'|^2/2\zeta^2}. \qquad (4)$$

Yani iki noktadaki potansiyel değerleri, aralarındaki uzaklık ζ'den küçükse güçlü biçimde ilişkilidir; çok daha büyükse
bağımsızdır.

- Aynı tohum her ızgara çözünürlüğünde aynı fiziksel alanı verir.
- Klasik çözücü ve PINN aynı alanı kullanır; iki hesap 10⁻⁹ bağıl düzeyde uyuşur (README §4; `tests/test_potential.py`).

### 5.5 İletim ve yansıma

Düzensiz bölge a < x < b aralığındadır. Ölçüm, paket bölgeyi geçtikten sonraki bir t_end anında yapılır (`transport.py`).

| Nicelik | Tanım | Anlamı |
|---|---|---|
| T | sağ uçta soğurulan olasılık + P(x > b) | paketin bölgeyi geçen kesri |
| R | sol uçta soğurulan olasılık + P(x < a) | geri dönen kesir |
| P_in | 1 − T − R | t_end'de hâlâ bölgede kalan kesir |

P_in küçükse T ve R, paketin enerji bandı üzerinden ortalanmış saçılma olasılıklarıdır.

### 5.6 Born saçılma kuramı ve boyutsuz parametreler

**Born yaklaşımı.** Zayıf düzensizlikte bir dalganın saçılma hızını veren birinci mertebe pertürbasyon kuramıdır. Fermi'nin
altın kuralına göre, sınırsız 2B Gauss-korelasyonlu düzensizlikte k dalga sayılı bir parçacık için iki oran tanımlanır
(`disorder.born_rates`; README §4):

- 1/τ_s, herhangi bir yöne saçılma oranıdır;
- 1/τ_tr, yönü değiştiren saçılmanın etkin oranıdır (taşınım oranı).

$$ \frac{1}{\tau_s}=2\pi W^2\zeta^2 e^{-a}I_0(a),\qquad \frac{1}{\tau_{tr}}=2\pi W^2\zeta^2e^{-a}\left[I_0(a)-I_1(a)\right],\qquad a = k^2\zeta^2 \qquad (5)$$

Burada I₀ ve I₁ değiştirilmiş Bessel fonksiyonlarıdır. Ortalama serbest yollar ℓ_s = kτ_s ve ℓ_tr = kτ_tr'dir. ℓ_tr,
dalganın yönünü "unutmadan" önce aldığı tipik yoldur.

**Kuramın ζ öngörüsü** (kuram, bu çalışmanın bulgusu değil). Sabit W'da ℓ_tr, ζ'nin monoton olmayan bir fonksiyonudur
(README §9.3; `disorder.born_rates` belgesi):

- kζ ≪ 1'de düzensizliğin toplam saçılma gücü ζ² ile küçülür ("beyaz gürültü" sınırı); saçılma zayıftır ve ℓ_tr büyüktür.
- kζ ≫ 1'de düzgün düzensizlik dalgayı ağırlıklı olarak ileri yönde saçar. Yön az değişir, ℓ_tr yine büyüktür.
- Arada, kζ ≈ 0.8 civarında ℓ_tr en küçüktür; saçılma en etkilidir.

**Basit difüzyon tahmini.** Uzunluğu L_d olan bölgenin iletimi için:

$$ T_{\mathrm{est}} \approx \frac{\ell}{L_d+\ell}, \qquad \ell = \frac{\pi}{2}\,\ell_{tr}. \qquad (6)$$

Bu tahminin sınırları:

- Yalnızca zayıf saçılmada (kℓ ≫ 1) geçerlidir.
- Dalga kılavuzunun ayrık kanal yapısını içermez.
- π/2 öneki konvansiyona bağlıdır.

Bu nedenle **mertebe kıyası** olarak kullanılır ve yalnızca eğilimleri yorumlanır.

**Boyutsuz parametreler:**

| Parametre | Fiziksel anlamı |
|---|---|
| kζ | tanecik boyu / dalga boyu |
| kℓ_tr | saçılma gücü; kℓ ≫ 1 zayıf, kℓ ~ 1 güçlü saçılma |
| W/E₀ | potansiyel genliğinin paket enerjisine oranı |
| L_d/ℓ, L_d/ξ | bölge uzunluğunun taşınım ölçeklerine oranı |

> **Taslak notu (Bölüm 5).** `README.md` §3, §4, §8, §9.3, §10; `config.py`; `disorder.py` (`born_rates`,
> `make_disorder`); `initial_condition.py`; `transport.py`.

---

## 6. Fiziksel model ve sayısal yöntem

### 6.1 Geometriler

| Geometri | Kullanım | Alan | Başlangıç durumu | Sınırlar |
|---|---|---|---|---|
| Açık dalga kılavuzu | ana taşınım deneyleri | L_y = 1, L_x = 2.9 + L_d | (2), x₀ = 1.2, σ_x = 0.2, k_x0 = 20 | y'de sert duvar; x uçlarında soğurucu katman |
| Kapalı kutu | sonlu-boyut testi | L × L, L ∈ {1, 2} | 2B Gauss, merkezde, σ = 0.08, k_x0 = 20 | dört duvar sert; düzensizlik her yerde |
| PINN problemi | PINN ile klasik çözücünün karşılaştırması | 1 × 1, t ∈ [0, 0.02] | 2B Gauss, (0.35, 0.5), σ = 0.07, k_x0 = 15 | dört duvar sert |

**Dalga kılavuzu yerleşimi** (x ekseni boyunca):

1. sol soğurucu katman, [0, 0.6];
2. paket merkezi, x₀ = 1.2;
3. düzensiz bölge, [1.9, 1.9 + L_d];
4. 0.4'lük boşluk;
5. sağ soğurucu katman, genişlik 0.6.

**Sayısal ayarlar:** ızgara aralığı dx = dy = 1/64, zaman adımı Δt = 2×10⁻⁴ (`config.py::waveguide_config`).

### 6.2 Soğurucu sınırlar

Açık bir kılavuzu sonlu bir alanda modellemek için uçlara karmaşık soğurucu potansiyel (CAP) konur. CAP, gelen dalgayı geri
yansıtmadan yutması amaçlanan sanal bir sürtünme terimidir:

- Denklem (1)'e −iW_CAP(x)ψ terimi eklenir.
- W_CAP = W₀ s(x)², W₀ = 300'dür. s(x) katman içindeki derinliğin katman genişliğine oranıdır, yani 0 ≤ s ≤ 1
  (`reference_solver.py::cap_profile`) [19].
- Her uçta soğurulan olasılık ayrı ayrı kaydedilir. Böylece T + R + P_in = 1 hesap hatası düzeyinde sağlanır.

### 6.3 DST tabanlı Strang split-step çözücü

**Temel fikir.** Kinetik ve potansiyel adımlar ayrı ayrı uygulanır; bu adım bölme yöntemidir [15–17].

- **Kinetik adım.** Sert duvarlı bir kutuda sinüs fonksiyonları kinetik enerji operatörünün kesin özfonksiyonlarıdır. Bu
  nedenle kinetik adım ayrık sinüs dönüşümüyle (DST-I) bu temsil için tam olarak uygulanır ve duvarda ψ = 0 kendiliğinden
  sağlanır.
- **Potansiyel adım.** Potansiyel adımı uzayda noktasaldır.
- **Strang sırası.** Adımlar yarım-tam-yarım sırasıyla uygulanır. Bu sıra Δt'de ikinci dereceden doğruluk verir [17, 18].

Birden fazla realizasyon tek seferde (yığın hâlinde) hesaplanır (`reference_solver.py`; README §5).

### 6.4 Doğrulama: üç ayrı tür

Bu raporda üç doğrulama türü ayrı tutulur. Sayısal ayrıntılar Ek A'dadır.

| Tür | Ne sınar? | Testler | Sonuç (özet) |
|---|---|---|---|
| **A. İç çözücü doğrulaması** | çözücü, kesin bilinen çözümleri ve korunum yasasını doğru üretiyor mu? | V1 özmod, V2 serbest Gauss paketi, V6 olasılık korunumu | hatalar 4.6×10⁻¹⁴, 9.4×10⁻¹⁰ ve ≤ 1.5×10⁻¹² |
| **B. Yakınsama ve ayrıklaştırma kontrolleri** | Δt, dx ve M inceltildiğinde sonuç ne kadar değişiyor? CAP ne kadar sahte yansıma üretiyor? | V3 zaman adımı, V4 ızgara, V5 CAP geri dönüşü, V7 düzensizlik mod sayısı | üretim ayarlarında iletimdeki ayrıklaştırma hatası ~10⁻⁵; CAP geri dönüşü k_x0 ≥ 20 için ≤ 5.2×10⁻⁴, k_x0 = 12'de 1.3×10⁻² |
| **C. Bağımsız durağan çapraz doğrulama** | kod paylaşmayan, sabit enerjide çalışan bir yöntem aynı iletimi veriyor mu? | B2 (düzensizlik **yok**): TDSE ile bağımsız 1B denklem çözümü; B4 (düzensizlik **var**) | B2 geçti; **B4 yapılmamıştır** |

**Bu doğrulamaların kapsamı.**

- A ve B, çözücünün **kendi** doğruluğunu ve ayrıklaştırma hatasını ölçer. Sonuç: üretim ayarlarındaki ayrıklaştırma
  hatası (~10⁻⁵), istatistiksel belirsizliğin (~10⁻²) yaklaşık üç mertebe altındadır.
- C türünde yalnızca düzensizlik içermeyen bir bariyer için karşılaştırma yapılmıştır.
- **Düzensizlik içeren taşınım sonuçları, durağan bir yöntemle henüz bağımsız olarak çapraz doğrulanmamıştır (B4 yapılmadı).**
- Raporun ana sonuçları bu nedenle A ve B türü doğrulamaya dayanır.
- Bilinen tek ihmal edilemez sistematik hata, düşük momentumda CAP'ın sahte geri dönüşüdür. k_x0 = 12'de bu geri dönüş R'yi
  yaklaşık 0.01 yukarı kaydırır.

### 6.5 İstatistiksel yöntem

**Realizasyonlar ve tohumlar.** Her parametre noktası n bağımsız düzensizlik realizasyonu (rastgele potansiyel örneği)
üzerinden ortalanır. Realizasyon i her noktada tohum i'yi kullanır.

- Bu **eşleştirilmiş tasarımda** farklı noktalar aynı rastgele çekilişleri paylaşır; yalnızca W ve ζ ile ölçeklenirler.
- Bu, komşu noktalar arasındaki farkların gürültüsünü azaltır. Bedeli şudur: bir tarama içindeki noktalar istatistiksel
  olarak bağımsız değildir.

**Raporlanan büyüklükler** (`stats_utils.py`):

- ortalama;
- örneklem standart sapması (sd);
- standart hata (SE = sd/√n);
- %95 güven aralığı (CI): yüzdelik bootstrap aralığı [31, 35].

Güven aralıkları her nokta için ayrı ayrı geçerlidir.

**Yeterli örnek sayısı.** W = 100, ζ = 0.05, L_d = 1'de 192 realizasyonla:

- ⟨T⟩ = 0.509 [0.500, 0.519];
- SE, 1/√n ile azalır: n = 48'de 0.010, n = 192'de 0.005;
- 48 realizasyonluk ortalama (0.494) tam ortalamayla 1.5 SE içinde uyumludur.

24–48 realizasyonlu taramalar ⟨T⟩'de ≳ 0.03'lük farkları ayırt edebilir (`results/transport/convergence/summary.csv`;
`figures/realization_convergence.png`; README §9.8).

**Analizin niteliği.** Bölüm 7'deki bütün analizler **keşifseldir** (*exploratory*).

- Önceden kaydedilmiş doğrulayıcı (*confirmatory*) bir kampanya tasarlanmış, ancak dondurulmamış ve çalıştırılmamıştır
  (`docs/PHYSICAL_CAMPAIGN_PLAN.md`, revizyon 3).
- Gereken kök entropi ve C0 tohumu kaydedilmemiştir (`docs/CAMPAIGN_STATUS_CONTRACT.md`).

### 6.6 PINN yöntemi (özet)

PINN, çözümü bir sinir ağının çıktısı olarak temsil eder. Ağ, denklemin kalıntısını (*residual*: denklemin ne kadar
sağlanmadığının ölçüsü) seçilmiş noktalarda küçültecek biçimde eğitilir [20, 21].

**Kalıntılar.** ψ = u + iv yazıldığında denklem (1) iki gerçel denkleme ayrılır:

$$ f_u \equiv \partial_t u + \tfrac{1}{2}\nabla^2 v - V v = 0, \qquad f_v \equiv \partial_t v - \tfrac{1}{2}\nabla^2 u + V u = 0. \qquad (7)$$

**Ağ ve kısıtlar:**

- Mimari: 4 gizli katman × 64 nöron.
- **Fourier öznitelikleri** [23]: ağın yüksek frekanslı salınımları öğrenmesini kolaylaştıran ek girdilerdir.
- **Sert kısıtlar:** başlangıç koşulu tam olarak, duvar koşulu yaklaşık olarak ağın yapısına gömülüdür (README §6;
  `model.py`).

**Eğitim:**

- Her adımda yeni bir 4000 noktalık rastgele örnekle 2000 Adam adımı atılır [37].
- Ardından 200 L-BFGS iterasyonu yapılır [30].

Ayrıntılar Ek C'dedir.

> **Taslak notu (Bölüm 6).** `config.py`; `reference_solver.py` (`cap_profile`); `validate_solver.py`;
> `results/validation/validation_summary.json`, `v*.csv`; `stats_utils.py`; `results/transport/convergence/summary.csv`;
> `docs/PHYSICAL_CAMPAIGN_PLAN.md`; `docs/CAMPAIGN_STATUS_CONTRACT.md`; `model.py`; `pinn_core.py`; `train.py`;
> `README.md` §5–§8; `docs/RGF_B2_VALIDATION.md`.

---

## 7. Sonuçlar

Ortak bilgiler:

- **Kaynak.** Bütün taşınım sonuçları DST çözücüsüyle üretilmiştir (`experiment_disorder.py`). Ham veriler
  `results/transport/<sweep>/<point>.npz`, istatistikler `summary.csv` dosyalarındadır.
- **Aralıklar.** Köşeli parantezler ansambl ortalamasının %95 bootstrap güven aralığıdır (sütunlar `T_ci_boot_lo`,
  `T_ci_boot_hi`).
- **Analiz türü.** Bütün analizler keşifseldir (§6.5).

Her alt bölümde dört tür ifade ayrılmıştır:

- **Gözlem:** simülasyonda ölçülen;
- **Kuramsal beklenti:** Born kuramının öngördüğü;
- **Önerilen açıklama:** akla yatkın ama ayrıca sınanmamış yorum;
- **Sınanmamış:** test edilmemiş mekanizmalar.

### 7.1 Ana sonuç: korelasyon uzunluğu ve iletim

Ayarlar: W = 100, L_d = 1, k_x0 = 20 (k ≈ 20.25), t_end = 1, n = 48 realizasyon.

| ζ | kζ | kℓ_tr (Born) | ⟨T⟩ [%95 CI] | ⟨R⟩ | sd(T) |
|---|---|---|---|---|---|
| 0.02 | 0.40 | 20.8 | 0.670 [0.656, 0.684] | 0.329 | 0.050 |
| 0.03 | 0.61 | 12.4 | 0.545 [0.526, 0.565] | 0.454 | 0.068 |
| 0.05 | 1.01 | 10.4 | 0.494 [0.475, 0.513] | 0.505 | 0.067 |
| 0.07 | 1.42 | 14.4 | 0.535 [0.506, 0.562] | 0.465 | 0.100 |
| 0.10 | 2.02 | 24.0 | 0.590 [0.560, 0.619] | 0.409 | 0.106 |
| 0.15 | 3.04 | 38.9 | 0.694 [0.650, 0.735] | 0.306 | 0.153 |
| 0.20 | 4.05 | 53.0 | 0.734 [0.687, 0.777] | 0.265 | 0.159 |

: Tablo 7.1. Kaynak: `results/transport/corr_length/summary.csv` (`T_mean`, `T_ci_boot_lo/hi`, `R_mean`, `T_std`);
kℓ_tr için `results/transport/analysis.json["corr_length"]` (`k_l_tr`). kζ değerleri k = 20.245 ile hesaplanmıştır
(`results/analysis/paired_corr_length.json`, alan `k`).

![Şekil 7.1. Sabit W = 100'de korelasyon uzunluğunun etkisi (n = 48, %95 bootstrap aralıkları). ⟨T⟩ kζ ≈ 1 civarında en küçüktür; minimum yalnızca yedi noktalı ızgara çözünürlüğünde belirlenmiştir. Kesikli çizgi: Born + difüzyon tahmini (uydurma değildir). Veri: `results/transport/corr_length/summary.csv`.](figures/transport_vs_correlation_length.png)

**Gözlemler.**

1. ⟨T⟩, ζ'nin **monoton olmayan** bir fonksiyonudur. Örneklenen yedi nokta içinde en küçük değer ζ = 0.05'tedir
   (kζ = 1.01). ζ = 0.02, 0.05 ve 0.20 noktalarının %95 aralıkları örtüşmez.
2. Realizasyondan realizasyona dağılım ζ ile büyür: sd(T) 0.050'den 0.159'a çıkar.

**Eşleştirilmiş yeniden analiz (POST-HOC).**

- ζ noktaları aynı tohumları paylaştığı için farklar tohum bazında eşleştirilmiş bootstrap ile yeniden hesaplanmıştır.
- Ayarlar: `paired_analysis.py`, 48 tohum, B = 20000 bootstrap örneği; üç birincil kontrastta Holm çoklu karşılaştırma
  düzeltmesi [38]; bootstrap p-değeri (b+1)/(B+1) biçimindedir [40, 41].
- **Kontrastlar, eşleştirilmemiş sonuçlar görüldükten sonra seçilmiştir.** Bu nedenle bağımsız bir tekrar ya da
  doğrulayıcı bir test değildir.

| Kontrast | Tanım | Ortalama fark | %95 CI (eşleştirilmiş bootstrap) | Holm-düzeltilmiş p |
|---|---|---|---|---|
| C1 | ⟨T⟩(0.05) − ⟨T⟩(0.02) | −0.176 | [−0.198, −0.154] | 1.5×10⁻⁴ |
| C2 | ⟨T⟩(0.05) − ⟨T⟩(0.20) | −0.240 | [−0.284, −0.193] | 1.5×10⁻⁴ |
| C3 | eğrilik: [T(0.03) + T(0.07)]/2 − T(0.05) | +0.046 | [0.023, 0.068] | 1.5×10⁻⁴ |
| C3b | ⟨T⟩(0.07) − ⟨T⟩(0.05) (destekleyici) | +0.041 | [0.006, 0.074]; düzeltilmiş %98.33 CI [−0.002, 0.081] | — |

: Tablo 7.2. Kaynak: `results/analysis/paired_corr_length.json` (`contrasts.*`: `mean`, `boot_ci_0.9500`,
`boot_ci_0.9833`, `p_boot_holm`). Girdi: `results/reproduction/transport/corr_length`; bu dizin
`results/transport/corr_length` ile bit düzeyinde aynıdır (`docs/REPRODUCIBILITY_AUDIT.md` §5).

**Analizin çıktıları:**

- **p değerleri.** Bootstrap çözünürlük sınırındadır (≈ 5×10⁻⁵; Holm sonrası 1.5×10⁻⁴). Bunlar tam p değerleri değildir.
- **Minimumun konumu.** Bootstrap örneklerinin %98.9'unda ζ = 0.05'te, %1.1'inde ζ = 0.07'dedir (`argmin_distribution`).
- **Parabolik tepe kestirimi.** kζ = 0.95 [0.88, 1.12] (`vertex_k_zeta`).
- **Çoklu karşılaştırma.** ζ = 0.07 ile fark, düzeltilmiş aralıkta sıfırı içerir.

**Ana sonucun kapsamı.**

- Bu geometri, bu paket ve W = 100 için ⟨T⟩(ζ)'nin kζ ≈ 1 civarında bir minimumu olduğu veriyle **keşifsel** düzeyde
  desteklenmektedir.
- Minimumun **konumu yalnızca kaba biçimde** (kζ ≈ 0.6–1.4) belirlenmiştir. Bunun dört nedeni vardır:
  1. yalnızca yedi ζ noktası örneklenmiştir;
  2. eşleştirilmiş kontrastlar post-hoc seçilmiştir;
  3. önceden kayıtlı fiziksel kampanya dondurulmamış ve çalıştırılmamıştır;
  4. ζ = 0.07 ile fark çoklu karşılaştırma düzeltmesinden sonra anlamlı değildir.
- Bu davranışın başka W, L_d veya kanal sayılarında sürüp sürmediği **sınanmamıştır.**

### 7.2 Born saçılma kuramıyla fiziksel yorum

**Kuramsal beklenti** (§5.6). Born kuramında ℓ_tr, kζ ≈ 0.8 civarında en küçüktür. Kısa ℓ_tr, daha fazla saçılma ve daha
düşük iletim demektir. Tablo 7.1'deki Born kℓ_tr sütunu, **örneklenen noktalar içinde** en küçük değerini ζ = 0.05'te alır
(10.4).

**Gözlem ile beklentinin karşılaştırılması.**

- Ölçülen ⟨T⟩ minimumu ile Born ℓ_tr minimumu niteliksel olarak aynı yerdedir.
- Bu bir **niceliksel doğrulama değildir:** Born ifadeleri sınırsız 2B ve zayıf saçılma içindir; dalga kılavuzunun kanal
  yapısını içermez.

**Born + difüzyon tahminiyle niceliksel karşılaştırma** (`results/transport/analysis.json`; alanlar `k_l_tr`, `ratio`,
`ratio_avg`):

- Tahmin (6), paketin momentum dağılımı üzerinden ortalanmıştır.
- **Zayıf saçılma.** kℓ_tr ≥ 10 olan bütün noktalarda (ζ, W ve k taramaları) ⟨T⟩_sim/T_est oranı **0.950–1.121**
  aralığındadır. Yani tahmin ölçümü yaklaşık %12 içinde yakalar.
- **ζ taraması.** Yumuşak düzensizlikte (kζ ≳ 2) oran ≈ 0.95–0.97'dir; tahmin T'yi yaklaşık %5 fazla öngörür.
- **Güçlü saçılma.** kℓ_tr azaldıkça oran **sistematik olarak artar:**

  | Nokta | kℓ_tr | oran |
  |---|---|---|
  | W = 200 | 2.6 | 1.213 |
  | k_x0 = 12 | 1.71 | 1.303 |
  | W = 300 | 1.16 | 1.309 |

  Üç taramanın noktaları kℓ_tr'ye karşı ortak bir eğilim izler.

![Şekil 7.2. Simüle edilen ⟨T⟩'nin Born + difüzyon tahminine oranı, kℓ_tr'ye karşı (üç tarama). kℓ_tr ≳ 10'da oran 0.95–1.12; kℓ_tr ≈ 1–3'te 1.2–1.3. Veri: `results/transport/analysis.json`.](figures/born_ratio_vs_kl.png)

[REVIEW: `figures/CAPTIONS.md` bu şekil için "≈ 1.05'e yaklaşır" der; yukarıdaki aralık `analysis.json`'dan alınmıştır
(K-2).]

**Önerilen açıklamalar.**

- **ζ ile artan dağılım.** sd(T)'nin ζ ile büyümesi şöyle açıklanabilir: bölgedeki bağımsız korelasyon alanlarının sayısı
  (∝ L_d L_y/ζ²) azalır ve iletim daha az "öz-ortalanır". Bu ayrıca sınanmamıştır.
- **Güçlü saçılmada artan oran.** Enerji ortalaması bunun bir kısmını açıklar (W = 300'de oran 1.48'den 1.31'e iner),
  tamamını açıklamaz.

**Sınanmamış mekanizmalar.** Güçlü saçılmadaki sapma için aday nedenler şunlardır. **Hiçbiri test edilmemiştir.**

- Born yaklaşımının çökmesi;
- dalga kılavuzunun 6 kanallı ayrık yapısı;
- ℓ/(L+ℓ) biçiminin kabalığı.

![Şekil 7.3. ζ taramasının tohum bazında eşleştirilmiş yeniden analizi (POST-HOC; 48 tohum; B = 20000). Veri: `results/analysis/paired_corr_length.json`.](figures/paired_corr_length.png)

### 7.3 Destekleyici bağlam: düzensizlik gücü W

Ayarlar: ζ = 0.05, L_d = 1, k_x0 = 20, t_end = 1, n = 48. W = 0 için tek bir deterministik koşu yapılmıştır.

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

: Tablo 7.3. Kaynak: `results/transport/strength/summary.csv`; kℓ_tr için `analysis.json["strength"]`.

**Gözlem.**

- ⟨T⟩ incelenen aralıkta W ile monoton ve düzgün biçimde azalır.
- P_in ≤ %0.32'dir, yani ölçüm anında olasılığın neredeyse tamamı bölgeyi terk etmiştir.

**Kuramsal beklenti.** Daha güçlü düzensizlik daha kısa ℓ_tr verir (denklem 5).

**Önceki sürümle ilişki.** Bu sonuç, ilk sürümdeki "düzensizlik iletimi artırır (0.424 ile 0.277)" iddiasının yerini alır.
O iddia bir ölçüm yapaylığıydı (`docs/AUDIT.md` 2.3; README §9.9).

![Şekil 7.4. Düzensiz bölge (L_d = 1, ζ = 0.05) boyunca taşınım; ansambl ortalamaları ve %95 bootstrap aralıkları (n = 48). Kesikli çizgi: Born + difüzyon tahmini. Veri: `results/transport/strength/summary.csv`.](figures/transport_vs_strength.png)

### 7.4 Gelen momentuma bağımlılık

Ayarlar: W = 100, ζ = 0.05, L_d = 1, n = 32; her k için temiz kılavuz (W = 0) referansı da hesaplanmıştır.

| k_x0 | ⟨T⟩ [%95 CI] | ⟨R⟩ | kℓ_tr | temiz kılavuz T / R | t_end |
|---|---|---|---|---|---|
| 12 | 0.251 [0.222, 0.280] | 0.743 | 1.71 | 0.984 / 0.014 | 1.667 |
| 16 | 0.360 [0.333, 0.389] | 0.638 | 4.31 | 0.998 / 0.001 | 1.25 |
| 20 | 0.505 [0.483, 0.529] | 0.494 | 10.4 | 1.000 / 0.000 | 1.0 |
| 25 | 0.691 [0.676, 0.707] | 0.309 | 29.7 | 1.000 / 0.000 | 1.0 |
| 30 | 0.835 [0.826, 0.843] | 0.165 | 77.4 | 0.999 / 0.001 | 1.0 |

: Tablo 7.4. Kaynak: `results/transport/momentum/summary.csv`; kℓ_tr için `analysis.json["momentum"]`.

**Gözlem.** ⟨T⟩ k ile monoton artar; daha hızlı paketler daha az saçılır.

**Kuramsal beklenti.** ℓ_tr, k ile büyür.

**Sistematik hata.** k_x0 = 12'de CAP, temiz kılavuzda bile olasılığın ≈ %1.3–1.4'ünü sahte olarak geri döndürür (temiz
R = 0.014). Bu, o noktada R'yi ≈ 0.01 yukarı kaydırır.

**Ayrıntı.** k_x0 = 20 noktası (n = 32; 0.505), W taramasındaki W = 100 noktasıyla (n = 48; 0.494) aralıklar içinde
uyumludur.

![Şekil 7.5. Gelen dalga sayısının etkisi (W = 100, ζ = 0.05, n = 32). k_x0 = 12'de soğurucu sınırın sahte geri dönüşü R'yi yaklaşık 0.01 yukarı kaydırır. Veri: `results/transport/momentum/summary.csv`.](figures/transport_vs_momentum.png)

### 7.5 Düzensiz bölge uzunluğuna bağımlılık

Ayarlar: ζ = 0.05, k_x0 = 20, n = 24, t_end = 1 + 0.5 L_d.

| L_d | W = 100: ⟨T⟩ [%95 CI] | W = 200: ⟨T⟩ [%95 CI] |
|---|---|---|
| 0.25 | 0.827 [0.801, 0.852] | 0.602 [0.544, 0.656] |
| 0.5 | 0.694 [0.665, 0.722] | 0.414 [0.368, 0.460] |
| 1 | 0.518 [0.490, 0.548] | 0.249 [0.221, 0.278] |
| 2 | 0.320 [0.302, 0.338] | 0.104 [0.092, 0.116] |
| 3 | 0.232 [0.214, 0.250] | 0.056 [0.050, 0.064] |

: Tablo 7.5. Kaynak: `results/transport/length/summary.csv`.

Tablo 7.6'daki iki fit 1/SE² ağırlıklıdır (`analyze_transport.py`):

- **difüzif fit** ℓ/(L + ℓ): her saçılmanın rastgele yürüyüşe katkı yaptığı rejimde beklenen biçim;
- **üstel fit** A e^{−L/ξ}: lokalize rejimde beklenen biçim.

χ²_red (indirgenmiş ki-kare), fitin verilere uyumunu ölçer; 1'e yakın değerler iyi uyum demektir.

| | W = 100 | W = 200 |
|---|---|---|
| Difüzif fit: ℓ, χ²_red | 1.00 ± 0.03, **3.6** | 0.24 ± 0.01, 16.1 |
| Üstel fit: ξ, χ²_red | 2.01 ± 0.06, 6.8 | **1.15 ± 0.04, 5.7** |
| Born ℓ = (π/2)ℓ_tr | 0.81 | 0.20 |
| var(ln T) / ⟨ln T⟩ (L_d = 3) | 0.042 / −1.48 | 0.133 / −2.93 |

: Tablo 7.6. Kaynak: `results/transport/analysis.json["length"]` (`ohmic_l`, `ohmic_l_err`, `ohmic_chi2_red`, `exp_xi`,
`exp_xi_err`, `exp_chi2_red`, `born_l`, `lnT_var`, `lnT_mean`).

**Gözlemler.**

- **W = 100.** Difüzif biçim daha iyi uyar. Uydurulan ℓ, Born değerine yakındır.
- **W = 200.** ⟨T⟩ difüzif biçimden daha hızlı azalır ve üstel biçim daha iyi uyar. ξ ≈ 1.15, 6 kanal için beklenen yarı-1B
  lokalizasyon uzunluğuyla (ξ ~ Nℓ ~ 1.2, mertebe kestirimi [6]) aynı mertebededir.

**Yorum: lokalizasyon kanıtı değildir.** W = 200'deki davranış lokalizasyonun başlangıcıyla **tutarlıdır.** Gerekçeler:

1. Yalnızca L_d/ξ ≤ 2.6'ya ulaşılmıştır; lokalizasyonu kanıtlamak için L_d ≫ ξ gerekir.
2. var(ln T), |⟨ln T⟩|'nin çok altındadır. Oysa derin lokalize rejimde ln T geniş dağılır. Paketin enerji bandı üzerinden
   ortalama alması bu dalgalanmaları bastırır.
3. İki fit de istatistiksel olarak tam değildir (χ²_red > 1). Fitler eşleştirilmiş tasarımın getirdiği korelasyonları
   hesaba katmaz.
4. W ≈ E₀'da düzensiz bölgenin bir kısmı klasik olarak yasaktır (V > E₀). Klasik engellenme de katkıda bulunabilir
   (Ek D).
5. Örnek sayısı daha küçüktür (n = 24) ve bu tarama yeniden üretilmemiştir (§12).

![Şekil 7.6. Düzensiz bölge uzunluğuna bağımlılık (n = 24). W = 200'de üstel biçim daha iyi uyar; bu lokalizasyon kanıtı değildir (L_d/ξ ≤ 2.6, χ²_red > 1). Veri: `results/transport/length/summary.csv`, `analysis.json`.](figures/transport_vs_length.png)

### 7.6 Kapalı kutu sonuçları ve sınırları

**Amaç.** Bir dalga paketi kapalı bir kutuda lokalize olursa, kapladığı alan kutu büyüse bile değişmez. Lokalize olmazsa
kutuyu doldurur. Bu, **katılım oranıyla** (PR) ölçülür; PR olasılık yoğunluğunun etkin olarak kapladığı alandır:

PR = (∫ρ)² / ∫ρ², ρ = |ψ|².

Ayarlar: L × L kutu, L ∈ {1, 2}, düzensizlik her yerde, ζ = 0.05, k_x0 = 20, t_end = 0.5, n = 16 (W = 0 için tek koşu).

| W | geç zaman PR/alan, L = 1 | geç zaman PR/alan, L = 2 | PR(L = 2)/PR(L = 1) |
|---|---|---|---|
| 0 | 0.375 | 0.366 | 3.9 |
| 100 | 0.415 ± 0.007 | 0.425 ± 0.022 | 4.1 |
| 200 | 0.334 ± 0.020 | 0.312 ± 0.031 | 3.7 |
| 300 | 0.267 ± 0.023 | 0.228 ± 0.032 | 3.4 |

: Tablo 7.7. ± = SE. Kaynak: `results/localization/summary.csv` (`PR_late_over_area_mean`, `PR_late_over_area_se`,
`PR_late_mean`).

**Gözlem.** Alan dört katına çıktığında PR de yaklaşık dört katına çıkar (oran 3.4–4.1). Bu boyutlarda **2B lokalizasyon
saptanmamıştır.**

**Kuramsal beklenti.** 2B lokalizasyon uzunluğu ξ ≈ ℓ e^{πkℓ/2}, W ≤ 150 için kutu boyutunu mertebelerce aşar [3]. Bu nedenle
negatif sonuç beklenen sonuçtur.

**Ek açıklamalar.**

- **PR/alan'ın mutlak değeri** (≈ 0.23–0.43) tam yayılmış bir durumda bile 1 değildir. `observables.py` belgesine göre
  rastgele dalga türünden bir durumun ortalama yoğunluğunda PR ≈ alan/2 ile alan/3 arasındadır. Bu yüzden test mutlak değere
  değil, L ile ölçeklemeye dayanır.
- **W = 300'deki düşüş.** Klasik olarak yasak bölgelerin katkısı olabilir. Kuantum lokalizasyonuna atfedilemez (Ek D).
- **Önerilen açıklama.** W = 100'de kutu, temiz kutudan daha düzgün dolmuştur. Bunun nedeni düzensizliğin momentum yönünü
  rastgeleleştirmesi olabilir. Bu ayrıca sınanmamıştır (README §9.6).

**Sınırlar.**

- Kutular en fazla 2 × 2'dir.
- Yalnızca iki boyut ve 16 realizasyon kullanılmıştır.
- Bu çalışma, savunulabilir bir 2B Anderson lokalizasyonu çalışması değildir.

![Şekil 7.7. Kapalı kutularda katılım oranı (L = 1, 2; n = 16). PR alanla orantılıdır; bu boyutlarda 2B lokalizasyon saptanmamıştır. Veri: `results/localization/`.](figures/localization_closed_box.png)

### 7.7 Diğer sonuçlar

Bariyer ile düzensizliğin birlikte etkisi (tünelleme sorusu) Ek C.1'dedir. Bu parametrelerde düzensizlik destekli tünelleme
gözlenmemiştir; düzensizlik iletimi yalnızca azaltmıştır.

**İllüstrasyon şekilleri (nicel analizde kullanılmamıştır):**

| Şekil | Dosya yolu | Not |
|---|---|---|
| Ş7.8 | `figures/disorder_examples.png` | aynı tohumla üç ζ değeri için düzensizlik alanı (W = 100). Deterministik olarak `disorder.py` ile üretilir; ayrı bir veri dosyası yoktur |
| Ş7.9 | `figures/waveguide_snapshots.png` | **illüstrasyon.** Tek bir realizasyon için (W = 150, ζ = 0.05, tohum 0) olasılık yoğunluğu. Saklı verisi yoktur; çözücü `make_figures.py::fig_snapshots` içinde yeniden çalıştırılarak üretilir |

> **Taslak notu (Bölüm 7).**
> - **Taramalar:** `results/transport/{corr_length,strength,momentum,length,convergence}/summary.csv`.
> - **Analizler:** `results/transport/analysis.json`; `results/analysis/paired_corr_length.json`;
>   `results/localization/summary.csv`.
> - **Kod ve belgeler:** `observables.py`; `analyze_transport.py`; `README.md` §8–§9; `figures/CAPTIONS.md`;
>   `make_figures.py`.

---

## 8. İkincil yöntemsel sonuç: PINN değerlendirmesi

### 8.1 Ne yapıldı?

- **Problem.** PINN, 1 × 1 kapalı kutuda kısa bir zaman penceresinde (t ≤ 0.02) bir Gauss paketinin yayılımını çözmüştür.
- **Karşılaştırma.** Her PINN çözümü, **aynı potansiyelle** çalışan klasik çözücünün sonucuyla karşılaştırılmıştır. Klasik
  çözücünün bu problemdeki ayrıklaştırma hatası ≈ 10⁻⁶'dır (V3, V4; Ek A).
- **Koşular.** 18 çalışma koşusu ve 3 teşhis koşusu vardır; donanım yalnızca CPU'dur.
- **Kapsam.** PINN, ana taşınım deneylerinin dalga kılavuzu geometrisinde **çalıştırılmamıştır.** Bu nedenle PINN sonuçları
  ana soruya doğrudan katkı sağlamaz; yöntemsel bir değerlendirmedir.

### 8.2 Klasik çözücüyle karşılaştırma

| Yapılandırma | bağıl L2 hatası ψ(T) | maks \|norm − 1\| |
|---|---|---|
| temiz, 3 ağ tohumu | 0.489, 0.487, 0.567 (ort. 0.514) | ≤ 0.061 |
| W = 28 (E₀/4) | 0.530 | 0.038 |
| W = 56 (E₀/2), 3 ağ tohumu | 0.586, 0.574, 0.610 (ort. 0.590) | ≤ 0.072 |
| W = 112 (E₀) | 0.713 | 0.053 |
| W = 56, ζ = 0.04 / 0.08 / 0.16 | 0.685 / 0.586 / 0.497 | ≤ 0.04 |
| Fourier özniteliği olmadan, temiz | **1.152** | 0.047 |
| zaman penceresi T = 0.01 / 0.02 / 0.04 (temiz) | **0.184** / 0.489 / **0.986** | 0.009 / 0.023 / 0.123 |
| 3 kat eğitim bütçesi (temiz) | 0.281 | 0.007 |

: Tablo 8.1. Bağıl L2 hatası, PINN çözümünün klasik çözümden farkının normunun klasik çözümün normuna oranıdır; 1 hata,
ψ = 0 tahmininin hatasıyla aynı büyüklüktedir. Kaynak: `results/pinn/summary.csv` (`rel_l2_final`, `norm_max_dev`). Tam
tablo Ek C'dedir.

**Sonuçlar:**

- **Doğruluk.** PINN bu bütçede taşınım fiziği için **yeterince doğru değildir.** En kolay ayarda bile hata ≈ %49'dur.
- **Maliyet.** Klasik çözücü 0.7–2.9 s sürerken PINN eğitimi aynı CPU'da 470–3100 s sürmüştür (`reference_time_s`,
  `train_time_s`).
- **Genel değerlendirme.** Bu ileri problem için PINN'in klasik çözücüye göre **bir üstünlüğü bulunmamıştır.**

### 8.3 Gözlenen başarısızlık modları

1. **Zaman penceresi.** Hata pencere uzadıkça büyür: T = 0.01'de %18, T = 0.04'te %99. PINN, paket duvarlara ulaşınca
   oluşan girişim saçaklarını yakalayamaz (README §9.10).
2. **Düzensizlik.**
   - Hata W ile artar: 0.53 → 0.59 → 0.71.
   - ζ küçüldükçe de artar: 0.50 → 0.59 → 0.69.
   - Temiz ile W = 56 farkı (3 + 3 tohum) yalnızca bir eğilimdir. Welch t testi: t = 2.70, p = 0.086 [39]; %5 düzeyinde
     anlamlı değildir.
3. **Spektral yanlılık.** Sinir ağları düşük frekansları daha kolay öğrenir [22]. Fourier öznitelikleri olmadan hata
   1.15'tir.
4. **Kolokasyon noktalarına aşırı uyum.** Sabit eğitim noktalarıyla düşük eğitim kaybına rağmen hata 1.22 olmuştur
   (Ek C.2).
5. **Kısmi iyileşmeler.**
   - 3 kat bütçe hatayı 0.49'dan 0.28'e düşürmüştür.
   - Kalıntı-uyarlamalı örnekleme (RAD) ılımlı bir iyileşme vermiştir; bu sonuç tek tohumludur [28].

Bu modlar literatürde bildirilen PINN güçlükleriyle niteliksel olarak benzerdir [24–26, 29]. Bu benzerliğin nedeni bu
çalışmada sınanmamıştır.

### 8.4 Sınırlar

- Yalnızca CPU bütçesi, tek bir mimari ve yapılandırma başına en fazla 3 ağ tohumu kullanılmıştır.
- Sonuç PyTorch iş parçacığı sayısına bağlıdır (Ek C.2).
- Daha büyük ağlar, daha uzun eğitim, GPU, zaman adımlı eğitim ve dalga kılavuzu geometrisi denenmemiştir.
- Sonuçlar **yalnızca bu bütçe ve bu standart teknikler** içindir; PINN'ler hakkında genel bir hüküm değildir.

![Şekil 8.1. Her PINN koşusunun t = T'deki bağıl L2 hatası, aynı yapılandırma ve aynı düzensizlik realizasyonuyla çalışan klasik çözücüye göre. Veri: `results/pinn/summary.csv`.](figures/pinn_errors.png)

| Şekil | Dosya yolu | Veri kaynağı | Not |
|---|---|---|---|
| Ş8.2 | `figures/pinn_error_vs_time.png` | `results/pinn/<run>/eval_timeseries.csv` | hatanın zamanla büyümesi |
| Ş8.3 | `figures/pinn_vs_reference_W56.png` | `results/pinn/W56_fourier_uniform_s0/` | tek koşu; temsili değil |

> **Taslak notu (Bölüm 8).** `results/pinn/summary.csv`; `results/pinn/*/train_info.json`, `eval_summary.json`,
> `eval_timeseries.csv`; `README.md` §6, §9.10, §10; `docs/REPRODUCIBILITY_AUDIT.md` §5. Welch testi CSV
> değerlerinden `scipy.stats.ttest_ind(equal_var=False)` ile yeniden hesaplanmıştır (R0 aşamasında).

---

## 9. Bağımsız durağan doğrulama hattı: kısa durum

> **Bu bölüm ana fiziksel araştırmadan ayrıdır** ve Bölüm 7'deki hiçbir sayıyı değiştirmez. Ayrıntılar Ek B'dedir.

**Neden var?**

- §6.4'teki A ve B türü testler çözücünün kendi doğruluğunu ölçer.
- Düzensiz taşınımı **kod paylaşmayan** bağımsız bir yöntemle karşılaştırmak için sabit enerjide çalışan bir yöntem
  seçilmiştir. Bu yöntem, mod uzayında özyinelemeli Green fonksiyonu (RGF) ve Fisher–Lee bağıntısıdır [32–34, 36].
- Önceden belirlenmiş B1–B4 ölçütleri tanımlanmıştır (`docs/PHASE2_METHOD_SELECTION.md` §2–§3). Belgenin kuralına göre
  taşınım sonuçları **yalnızca B1–B4'ün hepsi geçerse** bağımsız olarak doğrulanmış sayılır.

| Test | Kısa tanım | Durum |
|---|---|---|
| B1 | temiz kılavuz, 6 konfigürasyon | A–E geçti; **B1-F FAIL** (faz kriteri 10⁻¹⁰, değiştirilmedi) |
| B2 | y'den bağımsız bariyer (düzensizlik yok); RGF ve TDSE, bağımsız 1B denklem çözümüne karşı | **PASS** |
| B3 | sabit düzensizlik realizasyonları; RGF ile seyrek Green referansının karşılaştırması | **FAIL** (dx kriteri) |
| B3X | 6 konfigürasyon; bağımsız dalga fonksiyonu referansı | **FAIL** (dx kriteri) |
| B3Y | başarısız B3X noktalarında ince ızgara, 4 yöntem | **PASS** |
| B4 | sabit düzensizlik realizasyonunda TDSE ile durağan çözümün karşılaştırması | **YAPILMADI** |

: Tablo 9.1. Ölçülen değerler ve kaynaklar: Ek B.1.

**Sonuç.** Bağımsız durağan çapraz doğrulama **tamamlanmamıştır.** Düzensiz taşınım sonuçları durağan bir yöntemle
**henüz bağımsız olarak çapraz doğrulanmamıştır;** bunun için B4 gereklidir.

**RGF sayısal teşhis çalışması (W/Z).**

- B1-F başarısızlığını anlamak için RGF'nin cebirsel olarak eşdeğer bir yeniden yazımı (W/Z) incelenmiştir.
- Bu teşhis çalışmaları **taşınım sonuçlarını doğrulamaz**; yalnızca bir sayısal yöntemin yuvarlama davranışını inceler.

| Konu | Durum |
|---|---|
| W/Z | **NOT VALIDATED** |
| coupled-block deneyi | **DESIGN INSUFFICIENT** (W/Z FAIL anlamına gelmez) |
| ana coupling-scaling deneyi | **NOT RUN** (yalnızca girdi manifestleri üretilmiştir) |
| bağımsız uygulama incelemesi (R0) | **PENDING** |

> **Taslak notu (Bölüm 9).** `docs/PHASE2_METHOD_SELECTION.md` §2–§3; `docs/RGF_B1/B2/B3/B3X/B3Y_VALIDATION.md`;
> `results/rgf_b{1,2,3,3x,3y}_validation/summary.json`; `experiments/wz_rgf_diagnostic/`; `REPORT_READINESS_AUDIT.md` §4.

---

## 10. Tartışma

### 10.1 Ana bulgunun yorumu

Ölçülen ⟨T⟩(ζ) minimumu ile Born kuramının ℓ_tr minimumu niteliksel olarak aynı bölgededir (kζ ~ 1).

| Bulgu | Ne söylenebilir? |
|---|---|
| minimumun varlığı | keşifsel düzeyde desteklenmiştir (Tablo 7.1, 7.2) |
| minimumun konumu | yalnızca kaba biçimde belirlenmiştir (kζ ≈ 0.6–1.4) |
| Born ile niceliksel uyum | zayıf saçılmada ≈ %12; güçlü saçılmada sistematik sapma, nedeni belirlenmemiş |
| başka W, L_d, kanal sayısı | **sınanmamıştır** |

Bu sonuç doğrulayıcı bir sonuç değildir. Doğrulayıcı test önceden kayıtlı kampanya ile yapılabilirdi; kampanya
çalıştırılmamıştır.

### 10.2 Gözlem, beklenti, açıklama ve sınanmamış mekanizma

| Gözlem (simülasyon) | Kuramsal beklenti (Born) | Önerilen açıklama (sınanmamış) | Sınanmamış mekanizma |
|---|---|---|---|
| ⟨T⟩(ζ) minimumu kζ ≈ 1 | ℓ_tr minimumu kζ ≈ 0.8 | — | Born dışı etkilerin minimumu kaydırıp kaydırmadığı |
| sd(T) ζ ile artar | — | daha az bağımsız korelasyon alanı, daha zayıf öz-ortalama | — |
| güçlü saçılmada Born oranı > 1 | Born zayıf saçılma içindir | enerji ortalaması kısmen açıklar | Born'un çökmesi; kanal yapısı; ℓ/(L+ℓ) biçimi |
| ⟨T⟩ W ile azalır | ℓ_tr W ile kısalır | — | — |
| ⟨T⟩ k ile artar | ℓ_tr k ile uzar | — | — |
| W = 200'de üstel benzeri L_d bağımlılığı | yarı-1B lokalizasyon ξ ~ Nℓ | lokalizasyonun başlangıcı (kanıt değil) | klasik engellenmenin katkısı |
| kapalı kutuda PR ∝ alan | 2B ξ kutudan çok büyük | — | — |
| PINN hatası T ve W ile artar | — | girişim saçakları; spektral yanlılık | — |

### 10.3 Bu çalışmanın katkısı (abartısız)

Bu belirli geometri (6 kanallı 2B kılavuz, Gauss-korelasyonlu düzensizlik) için:

- ⟨T⟩'nin ζ'ye bağımlılığının güven aralıklı sayısal haritası ve kζ ≈ 1 civarındaki minimumun keşifsel gösterimi;
- Born + difüzyon tahminine göre kℓ'ye bağlı sistematik sapmanın gözlenmesi;
- W, k ve L_d bağımlılıkları;
- erişilebilir boyutlarda 2B lokalizasyonun saptanmaması;
- PINN'in bu bütçede yetersiz kaldığının nicel gösterimi.

Yeni bir kuram, yeni bir PINN mimarisi, Anderson lokalizasyonunun kanıtı veya deneysel doğrulama **iddia edilmez.**

### 10.4 Açık tutarsızlıklar [REVIEW]

| # | Tutarsızlık | Kaynaklar |
|---|---|---|
| K-1 | README'nin depo yapısı satırı B3Y'den (PASS) söz etmez | `README.md` §12; `docs/RGF_B3Y_VALIDATION.md` |
| K-2 | Born oranı şeklinin açıklaması "≈ 1.05'e yaklaşır" der; veri 0.95–1.12 aralığını verir | `figures/CAPTIONS.md`; `analysis.json` |
| K-3 | Dondurulmuş coupling-scaling spec'inin başlığı "nothing implemented" der; altyapı ve manifestler mevcuttur | `PREDICTIONS_coupling_scaling_v2.md`; `fac5f5e` |
| K-4 | B3Y protokol JSON'unda tarihsel "not run" metni (belgede açıklanmış) | `docs/rgf_b3y_protocol.json` |
| K-6 | Paired analizin girdisi `results/reproduction/...` dizinidir (bit düzeyinde aynı) | `paired_corr_length.json` |
| K-7 | README'deki "kalıntı 0.22" iddiasının kaynağı yok; bu taslakta kullanılmaz | `README.md` §6 |
| K-9 | "Yüksek boyutlu ML" motivasyonu depoda belgelenmemiştir | `git log`; §4.5 |
| K-10 | README'deki "%15 / %25 yasak alan" değerlerinin kaynağı yok; Ek D | `README.md` §9.5–§9.6 |

> **Taslak notu (Bölüm 10).** `README.md` §9.2–§9.11; `REPORT_READINESS_AUDIT.md` §7–§8; R0.1 §10.4.

---

## 11. Sınırlamalar

1. **Geometri ve boyut.**
   - 6 açık kanal vardır ve L_d ≤ 3'tür; kapalı kutular en fazla 2 × 2'dir.
   - Bu boyutlar savunulabilir bir 2B Anderson lokalizasyonu çalışmasına izin vermez.
2. **Enerji ortalaması.** Paket σ_k ≈ 3.5 genişliğinde bir bant taşır. T, R ve dalgalanmaları enerji ortalamalıdır.
3. **Sonlu gözlem süresi.** t_end'de P_in en fazla %1.15'tir (W = 200, L_d = 3; `length/summary.csv`, `P_in_mean`).
4. **Soğurucu sınırlar.** Düşük momentumda sahte geri dönüş vardır: k_x0 = 12'de 1.3×10⁻².
5. **Kuram karşılaştırması.** Born ifadeleri sınırsız 2B ve zayıf saçılma içindir; ℓ/(L+ℓ) bir mertebe tahminidir.
6. **İstatistik.**
   - Nokta başına 24–48 realizasyon vardır.
   - Eşleştirilmiş tasarım noktaları korele kılar; aralıklar nokta başınadır.
   - Analizlerin tümü keşifseldir; ζ kontrastları post-hoc'tur.
   - Yalnızca yedi ζ noktası örneklenmiştir.
7. **Düzensizlik modeli.** Yalnızca Gauss-korelasyonlu düzensizlik incelenmiştir. Tek-nokta istatistiği yaklaşık olarak
   Gauss'tur.
8. **Bağımsız doğrulama.**
   - Düzensiz taşınım için durağan bağımsız çapraz doğrulama (B4) **yapılmamıştır.**
   - Bağımsız yöntemle yapılan tek TDSE karşılaştırması (B2) düzensizlik içermez.
   - RGF hattında B1-F, B3 ve B3X başarısızdır.
9. **Tekrar üretilebilirlik.** Yalnızca ζ ve W taramaları, kapalı kutu, çözücü doğrulaması ve bir PINN koşusu bit
   düzeyinde yeniden üretilmiştir (§12).
10. **PINN.** Yalnızca CPU bütçesi, tek mimari ve en fazla 3 tohum kullanılmıştır. Yalnızca kutu problemi incelenmiş,
    dalga kılavuzu incelenmemiştir.
11. **Modellenmeyenler:** etkileşimler, zamana bağlı düzensizlik, faz bozulması, üç boyut.

---

## 12. Sonuç ve gelecek çalışma

### 12.1 Sonuç

Bu çalışma, 6 kanallı 2B sert duvarlı bir dalga kılavuzunda Gauss-korelasyonlu statik düzensizliğin korelasyon uzunluğunun
dalga paketi iletimine etkisini incelemiştir. Kullanılan DST split-step çözücü iç doğrulama ve yakınsama testlerinden
geçmiştir.

**Ana sonuç.** Sabit W = 100'de ⟨T⟩, ζ'nin monoton olmayan bir fonksiyonudur ve kζ ≈ 1 civarında bir minimum gösterir.
Bu, Born kuramının taşınım ortalama serbest yolu için öngördüğü minimumla niteliksel olarak uyumludur. Kanıt keşifsel ve
post-hoc düzeydedir; minimumun konumu kaba çözünürlüktedir.

**Destekleyici sonuçlar:**

- ⟨T⟩ W ile azalır ve k ile artar.
- Born + difüzyon tahmini zayıf saçılmada ≈ %12 içinde uyar; güçlü saçılmada mekanizması belirlenmemiş bir sapma vardır.
- W = 200'de L_d bağımlılığı üstel benzeridir; bu lokalizasyon kanıtı değildir.
- Kapalı kutularda 2B lokalizasyon saptanmamıştır.

**İkincil yöntemsel sonuç.** PINN, bu bütçede klasik çözücüye yakın doğruluğa ulaşmamıştır.

**Doğrulama durumu.** Düzensiz taşınım sonuçları durağan bir yöntemle henüz bağımsız olarak çapraz doğrulanmamıştır
(**B4 yapılmadı**). RGF/W-Z durumları Bölüm 9'dadır.

### 12.2 Gelecek çalışma

Bunların hiçbiri henüz yapılmamıştır.

1. **B4.** Seçilmiş düzensizlik realizasyonlarında TDSE ile durağan RGF çözümünün enerji genişliği hesaba katılarak
   karşılaştırılması. Bağımsız çapraz doğrulama için gereklidir.
2. **Önceden kayıtlı ζ kampanyası.** Fiziksel kampanyanın dondurulması (kök entropi ve C0 tohumunun kaydı) ve
   çalıştırılması. Daha sık bir ζ ızgarası ile minimumun konumunu daraltmak ve doğrulayıcı bir test yapmak bunun amacıdır.
3. **ζ bulgusunun genelliği.** Başka W, L_d ve kanal sayılarında sınanması.
4. **Daha uzun bölgeler ve dar paketler.** L_d ≈ 6–8 ve ≥ 64 realizasyon; ln T dağılımının sınanması.
5. **Kanal sayısı taraması.** ξ ∝ N bağıntısının sınanması.
6. **Daha büyük kapalı kutular.** L = 4.
7. **PINN.** Daha fazla tohum, daha uzun eğitim, zaman adımlı eğitim ve dalga kılavuzu geometrisi.
8. **W/Z.** Bağımsız R0 incelemesi; ardından coupling-scaling deneyi.
9. **Yüksek boyut.** [TO BE COMPLETED: yazarlar hedefi ve yöntemi tanımlamalıdır; mevcut veriyle desteklenmez.]

---

## 13. Tekrar üretilebilirlik ve veri erişimi

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
| momentum, length, barrier, convergence taramaları; diğer 17 PINN koşusu | yeniden çalıştırılmadı (**UNKNOWN**) [TO BE COMPLETED: `--strict` ile yeniden üretim] |

**Meta veri notu.** 44 transport `.npz` dosyasının 22'sinde `config_json` T = 1.0 yazar, gerçek t_end ise 1.25–2.5'tir.
Sayılar bundan etkilenmez (`docs/REPRODUCIBILITY_AUDIT.md` §2).

---

## 14. Kaynaklar

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
39. B. L. Welch, "The generalization of 'Student's' problem when several different population variances are involved," *Biometrika* **34**(1/2), 28–35 (1947). [REFERENCE DETAILS TO VERIFY: DOI doğrulanmamıştır; doğrulanana kadar DOI verilmemiştir.]
40. B. V. North, D. Curtis, P. C. Sham, "A note on the calculation of empirical P values from Monte Carlo procedures," *Am. J. Hum. Genet.* **71**, 439 (2002). PMID 12111669
41. B. Phipson, G. K. Smyth, "Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn," *Stat. Appl. Genet. Mol. Biol.* **9**(1), Article 39 (2010). doi:10.2202/1544-6115.1585

**Kaynakların doğrulanma yolları.**

- [1]–[31]: `docs/LITERATURE.md`'den alınmıştır.
- [32]–[36], [40], [41]: `docs/PHASE2_METHOD_SELECTION.md` §8, §8.1'den alınmıştır.
- [37]–[39]: R0.1 sırasında yayın kayıtları üzerinden doğrulanmıştır; [39]'un DOI'si doğrulanmamıştır.

**Metindeki kullanım yerleri:**

| Kaynak | Kullanıldığı yer |
|---|---|
| [1]–[5], [8] | §4.1 |
| [2], [3], [6] | §4.2 |
| [9], [13] | §4.1 |
| [10]–[12] | §4.1 |
| [14] | §5.4 |
| [15]–[18] | §6.3 |
| [19] | §6.2 |
| [20], [21], [23], [30], [37] | §6.6 |
| [22], [24]–[26], [28], [29], [39] | §8.3 |
| [27] | Ek C.2 |
| [31], [35] | §6.5 |
| [32]–[34], [36] | §9 |
| [38], [40], [41] | §7.1 |
| [7] | Ek D (2B taşınım için arka plan kaynağı) |

---

## Ek A. Çözücü doğrulamasının ayrıntıları

Kaynaklar: `results/validation/validation_summary.json`, `v1_eigenmode.csv` … `v7_n_modes.csv`; `validate_solver.py`;
`figures/validation_solver.png`.

**A.1 İç çözücü doğrulaması (tür A)**

| Kontrol | Ayar | Sonuç |
|---|---|---|
| V1 tam özmod (sabit V = 50) | 63², Δt = 10⁻⁴, t ≤ 0.1 | max \|hata\| = 4.58×10⁻¹⁴ (`v1_eigenmode_max_abs_err`) |
| V2 analitik serbest Gauss | 127² ve 255² | bağıl L2 ≤ 9.39×10⁻¹⁰ (`v2_free_gaussian_max_rel_l2`) |
| V6 korunum, kapalı kutu | PINN problemi, W = 112.5 (= E₀), CAP yok | max \|N(t) − 1\| = 1.51×10⁻¹² (`v6.closed_norm_max_dev`) |
| V6 korunum, açık kılavuz | W = 200, t ≤ 1 | max \|N + soğurulan_sol + soğurulan_sağ − N(0)\| = 1.60×10⁻¹³ (`v6.open_bookkeeping_max_dev`) |

V6'nın tanımı `validate_solver.py::v6_conservation` fonksiyonundadır. Değerler `results/validation_log.txt` dosyasında da
kayıtlıdır.

**A.2 Yakınsama ve ayrıklaştırma kontrolleri (tür B)**

| Kontrol | Ayar | Sonuç |
|---|---|---|
| V3 Δt, PINN problemi (W = E₀) | 127², Δt = 1.25×10⁻⁶'ya göre | 2. derece; Δt = 10⁻⁵'te bağıl L2 5.29×10⁻⁷ (`v3`, `pinn_W112`) |
| V3 Δt, dalga kılavuzu (W = 200) | dx = 1/64, t = 1 | 2. derece; Δt = 2×10⁻⁴'te \|ΔT\| = 1.96×10⁻⁵, \|ΔR\| = 2.25×10⁻⁵ (`v3`, `waveguide_W200`) |
| V4 ızgara, dalga kılavuzu (W = 200; ζ = 0.02 ve 0.05) | dx = 1/32, 1/64, 1/128 | dx = 1/64'te \|ΔT\| ≤ 1.4×10⁻⁵; dx = 1/32'de ≤ 2×10⁻³ (README §7; ham T değerleri `v4`) |
| V4 ızgara, PINN problemi | dx = 1/32 … 1/128 | dx = 1/128'de (127²) bağıl L2 1.90×10⁻⁶ (`v4_space_convergence.csv`, `pinn_W112_psi`; değer `T` sütununda) |
| V5 CAP sahte geri dönüşü | temiz kılavuz, k_x0 = 10 … 30 | 20: 1.0×10⁻⁴; 25: 1.6×10⁻⁴; 30: 5.19×10⁻⁴; 16: 1.3×10⁻³; **12: 1.33×10⁻²; 10: 3.07×10⁻²** (`v5_cap_reflection.csv`, `spurious_return`) |
| V7 düzensizlik mod sayısı | W = 100, 24 realizasyon | M = 64: 0.488 ± 0.016; 256: 0.515 ± 0.018; 512: 0.518 ± 0.015; 2048: 0.501 ± 0.014 (`v7_n_modes.csv`) |

**A.3 Bağımsız durağan çapraz doğrulama (tür C)**

| Karşılaştırma | Kapsam | Sonuç |
|---|---|---|
| B2: TDSE paket iletimi ile enerji-ortalamalı bağımsız 1B denklem çözümü | düzensizlik yok; y'den bağımsız yumuşak bariyer, V_b = 150 ve 250 | uzaklık 3.2×10⁻⁴ / 7.0×10⁻⁴; kriter 10⁻³; geçti (`docs/RGF_B2_VALIDATION.md`) |
| B4: TDSE ile durağan RGF, sabit düzensizlik realizasyonunda | Bölüm 7'deki düzensiz taşınım | **YAPILMADI** |

**README ile küçük fark.** README, CAP geri dönüşünü k_x0 ≥ 20 için "≤ 5×10⁻⁴" olarak verir. CSV'deki en büyük değer k_x0 = 30'da
5.19×10⁻⁴'tür; bu taslak CSV değerini kullanır.

---

## Ek B. RGF doğrulama hattı ve W/Z teşhis çalışmalarının ayrıntıları

> Buradaki sonuçlar **fiziksel taşınım sonuçlarını doğrulamaz** ve Bölüm 7'yi değiştirmez.

### B.1 B1–B4 ölçümleri

| Test | Amaç | Durum | Ölçülen | Kaynak |
|---|---|---|---|---|
| B1 | temiz kılavuz, 6 konfigürasyon, 50 enerji; analitik kafes cevabı | A–E geçti; **B1-F FAIL** (faz kriteri 10⁻¹⁰, değiştirilmedi) | F (N_y = 63, dx = 1/1024, L_d = 3, N = 3483): faz büyüklüğü 1.76851694417934×10⁻¹⁰; 100 basamaklı kafes fazına göre 1.5424928001550615×10⁻¹⁰ | `results/rgf_b1_validation/summary.json`, `phase_audit_hp_summary.json`; `docs/RGF_B1_VALIDATION.md` |
| B2 | y'den bağımsız bariyer; RGF ile bağımsız 1B ODE | **PASS** | dx = 1/512'de maks \|ΔT\| = 3.59×10⁻⁵ (V_b = 150), 5.95×10⁻⁵ (V_b = 250); kriter 2×10⁻⁴ | `results/rgf_b2_validation/summary.json`; `docs/RGF_B2_VALIDATION.md` |
| B3 | sabit düzensizlik realizasyonları; seyrek Green referansı | **FAIL: dx_1_512_vs_1_1024** (10 grubun 9'u geçti) | S uyumu ≤ 4.9×10⁻¹²; dx iyileştirmesinde maks \|ΔT₁\| = 5.97×10⁻⁴ > 2×10⁻⁴ (W = 200) | `results/rgf_b3_validation/summary.json`; `docs/RGF_B3_VALIDATION.md` |
| B3X | 6 konfigürasyon; dalga fonksiyonu referansı | **FAIL: C11_dx_refinement** (14 kriterin 13'ü geçti) | S uyumu ≤ 5.3×10⁻¹²; C11: 9.17×10⁻⁴ > 2×10⁻⁴ | `results/rgf_b3x_validation/summary.json`; `docs/RGF_B3X_VALIDATION.md` |
| B3Y | başarısız B3X noktalarında ince ızgara; 4 yöntem | **PASS** | S'de ≤ 9.1×10⁻¹¹ uyum; dx değişimleri dört yöntemde aynı, yani sonlu-fark kafes probleminin özelliği | `results/rgf_b3y_validation/summary.json`; `docs/RGF_B3Y_VALIDATION.md` Part II |
| B4 | sabit düzensizlik realizasyonunda TDSE ile durağan çözüm | **YAPILMADI** | — | `docs/PHASE2_METHOD_SELECTION.md` §3 |

**Kapsam:**

- B1 bir iç tutarlılık testidir.
- B2 kanal bağlaşımını sınamaz.
- B3Y sürekli limite yakınsamayı göstermez ve tek bir realizasyonla sınırlıdır.

### B.2 W/Z teşhis çalışmaları

**W/Z değişkenleri.** RGF'nin küçük büyüklüklü değişkenlerle yazılmış cebirsel bir eşdeğeridir:

$$ \Delta_j = (E\,I - \mathrm{diag}(E_\perp) - V_j)\,dx^2, \qquad W_j = D_j/t + I. $$

| Çalışma | Durum (aynen) | Başlıca ölçüm | Sınırlama | Kaynak |
|---|---|---|---|---|
| B1-F en kötü nokta ayrıştırması | — | 110 basamaklı tam ayrıştırma: hatanın +1.6108×10⁻¹⁰'u (≈ %104) saklanan köşegenin binary64 temsilinden gelir | tek nokta | `experiments/wz_rgf_diagnostic/REPORT.md` §4 |
| Scalar W/Z (tek nokta, 1×1 kanal) | **INCONCLUSIVE** | faz −1.19×10⁻¹⁵ (orijinal +1.5425×10⁻¹⁰); genlik +4.82×10⁻¹⁴ (orijinal −1.08×10⁻¹⁵) | yalnızca kanal 5 | `REPORT.md`; `ERRATUM_worst_point.md` |
| 63 kanallı genlik yuvarlaması | H_A: SUPPORTED (kanal 5'te bağımsız olarak ayırt edilemedi) | W/Z fazı ≤ 3.24×10⁻¹⁵ (orijinal 1.0×10⁻¹¹–1.54×10⁻¹⁰); W/Z genliği maks 1.93×10⁻¹³ (≈ 3.4 kat); \|P\| − 1 = 1.37×10⁻¹⁷ → 4.757×10⁻¹⁴ | temiz cihaz | `amplitude_roundoff_report.md` |
| Coupled-block | **DESIGN INSUFFICIENT** | 76 execution / 72 konfigürasyon; referans ≤ 2.0×10⁻⁶¹ / 3.2×10⁻⁶¹; uç bloklar ≤ 2.5×10⁻¹⁶; 480/480 mutant; eğim 0.78–0.79 < 0.8 | N ≤ 4; **W/Z FAIL değil** | `results_coupled_blocks/coupled_blocks_report.md` |
| Coupling-scaling v2 | manifestler **GENERATED AND FROZEN** (`fac5f5e`); deney **NOT RUN** | XP 60, HO 120 konfigürasyon | "duplicate set" tanımsız (G-1); R0 imzası yok | `PREDICTIONS_coupling_scaling_v2.md`; `R0_REVIEW_PACKET.md` |
| Uygulama düzeyi R0 | **PENDING** | — | `REVIEW_coupling_scaling_v2.md` yok | `R0_REVIEW_PACKET.md` |
| T1–T6 yazılım testleri | PASS (yalnızca yazılım ve protokol) | 54/507/28/468/3528/36 kontrol, 0 başarısız; pytest 691 passed + 2 xfailed | fiziksel doğrulama değildir | `R0_REVIEW_PACKET.md` I.7 (`1668540`) |

**Özet:** W/Z **NOT VALIDATED**; coupled-block **DESIGN INSUFFICIENT**; coupling-scaling **NOT RUN**; R0 **PENDING**.

**Ek şekiller:**

| Şekil | Dosya yolu |
|---|---|
| B1 faz hataları | `docs/turkce_ozet_figurler/sekil1_b1_faz.png` |
| 63 kanal faz/genlik | `docs/turkce_ozet_figurler/sekil2_63kanal.png` |
| coupled-block uç blok hatası | `docs/turkce_ozet_figurler/sekil3_coupled.png` |
| B2 hata ve yakınsama | `results/rgf_b2_validation/{error_vs_dx,error_vs_energy,order_vs_energy}.png` |

---

## Ek C. Ek sonuçlar ve PINN ayrıntıları

### C.1 Bariyer ve düzensizlik (tünelleme sorusu)

Ayarlar: genişliği 0.1 olan, bölgenin ortasında yumuşak dikdörtgen bariyer; ζ = 0.05, L_d = 1, n = 24, t_end = 1.5. Bariyer
yüksekliği V_b, boylamasına enerji k_x0²/2 = 200'ün altında (150) veya üstündedir (250).

| V_b | W = 0 | W = 50 | W = 100 | W = 200 |
|---|---|---|---|---|
| 150 | 0.708 | 0.587 [0.561, 0.613] | 0.402 [0.365, 0.440] | 0.203 [0.174, 0.234] |
| 250 | 0.355 | 0.315 [0.296, 0.333] | 0.247 [0.218, 0.277] | 0.160 [0.132, 0.190] |

: Tablo C.1. ⟨T⟩ [%95 CI]. Kaynak: `results/transport/barrier/summary.csv`; şekil `figures/transport_barrier.png`.

**Gözlem.** Düzensizlik iletimi yalnızca azaltmıştır. Bu parametrelerde düzensizlik destekli tünelleme gözlenmemiştir. Bu,
iki bariyer yüksekliği ve üç düzensizlik gücüyle sınırlı bir negatif sonuçtur.

### C.2 PINN eğitim ayrıntıları

**Tam sonuç tablosu.**

| Yapılandırma | koşu(lar) | bağıl L2 ψ(T) | bağıl L2 \|ψ\|²(T) | maks \|norm − 1\| |
|---|---|---|---|---|
| temiz, 3 tohum | `clean_fourier_uniform_s{0,1,2}` | 0.489, 0.487, 0.567 | 0.57–0.76 | ≤ 0.061 |
| W = 28 | `W28_fourier_uniform_s0` | 0.530 | 0.607 | 0.038 |
| W = 56, 3 tohum | `W56_fourier_uniform_s{0,1,2}` | 0.586, 0.574, 0.610 | 0.67–0.75 | ≤ 0.072 |
| W = 56, farklı düzensizlik realizasyonu | `W56_d1_fourier_uniform_s0` | 0.582 | 0.760 | 0.038 |
| W = 112 | `W112_fourier_uniform_s0` | 0.713 | 0.957 | 0.053 |
| W = 56, ζ = 0.04 / 0.08 / 0.16 | `W56_zeta0.04…`, `W56_fourier_uniform_s0`, `W56_zeta0.16…` | 0.685 / 0.586 / 0.497 | 0.90 / 0.69 / 0.58 | ≤ 0.04 |
| Fourier yok: temiz / W = 56 | `*_nofourier_uniform_s0` | 1.152 / 0.724 | 1.42 / 0.78 | ≤ 0.05 |
| RAD örnekleme: temiz / W = 56 | `*_fourier_rad_s0` | 0.469 / 0.525 | 0.45 / 0.56 | ≤ 0.025 |
| T = 0.01 / 0.02 / 0.04 | `clean_T0.01…`, `clean_fourier_uniform_s0`, `clean_T0.04…` | 0.184 / 0.489 / 0.986 | 0.15 / 0.57 / 1.02 | 0.009 / 0.023 / 0.123 |
| 3 kat bütçe | `clean_long_fourier_uniform_s0` | 0.281 | 0.322 | 0.007 |

Kaynak: `results/pinn/summary.csv`.

**Ağ ve kısıt ayrıntıları** (`model.py`, `pinn_core.py`, `train.py`):

- Girdiler [0, 1]'e normalize edilir: x/L_x, y/L_y, t/T.
- Fourier öznitelikleri: birim başına 1–6 çevrim aralığında 6 frekans.
- Sert kısıtlar: u = u₀ + A·B·τ·N_u ve v = v₀ + A·B·τ·N_v; burada B = 16x̂(1−x̂)ŷ(1−ŷ) ve τ = t/T.
- Kayıp: ölçeklenmiş kalıntının (T/A)·(f_u, f_v) ortalama karesi, artı normalizasyon cezası (∫|ψ|² − 1)².

**Optimizasyon:**

- Adam [37]: öğrenme oranı 10⁻³, her 1000 adımda yarıya iner; PDE ağırlığı 300 adımda ısıtılır; gradyan kırpma uygulanır.
- Ardından sabit 16000 noktada L-BFGS [30].

**Teşhis koşuları:**

- **Sabit kolokasyon ile L-BFGS.** Son eğitim kaybı 9.07×10⁻⁴'tür (`diag_fixed_points_lbfgs/train_info.json`,
  `final_loss`). Buna karşın çözüm hatası 1.218'dir (`eval_summary.json`). Bu aşırı uyumdur. Bu koşu için ayrılmış noktalarda
  (*held-out*) bir kalıntı değeri saklanmamıştır.
- **Her adımda taze örnek.** Eğitim ve ayrılmış nokta kalıntıları 5.67×10⁻³ ve 9.28×10⁻³'tür
  (`diag_stochastic/train_info.json`, `final_train_pde`, `holdout_pde`). Çalışma koşularında bu protokol kullanılmıştır.
- **Nedensel ağırlıklandırma** (ε = 1) [27]. Son hata 0.565'tir; ağırlıklandırma olmadan 0.4885'tir (`diag_stochastic_causal`,
  `diag_stochastic`).

**İş parçacığı bağımlılığı.** 2 iş parçacığıyla ağırlıklar bit düzeyinde yeniden üretilmiştir. 4 iş parçacığıyla son hata
0.48850'den 0.48987'ye değişmiştir (`docs/REPRODUCIBILITY_AUDIT.md` §5).

---

## Ek D. Klasik olarak yasak alan oranı üzerine kısa not

**Bu bir ölçüm değil, analitik bir kestirimdir.** Ana soru için gerekli değildir. Yalnızca §7.5 ve §7.6'daki "klasik
engellenme de katkıda bulunabilir" uyarısının mertebesini göstermek için verilmiştir.

**Kestirim.** README §9.5–§9.6, düzensiz bölgede potansiyelin paket enerjisini aştığı (V > E₀) alan oranını W = 200 için
"≈ %15", W = 300 için "≈ %25" olarak verir. README bu değerlerin kaynağını belirtmez.

- Tek-nokta dağılımı Gauss kabul edilirse kestirim P(V > E₀) = ½ erfc(E₀/(√2 W)) olur.
- E₀ = 204.93 ile sonuçlar W = 200 için 0.153, W = 300 için 0.247'dir; README değerleriyle örtüşür. Bu taslak için yalnızca
  formülle aritmetik olarak hesaplanmıştır.

**Bunlar ölçülmüş alan oranları değildir:**

- tek-nokta istatistiği yalnızca yaklaşık olarak Gauss'tur (§5.4);
- S(x) penceresinin kenar etkisi hesaba katılmamıştır;
- realizasyonlar üzerinde bir ölçüm yapılmamıştır.

[REVIEW: README'nin bu değerleri nasıl elde ettiği belgelenmelidir.] Yarı-klasik ve dalga taşınımı arka planı için [7]'ye
bakılabilir.

---

## Ek E. Kısaltmalar

| Kısaltma | Açılımı |
|---|---|
| TDSE | zamana bağlı Schrödinger denklemi (*time-dependent Schrödinger equation*) |
| DST | ayrık sinüs dönüşümü (*discrete sine transform*) |
| CAP | karmaşık soğurucu potansiyel (*complex absorbing potential*) |
| PINN | fizik-bilgili sinir ağı (*physics-informed neural network*) |
| RGF | özyinelemeli Green fonksiyonu (*recursive Green's function*) |
| PR | katılım oranı (*participation ratio*) |
| CI | güven aralığı (*confidence interval*) |
| SE | standart hata |
| sd | standart sapma |
| RAD | kalıntı-uyarlamalı örnekleme (*residual-based adaptive sampling*) |
| L-BFGS | sınırlı bellekli Broyden–Fletcher–Goldfarb–Shanno optimizasyon yöntemi |
| W/Z | RGF'nin küçük büyüklüklü değişkenlerle yazılmış biçimi (Ek B.2) |
| R0 | coupling-scaling deneyi için bağımsız uygulama incelemesi |

---

## Ek F. R0.1 → R0.2 revizyon kaydı (iç not)

### F.1 Yapısal değişiklikler

- **Sonuçların sırası** ana soru etrafında yeniden düzenlendi:
  1. ζ (§7.1);
  2. Born yorumu (§7.2);
  3. W (§7.3);
  4. k (§7.4);
  5. L_d (§7.5);
  6. kapalı kutu (§7.6);
  7. PINN, ayrı bir ikincil bölüm olarak (§8).
- **Bölüm 6'nın yeni düzeni:**
  - doğrulama, üç türü özetleyen tek bir tabloya indirildi (§6.4); ayrıntılı tablolar Ek A'ya taşındı;
  - istatistiksel yöntem ve realizasyon yakınsaması tek bir alt bölümde toplandı (§6.5);
  - PINN yöntemi özetlendi (§6.6); ayrıntılar Ek C.2'ye taşındı.
- **Ek bölümlere taşınanlar:**
  - bariyer/tünelleme sonucu → Ek C.1;
  - "%15 / %25" kestirimleri → Ek D, ana metinden çıkarıldı;
  - RGF/W-Z ayrıntıları → Ek B;
  - kısaltmalar → Ek E (yeni).
- **Kaynakça.** Kaynakça korundu (41 kaynak). Metin içi kullanım tablosu eklendi; daha önce metinde atıf almayan [4], [5],
  [7], [8], [10]–[13], [18], [21], [24], [25] uygun yerlerde atıf aldı.

### F.2 Yeniden ifade edilen iddialar

| R0.1 ifadesi | R0.2 ifadesi |
|---|---|
| "doğrulanmış bir klasik çözücüyle" (§4.3, §4.4, §4.5, §12.1) | "iç doğrulama ve yakınsama testlerinden geçmiş çözücüyle" / "klasik çözücüyle" |
| "Çözücü … doğrulanmıştır" (özet) | üç doğrulama türü açıkça ayrıldı; B4 yapılmadığı açıkça belirtildi |
| "kζ ≈ 1 civarında bir minimumu olduğu verilerle desteklenmektedir" | "keşifsel düzeyde desteklenmektedir; konumu yalnızca kaba biçimde belirlenmiştir"; dört gerekçe listelendi |
| Born ℓ_tr minimumu ile ölçülen minimumun karşılaştırması | "niteliksel uyum; niceliksel doğrulama değil" |
| %15 / %25 yasak alan (ana metinde) | Ek D'de "analitik kestirim, ölçüm değil" |
| "W/Z teşhisi…" | "teşhis çalışmaları taşınım sonuçlarını doğrulamaz" |

> **Taslak notu (Ek F).** Kaynak: `RESEARCH_REPORT_DRAFT_R0_1.md` (`509011c`). Yeni sayısal değer eklenmemiştir.
