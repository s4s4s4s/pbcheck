> **SUPERSEDED WORKING DOCUMENT — NOT A LIVE REQUEST FOR A DECISION, AND NOT THE PRE-REGISTRATION.**
>
> This is the proposal the Phase 0 stratum-list freeze was decided from, committed for auditability
> and nothing else. The binding act is [`docs/PREREGISTRATION_STRATUM_LIST.md`](../../docs/PREREGISTRATION_STRATUM_LIST.md),
> **which governs wherever the two disagree**, and the freeze has already happened: the decision this
> document asks for was taken on 2026-08-16 and is recorded there.
>
> **It has known errors.** They are enumerated, with the correct values, in §10 of the governing
> document under "The discrepancies established are". Do not quote a number from this file.
>
> **Redaction, and its exact extent.** The circulated copy carried an absolute filesystem path from
> the author's Windows account in six places. Every occurrence of the repository root — one in the
> backslash form and five in the forward-slash form, all of them the same path — is replaced here by
> the literal `<REPO>`. **Nothing else is altered**: no prose, no number, no code block, no line
> break. The sha256 of the circulated copy is recorded in §2 of the governing document beside the
> sha256 of this redacted one, so the substitution is checkable by anyone holding the original.
>
> Circulated copy: sha256 `50872414b0727c129a824b0c65ed179674ac5d6c9ecaac53327568b3eae6fb48`, 92589 bytes.

---

# pbcheck Phase 0 — ПРОЕКТ пререгистрируемого стратум-листа

**Статус: ПРЕДЛОЖЕНИЕ. Не коммит, не правка репозитория, не допуск в sweep.**
Решение о заморозке принимает Александр. Репозиторий `<REPO>` при подготовке этого документа только читался.

| | |
|---|---|
| Дата подготовки | 2026-08-16 |
| Манифест-источник | `census_candidates_full.json`, `generated_utc = 2026-08-15T22:18:37+00:00` |
| Census | `2025-01-30` (пин из §1; алиасы `latest`/`stable` отвергаются `open_census()`) |
| Спека | `docs/PHASE0_SPEC.md` §1, §8(d) |
| Амендменты | `docs/AMENDMENTS.md`, Amendment 3 Change 1 (operating envelope) |
| Метаданные датасетов | CELLxGENE Discover curation API, `GET /curation/v1/datasets`, выгружено 2026-08-16 |
| Предложено датасетов | **12** (+2 «сиблинга» из тех же коллекций) |
| Запасных | **5** |
| Mathys 2019 (§8(d)) | **НЕ НАЙДЕН** — ни в кандидатах, ни во всём индексе Discover. См. раздел 5. |

---

## 1. Сводка манифеста

### 1.1 Что в артефакте

| Величина | Значение |
|---|---|
| Всего строк (контрастов `dataset_id × cell_type × disease`) | 2190 |
| `gate_status == 'candidate'` | **1197** |
| `excluded_inclusion_gate` | 981 |
| `excluded_confound` | 12 |
| Уникальных датасетов среди кандидатов | **68** |
| `admitted_to_sweep == True` | **0** (у всех 2190 строк) |

Все 68 датасетов-кандидатов найдены в текущем индексе Discover — «мёртвых»/tombstoned идентификаторов нет, метаданные восстанавливаются полностью.

### 1.2 Почему допуск равен нулю

`admission_blockers` присутствуют у **всех** 2190 строк, одинаковым набором из четырёх:

```
integer_check        2190   — §1 п.4, свойство X, а не obs; считается в io_counts.py
frozen_universe_size 2190   — §1 п.5 / C5; считается в gene_universe.py
sigma_donor_estimate 2190   — OPEN (Amendment 3)
envelope_membership  2190   — производная от sigma_donor_estimate
```

Заголовок манифеста фиксирует это буквально:

> `operating_envelope_source`: "docs/AMENDMENTS.md Amendment 3 Change 1. **SYNTHETIC**: sigma_donor is an unanchored simulator knob, so envelope membership is **PENDING for every row here** and no row may be admitted to the sweep on the strength of the envelope entry alone."

То есть манифест сам себя объявляет **списком предложений, а не списком допущенных**. Настоящий документ — предложение поверх предложения.

### 1.3 Кто отсеялся на inclusion gate (981 строка)

Топ-причины (одна строка может нести несколько):

| Причина | Строк |
|---|---|
| В группе A ровно 1 донор — донор коллинеарен условию | 268 |
| В группе A 1 донор после thin-donor drop (< 3) | 268 |
| В группе A 0 доноров после thin-donor drop | 233 |
| В группе A 2 донора после thin-donor drop | 189 |
| `donor_id` не вложен в condition (доноры несут клетки в обеих группах) | 279 (сумма по вариантам 1–15 доноров) |
| В группе B < 3 доноров после drop | 220 (сумма по вариантам 0/1/2) |
| Ни один донор не проходит правило ≥ 10 клеток/донор | 38 |
| `donor_id` константен (1 донор на весь стратум) | 22 |

Это и есть D4-материал (excluded-fraction ≈ 45% всех контрастов), его надо отчитывать вместе с GO.

### 1.4 Форма пространства кандидатов

**По донорам** (`min(n_donors_A, n_donors_B)` — то, что решает членство в envelope):

| Тир | Условие | Строк-кандидатов | Доля |
|---|---|---|---|
| σ_donor ≤ 0.2 недостижим | ровно 3 | 180 | 15% |
| σ_donor = 0.2 (≥ 4) | 4–7 | 463 | 39% |
| σ_donor = 0.35 (≥ 8) | 8–12 | 243 | 20% |
| σ_donor = 0.5 (≥ 13) | ≥ 13 | 311 | 26% |

Ключевое следствие: **если реальный σ_donor окажется ≈ 0.5, 74% кандидатных стратумов выпадают; если ≈ 0.7 (нужно ≥ 23/группу) — выпадают почти все.** Это ровно тот исход, который Amendment 3 называет «живым результатом исследования, а не сбоем».

**По ассеям** (число датасетов-кандидатов, где ассей встречается):

`10x 3' v3` 48 · `10x 3' v2` 27 · `10x 5' v1` 12 · `10x 3' v1` 6 · `10x 5' v2` 4 · `10x 3' transcription profiling` 2 · `10x 5' transcription profiling` 2 · `10x multiome` 2 · `Drop-seq` 2 · `Seq-Well` 2 · `Seq-Well S3` 2 · `MARS-seq` 2 · `TruDrop` 1

Обе ветки 3′ и 5′ представлены с запасом — требование §1 (iii) выполнимо.

**По суспензии:** `cell` — 50 датасетов, `nucleus` — 19 (один датасет несёт обе).

**По клеткам на донора** (медианы по группам, все 2394 групповых значения): min 10 · p25 36 · медиана 96.5 · p75 304 · max 12 241. Разброс — три порядка; ось D1 закрывается легко.

**Residual df:** min 4 · медиана 24 · max 259.
**Пермутаций < 1000** (полное перечисление обязательно, §4): 233 строки.

### 1.5 Конфаунд-флаги: что стоит знать до чтения карточек

| Флаг | Строк-кандидатов |
|---|---|
| `pooled: unresolved — no pool/library id in obs` | **1197 / 1197** |
| `sequencing_depth_bin` (любой уровень) | 1179 |
| `assay` | 547 |
| `tissue_general` | 299 |
| `suspension_type` | 173 |

**Флаг `pooled` стоит на 100% кандидатов** — Census obs в этом пине не выдаёт pool/library id. По D3 это означает, что донор-псевдобалк для **каждого** предложенного стратума остаётся *нижней оценкой* правильной единицы репликации, и золотой стандарт «донор-псевдобалк калиброван» ни на одном из них не может быть заявлен без оговорки. Это свойство пина Census, а не выбора датасетов, и никакой выбор его не чинит.

---

## 2. Методика отбора

### 2.1 Точные цитаты §1, под которые собирался список

**Единица анализа** (`PHASE0_SPEC.md:62`):

> **Unit of analysis.** The `(dataset_id × cell_type)` stratum, never a pooled population. Group A = one disease term; Group B = `normal`. If a dataset has > 2 disease terms, run one binary `disease-vs-normal` contrast per term; never pool into "any disease".

**Inclusion gate** (`PHASE0_SPEC.md:64–70`):

> **Inclusion gate (ALL must hold per stratum):**
> 1. ≥ 3 distinct `donor_id` in EACH group after filtering.
> 2. ≥ 10 cells per donor in that cell_type (donors below threshold are dropped, not merged).
> 3. `donor_id` present, non-constant, and nested within condition — reject designs where a "donor" appears in only one condition by construction with n=1 (donor collinear with condition).
> 4. Raw integer counts confirmed at load (`io_counts.py` assertion; non-integer strata dropped — Report 3/Census raw not guaranteed integer).
> 5. Frozen universe size ≥ minimum gene count (C5; default 200 genes) — else SKIP.

**Конфаунд-прескрин** (`PHASE0_SPEC.md:71`):

> **Confound pre-screen (`census_select.py`):** Cramér's V + perfect-separation check between condition and each of {assay, suspension_type, tissue_general, sequencing-depth bin, library/pool id}. If condition is perfectly separated by assay/suspension/pool (V ≈ 1) → **EXCLUDE** (inflation uninterpretable). If partially confounded → retain, carry the covariate into the DESeq2 design *only if* C4's df rule allows, and tag. Log the excluded fraction and its characteristics (D4).

**Флаг пулинга** (`PHASE0_SPEC.md:73`):

> **Pooling flag (D3).** If a stratum's donors share a pool/library id, tag `pooled=True`; such strata are usable for the permutation-null floor but excluded from the donor-pseudobulk-is-calibrated gold-standard claim.

**Собственно требование к первому проходу** (`PHASE0_SPEC.md:75`) — центральная цитата этого документа:

> **First pass = 8–12 datasets chosen to SPAN the outcome space** (not cherry-pick wins): (i) 2–3 with a biologically strong expected effect (pseudobulk shown non-null), (ii) 2–3 subtle/low-effect, (iii) deliberate variation in assay (10x 3′ vs 5′), tissue, donor count (some exactly 3v3, some ≥ 8v8), and **cells-per-donor** spanning the pre-registered bins (D1). **Pre-register the stratum list before computing any metric.**

**Состав манифеста** (`PHASE0_SPEC.md:77`):

> **Manifest (auditable, one row per stratum):** `dataset_id, cell_type, cell_type_ontology_depth, n_donors_A/B, n_cells, cells_per_donor_by_group, median_counts_per_cell_by_group, C(D,n_A) permutation count, confound_flags, pooled_flag, residual_df, integer_check, frozen_universe_size, fitType`.

### 2.2 Правила, которые я применял поверх §1

1. **Ось «сильный/тонкий» — литературное суждение, не измерение.** Основание для каждого датасета — его собственная публикация (DOI ниже) плюс общее знание о величине эффекта в данной ткани/болезни. Ни один эффект здесь не измерен нами; это критерий *отбора для покрытия*, а не предсказание результата. См. §7.
2. **Приоритет донор-богатым при прочих равных.** Только они переживут envelope при σ_donor = 0.35 и выше, и только они по A1 несут headline-floor и kill-switch (§4 «high-donor strata (≥ 8v8) carry the headline floor»).
3. **Штраф за косые группы.** Дизайны вида 3v68 или 100v10 отмечены как риск: baseline-группа в 3–6 доноров делает контраст заложником одной группы, а `permutation_count` при этом обманчиво велик (он считается от общего D).
4. **Штраф за высокий Cramér's V по `assay`.** §1 исключает только V ≈ 1, но V ≈ 0.55–0.65 (мета-атласы) означает, что ковариату придётся тащить в дизайн, что съедает df и путает интерпретацию инфляции.
5. **Внутридатасетные контрастные пары ценятся выше.** Датасет, где к одной и той же `normal`-группе прикладываются и сильный, и тонкий disease-термин (Yoshida: COVID-19 vs post-COVID-19; Rexach: AD vs PSP vs Pick; COMBAT: COVID-19 vs influenza; KPMP: AKI vs CKD), даёт ось «сила эффекта» при **зафиксированных** ассее, ткани, лаборатории и глубине. Это самая дешёвая из доступных нам форм контроля.
6. **D2 (кластеризация по датасету) учтён на входе.** Два датасета из одной коллекции (SEA-AD DLPFC/MTG, Emphysema immune/non-immune, Kong TI/colon) — это **не** независимые датасеты для решающего правила. Где я беру сиблинга, он помечен как сиблинг и не считается за второй независимый.

### 2.3 Чего методика сознательно НЕ делает

- Не ранжирует по ожидаемому λ_naive и не оптимизирует под GO. §1: «not cherry-pick wins».
- Не использует `envelope_membership` как критерий отбора (он `pending` у всех). Тир envelope в карточках — **описание, а не допуск**.
- Не проверяет `integer_check` и `frozen_universe_size`: оба вычисляются на загрузке X, а не из obs, и в манифесте стоят `pending`. Любой из 12 датасетов может отвалиться на них.

---

## 3. Предложенный список (12 датасетов)

Колонка «Env» = максимальный тир operating envelope, достижимый хотя бы одним стратумом датасета: `0.5` = есть стратумы ≥ 13v13; `0.35` = есть ≥ 8v8; `0.2` = есть ≥ 4v4; `—` = только 3v3.

| # | dataset_id | Публикация | Ассей | Susp. | Ткань | Болезнь vs normal | Стратумов-кандидатов | Env | Роль в покрытии |
|---|---|---|---|---|---|---|---|---|---|
| 1 | `6f7fd0f1` | Gabitto 2024 Nat Neurosci (SEA-AD) | 10x 3′ v3 + multiome | nucleus | DLPFC | dementia | 18 (все ≥13v13) | **0.5** | (i) сильный; донор-максимум 39v44; cells/donor до 6071 |
| 2 | `ac0c6561` | Rexach 2024 Cell | 10x 3′ v2 + v3 | nucleus | BA4 / insula / V1 | AD, PSP, Pick disease | 27 (23 ≥8v8) | 0.35 | (i) сильный; 3 болезни к одной norm-группе |
| 3 | `d8da613f` | Melms 2021 Nature | 10x 3′ v3 | nucleus | lung | COVID-19 (летальный) | 28 (0 ≥8v8) | 0.2 | (i) сильный максимум эффекта; косые группы 20v7 |
| 4 | `2a498ace` | Yoshida 2022 Nature | **10x 5′ v1** | cell | blood (PBMC) | COVID-19 **и** post-COVID-19 disorder | 47 (27 ≥13v13) | **0.5** | (ii) тонкий + (i) сильный в одном датасете; 5′ |
| 5 | `ebc2e1ff` | Ahern 2022 Cell (COMBAT) | **10x 5′ v1** | cell | blood (PBMC) | COVID-19 **и** influenza | 25 (21 ≥8v8) | 0.35 | 5′; сильный/умеренный контраст; A=100 |
| 6 | `f1606894` | Linna-Kuosmanen 2024 Cell Rep Med | 10x 3′ v3 | nucleus | right atrium | atrial fibrillation | 11 (все ≥13v13) | **0.5** | (ii) тонкий; сердце; сбалансированные 25v29 |
| 7 | `d18736c3` | Binvignat 2024 JCI Insight | 10x 3′ v3 | cell | blood | rheumatoid arthritis | 15 (13 ≥13v13) | **0.5** | (ii) тонкий; **экстремум глубины** (≈300–800 counts/cell) |
| 8 | `a12ccb9b` | KPMP kidney atlas v1.5 | 10x 3′ v3 | nucleus | kidney | acute kidney failure **и** CKD | 37 (25 ≥13v13) | **0.5** | почка; острое vs хроническое к одной norm |
| 9 | `19e46756` | Heimlich 2024 Blood Adv | 10x 3′ v3 | cell | blood (PBMC) | clonal hematopoiesis | 7 (0 ≥8v8) | 0.2 | (ii) самый тонкий ожидаемый эффект |
| 10 | `c893ddc3` | Phan 2024 Nat Commun | 10x 3′ v3 | nucleus | caudate / putamen | opiate dependence | 10 (0 ≥8v8) | 0.2 | (ii) тонкий; **экстремум глубины сверху** (до 57k counts/cell) |
| 11 | `8e47ed12` | Elmentaite 2020 Dev Cell | **10x 3′ v2** | cell | ileal mucosa | Crohn disease | 18 (0 ≥8v8) | 0.2 | кишечник; 3′ v2; низкие cells/donor (18–190) |
| 12 | `4b6af54a` | Wang 2023 Immunity (Emphysema Atlas) | 10x 3′ v3 | cell | alveolus of lung | pulmonary emphysema | 8 (**все ровно 3v3**) | **—** | (iii) обязательный 3v3-якорь §1 |
| 12b | `1e5bd3b8` | *сиблинг* #12, та же коллекция | 10x 3′ v3 | cell | alveolus of lung | pulmonary emphysema | 9 (все 3v3) | — | сиблинг, **не** независимый датасет (D2) |
| 1b | `c2876b1b` | *сиблинг* #1, SEA-AD MTG | 10x 3′ v3 + multiome | nucleus | MTG | dementia | 18 (все ≥13v13) | 0.5 | сиблинг, **не** независимый датасет (D2) |

**Независимых датасетов для D2-кластеризации: 12.** Сиблинги 1b и 12b — опционально, как внутриколлекционный контроль воспроизводимости, но в знаменатель «majority of independent datasets» не идут.

### 3.1 Проверка покрытия по §1 (iii)

| Ось | Требование §1 | Как закрыта |
|---|---|---|
| Сильный эффект | 2–3 | **3**: #1 SEA-AD, #2 Rexach, #3 Melms (+ COVID-плечо #4 и #5) |
| Тонкий эффект | 2–3 | **5**: #9 CHIP, #10 opioid, #6 AF, #7 RA, post-COVID-плечо #4 |
| Ассей 3′ | да | #1,2,3,6,7,8,9,10,11,12 |
| Ассей 5′ | да | **#4 (5′ v1), #5 (5′ v1)** |
| Другие ассеи | не требуется | `10x multiome` подмешан в #1; в резерве есть Drop-seq/Seq-Well/MARS-seq |
| Ткань | вариация | мозг (#1,2,10), лёгкое (#3,12), кровь (#4,5,7,9), сердце (#6), почка (#8), кишечник (#11) — **6 систем органов** |
| Суспензия | не в §1 явно, но D-риски требуют | `nucleus` #1,2,3,6,8,10 · `cell` #4,5,7,9,11,12 — **6 / 6** |
| Ровно 3v3 | «some exactly 3v3» | **#12 (+12b)** — все 8+9 стратумов ровно 3v3 |
| ≥ 8v8 | «some ≥ 8v8» | **#1,2,4,5,6,7,8** — 7 датасетов |
| ≥ 13v13 (единственные, кто выживет при σ=0.5) | не в §1, следствие Amendment 3 | **#1,4,6,7,8** — 5 датасетов |
| cells-per-donor | «spanning the pre-registered bins» | от 18 (#11 M-клетки, #12 ciliated) до 6071 (#1 L2/3-6 IT) и 5872 (#10 олигодендроциты) — **более двух порядков**. ⚠️ Сами bins в спеке не определены, см. §8. |
| counts-per-cell | не в §1, но D1-смежная | от ≈ 300 (#7 RA) до ≈ 57 000 (#10 MSN) — **более двух порядков** |

---

## 4. Карточки по датасетам

Во всех таблицах: `A` = группа болезни, `B` = `normal`; `cpd` = медиана клеток на донора (A/B); `cnts` = медиана counts на клетку (A/B); `rdf` = residual_df; `perm` = `permutation_count`. Флаг `pooled: unresolved` опущен — он стоит на всех строках без исключения.

---

### #1 — `6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3` · SEA-AD DLPFC

**Название:** «Whole Taxonomy - DLPFC: Seattle Alzheimer's Disease Atlas (SEA-AD)»
**Коллекция:** SEA-AD: Seattle Alzheimer's Disease Brain Cell Atlas · Gabitto et al. (2024) Nat Neurosci · `10.1038/s41593-024-01774-5`
**Ассей:** 10x 3′ v3 + 10x multiome · **suspension:** nucleus · **ткань:** dorsolateral prefrontal cortex · 1 395 601 клетка
**Контраст:** `dementia` vs `normal` (термин Census — клинический «dementia», не `Alzheimer disease`)

**Рекомендуемые стратумы:**

| cell_type | A | B | n_cells | cpd A/B | cnts A/B | rdf | Флаги |
|---|---|---|---|---|---|---|---|
| microglial cell | 39 | 44 | 40 625 | 516 / 413 | 4 962 / 4 734 | 81 | depth v=0.074 |
| astrocyte of the cerebral cortex | 39 | 44 | 82 936 | 1 008 / 816 | 9 164 / 9 068 | 81 | depth v=0.046 |
| L2/3-6 intratelencephalic projecting glutamatergic neuron | 39 | 44 | 547 665 | 6 071 / 6 672 | 30 421 / 32 772 | 81 | depth v=0.324 |
| oligodendrocyte | 39 | 44 | 136 076 | 1 147 / 1 758 | 6 040 / 5 921 | 81 | depth v=0.133 |
| cerebral cortex endothelial cell | 35 | 40 | 2 437 | 29 / 29 | 6 069 / 5 901 | 73 | depth v=0.038 |

**Закрывает оси:** (i) сильный эффект · донор-максимум всего списка · nucleus · мозг · cells-per-donor от 29 до 6 071 **внутри одного датасета** (это идеальная площадка для кривой D1 при фиксированных ассее и лаборатории).

**Envelope:** σ=0.35 — **все 18 стратумов проходят**. σ=0.5 (≥13) — **все 18 проходят**. σ=0.7 (≥23) — проходят 18 из 18 по группе A (33–39) и B (34–46). Это **единственный датасет в списке, живущий даже в тире σ=0.7**, и потому самый ценный, если σ_donor окажется высоким.

**Риски:**
- Термин `dementia`, а не `Alzheimer disease`: группа A определена клинически (CDR), гистопатологический континуум внутри неё широк → эффект размывается относительно «чистого» AD-дизайна. Ожидаемая сила эффекта: **умеренно-сильная**, не максимальная.
- Ассей смешанный (3′ v3 + multiome), но флага `assay` в конфаунд-скрине **нет** ни на одном стратуме → ассей не разделяет условия. Это хорошо и подтверждено, а не предположено.
- `cell_type_ontology_depth = pending`, а метки очень мелкие (`L2/3-6 IT`, `sncg GABAergic`, `chandelier pvalb`). Против крупнозернистых меток Rexach (`astrocyte`, `glutamatergic neuron`) это прямой D5-конфликт: **пулить headline через эти два датасета нельзя**.
- Огромная глубина (30–55k counts/cell у нейронных классов) — при матчинге по cells-per-donor надо помнить, что глубина тоже не сматчена.

---

### #2 — `ac0c6561-7a48-4185-af6f-af799f699172` · Rexach 2024 cross-dementia

**Название:** «All Cells - snRNA-seq»
**Коллекция:** Cross-dementia human brain snRNA-seq (Rexach et al 2024) · Rexach et al. (2024) Cell · `10.1016/j.cell.2024.08.019`
**Ассей:** 10x 3′ v2 + v3 · **suspension:** nucleus · **ткань:** Brodmann area 4 / insular cortex / primary visual cortex · 432 555 клеток
**Контрасты:** `Alzheimer disease`, `Pick disease`, `progressive supranuclear palsy` — каждый vs общая `normal`-группа (≈ 9–10 доноров). Описание коллекции: «8-11 cases per diagnosis and brain region».

**Рекомендуемые стратумы:**

| Болезнь | cell_type | A | B | n_cells | cpd A/B | cnts A/B | rdf | perm | Флаги |
|---|---|---|---|---|---|---|---|---|---|
| AD | oligodendrocyte | 10 | 10 | 76 181 | 4 058 / 3 351 | 3 276 / 3 246 | 18 | 1.85e5 | assay v=0.056; depth v=0.12 |
| AD | astrocyte | 10 | 10 | 32 391 | 1 598 / 1 416 | 3 492 / 3 049 | 18 | 1.85e5 | assay v=0.056; depth v=0.324 |
| AD | microglial cell | 10 | 9 | 10 817 | 716 / 448 | 2 776 / 2 402 | 17 | 9.24e4 | assay v=0.029; depth v=0.317 |
| PSP | oligodendrocyte | 11 | 10 | 68 607 | 2 852 / 3 351 | 2 796 / 3 246 | 19 | 3.53e5 | assay v=0.056; depth v=0.135 |
| Pick | glutamatergic neuron | 9 | 10 | 57 685 | 4 060 / 2 662 | 4 694 / 4 405 | 17 | 9.24e4 | assay v=0.085; depth v=0.317 |
| AD | pericyte | 9 | 9 | 916 | 53 / 30 | 2 392 / 1 929 | 16 | 4.86e4 | assay v=0.062; depth v=0.272 |

**Закрывает оси:** (i) сильный эффект · nucleus · мозг · **три болезни к одной norm-группе** (ось «сила эффекта» при полностью зафиксированных ассее/лаборатории/ткани) · cells-per-donor 53 → 4 060 внутри датасета.

**Envelope:** σ=0.35 — **23 из 27 стратумов** (≥ 8v8). σ=0.5 — **0 стратумов**. Потолок `min(A,B) = 11`. Датасет полностью выпадает, если σ_donor окажется ≈ 0.5.

**Риски:**
- Три `T cell`-стратума (3v3, 73–117 клеток, `perm = 20/35`) — на грани, стоит включать только как заведомо шумовые.
- Ассей-конфаунд слабый (V ≤ 0.085) — не проблема.
- Тканей три (BA4/insula/V1) внутри одного `dataset_id`; `tissue_general` для мозга схлопывается, так что скрин его не поймал. **Регион — неучтённый батч.** При наличии `residual_df` 16–19 добавить регион ковариатой можно не везде (C4-правило df ≥ 3 выполняется, но степени свободы заметно съедаются).
- Крупнозернистые метки типов (`astrocyte`, `glutamatergic neuron`) — D5-несовместимость с SEA-AD (см. #1).

---

### #3 — `d8da613f-e681-4c69-b463-e94f5e66847f` · Melms 2021 lethal COVID-19 lung

**Название:** «A molecular single-cell lung atlas of lethal COVID-19»
**Публикация:** Melms et al. (2021) Nature · `10.1038/s41586-021-03569-1`
**Ассей:** 10x 3′ v3 · **suspension:** nucleus · **ткань:** lung · 116 313 клеток
**Контраст:** `COVID-19` vs `normal`

**Рекомендуемые стратумы:**

| cell_type | A | B | n_cells | cpd A/B | cnts A/B | rdf | perm | Флаги |
|---|---|---|---|---|---|---|---|---|
| alveolar macrophage | 20 | 7 | 12 511 | 459 / 266 | 1 768 / 2 330 | 25 | 8.88e5 | depth v=0.316 |
| fibroblast | 20 | 7 | 15 973 | 566 / 411 | 1 269 / 1 002 | 25 | 8.88e5 | depth v=0.316 |
| pulmonary alveolar type 2 cell | 20 | 7 | 20 949 | 371 / **1 834** | 1 538 / 1 290 | 25 | 8.88e5 | depth v=0.316 |
| monocyte | 20 | 7 | 7 379 | 272 / 121 | 697 / 765 | 25 | 8.88e5 | depth v=0.316 |
| CD4-positive, alpha-beta T cell | 20 | 7 | 7 586 | 207 / 325 | 642 / 621 | 25 | 8.88e5 | depth v=0.262 |

**Закрывает оси:** (i) **максимум ожидаемого эффекта во всём списке** — аутопсийное лёгкое при летальном COVID-19 · nucleus · лёгкое.

**Envelope:** σ=0.2 — 25 из 28. σ=0.35 — **0 стратумов**. Потолок `min(A,B) = 7`, поскольку контрольная группа во всём датасете — 7 доноров. Датасет выпадает уже при σ_donor = 0.35.

**Риски:**
- **Косые группы 20v7** во всех рекомендованных стратумах. `perm = 8.88e5` выглядит комфортно, но он считается от C(27,20); реальная разрешающая способность ограничена группой B из 7 доноров.
- У `pulmonary alveolar type 2 cell` cells-per-donor различается в 5 раз между группами (371 vs 1 834) — при том, что cells-per-donor это **первичная ось** headline по D1. Такой стратум нельзя читать без матчинга.
- Ткань аутопсийная: post-mortem interval и агональные изменения — неизмеренный систематический сдвиг между группами, не покрываемый ни одной из ковариат скрина.

---

### #4 — `2a498ace-872a-4935-984b-1afa70fd9886` · Yoshida 2022 PBMC ⭐

**Название:** «PBMC»
**Коллекция:** Local and systemic responses to SARS-CoV-2 infection in children and adults · Yoshida et al. (2022) Nature · `10.1038/s41586-021-04345-x`
**Ассей:** **10x 5′ v1** · **suspension:** cell · **ткань:** blood · 422 220 клеток
**Контрасты:** `COVID-19` vs `normal` **и** `post-COVID-19 disorder` vs `normal` — к одной и той же контрольной группе из 35 доноров.

**Рекомендуемые стратумы — плечо «тонкий» (post-COVID-19 disorder):**

| cell_type | A | B | n_cells | cpd A/B | cnts A/B | rdf | Флаги |
|---|---|---|---|---|---|---|---|
| CD4-positive helper T cell | 20 | 35 | 25 328 | 448 / 269 | 3 539 / 3 305 | 53 | depth v=0.446 |
| classical monocyte | 20 | 35 | 37 297 | 490 / 508 | 2 303 / 3 016 | 53 | depth v=0.263 |
| natural killer cell | 18 | 35 | 35 723 | 515 / 575 | 1 916 / 2 379 | 51 | depth v=0.413 |
| naive B cell | 16 | 35 | 26 189 | 229 / 371 | 2 594 / 2 576 | 49 | depth v=0.418 |

**Рекомендуемые стратумы — плечо «сильный» (COVID-19, острый):**

| cell_type | A | B | n_cells | cpd A/B | cnts A/B | rdf | Флаги |
|---|---|---|---|---|---|---|---|
| classical monocyte | 13 | 35 | 28 868 | 234 / 508 | 4 446 / 3 016 | 46 | depth v=0.518 |
| naive thymus-derived CD4-positive T cell | 13 | 35 | 59 253 | 1 269 / 867 | 3 765 / 3 084 | 46 | depth v=0.518 |
| natural killer cell | 13 | 35 | 26 255 | 203 / 575 | 3 019 / 2 379 | 46 | depth v=0.464 |

**Закрывает оси:** **(ii) тонкий И (i) сильный в одном датасете, при одной и той же контрольной группе, одном ассее (5′ v1), одной лаборатории и одном типе суспензии.** Это самый информативный дизайн в предложении: если λ_naive высок на обоих плечах одинаково — это признак того, что инфляция от псевдорепликации, а не от биологии. Плюс: **обязательное 5′-плечо** оси (iii).

**Envelope:** σ=0.35 — 34 из 47. σ=0.5 — **27 из 47** (в т.ч. **все 4 рекомендованных post-COVID-стратума и 0 из рекомендованных COVID-стратумов при строгом ≥13v13**; у COVID-плеча `min = 13` ровно, т.е. на самой границе). Потолок `min(A,B) = 20`.

**Риски:**
- `sequencing_depth_bin` V = 0.26–0.57 — самый высокий среди донор-богатых кандидатов списка; ковариату придётся тащить, df хватает (rdf 46–53).
- Когорта смешана по возрасту (children and adults) — возраст в дизайне не учтён; в obs он есть (`development_stage`) и должен быть проверен на баланс до заморозки.
- `post-COVID-19 disorder` — клинический термин с очень широкой дисперсией фенотипа. Это ровно та причина, по которой я отношу его к «тонким», но она же означает, что нулевой результат здесь ничего не докажет про инструмент.

---

### #5 — `ebc2e1ff-c8f9-466a-acf4-9d291afaf8b3` · COMBAT blood atlas

**Название:** «COMBAT project: single cell gene expression data from COVID-19, sepsis and flu patient PBMCs»
**Коллекция:** A blood atlas of COVID-19 defines hallmarks of disease severity and specificity · Ahern et al. (2022) Cell · `10.1016/j.cell.2022.01.012`
**Ассей:** **10x 5′ v1** · **suspension:** cell · **ткань:** blood · 836 148 клеток
**Контрасты:** `COVID-19` vs `normal` **и** `influenza` vs `normal` (общая norm-группа = 10 доноров).

**Рекомендуемые стратумы:**

| Болезнь | cell_type | A | B | n_cells | cpd A/B | cnts A/B | rdf | Флаги |
|---|---|---|---|---|---|---|---|---|
| COVID-19 | classical monocyte | 100 | 10 | 208 034 | 1 444 / 1 580 | 3 379 / 2 874 | 108 | depth v=0.162 |
| COVID-19 | CD4-positive, alpha-beta T cell | 100 | 10 | 264 564 | 1 873 / 3 092 | 3 730 / 3 805 | 108 | depth v=0.049 |
| COVID-19 | plasmablast | 98 | 10 | 8 148 | 51 / 18 | 23 907 / 22 309 | 106 | depth v=0.045 |
| influenza | non-classical monocyte | 10 | 10 | 5 368 | 66 / 410 | 6 125 / 5 317 | 18 | depth v=0.12 |
| influenza | B cell | 10 | 10 | 3 900 | 39 / 309 | 3 766 / 3 403 | 18 | depth v=0.12 |

**Закрывает оси:** второе **5′**-плечо · кровь · **два инфекционных контраста к одной norm-группе** (COVID-19 — сильный, influenza — умеренный) · экстремум counts/cell сверху у плазмобластов (≈ 23k).

**Envelope:** σ=0.35 — 21 из 25. σ=0.5 — **0**. Потолок `min(A,B) = 10`, задан контрольной группой. Как и Rexach, выпадает при σ_donor = 0.5.

**Риски:**
- **Крайне косые группы 100v10.** `perm = 4.69e13` вводит в заблуждение: реальная разрешающая способность определяется 10 контролями. Это же делает датасет отличной иллюстрацией к A1/риску 3 — стоит явно использовать его как демонстрацию того, что `permutation_count` ≠ информативность.
- У influenza-стратумов cells-per-donor различается в 5–10 раз между группами (66 vs 410; 39 vs 309) — матчинг обязателен.
- Публикация — «COVID-19, sepsis and flu», но `sepsis` в кандидатах не всплыл; проверить, не отфильтрован ли он на inclusion gate, до заморозки.

---

### #6 — `f1606894-59df-4794-a37f-baa7c6fb6de1` · PERIHEART, right atrium

**Название:** «PERIHEART»
**Коллекция:** Transcriptomic and spatial dissection of human ex vivo right atrial tissue... · Linna-Kuosmanen et al. (2024) Cell Reports Medicine · `10.1016/j.xcrm.2024.101556`
**Ассей:** 10x 3′ v3 · **suspension:** nucleus · **ткань:** right atrium auricular region · 392 819 клеток
**Контраст:** `atrial fibrillation` vs `normal`

**Рекомендуемые стратумы:**

| cell_type | A | B | n_cells | cpd A/B | cnts A/B | rdf | Флаги |
|---|---|---|---|---|---|---|---|
| cardiac muscle cell | 25 | 29 | 86 194 | 1 478 / 1 610 | 3 626 / 2 959 | 52 | depth v=0.368 |
| fibroblast | 25 | 29 | 73 117 | 1 065 / 1 512 | 1 932 / 1 758 | 52 | depth v=0.189 |
| endocardial cell | 25 | 29 | 118 154 | 1 707 / 2 478 | 2 412 / 2 054 | 52 | depth v=0.189 |
| macrophage | 23 | 29 | 14 188 | 213 / 212 | 1 773 / 1 755 | 50 | depth v=0.241 |
| smooth muscle cell | 23 | 29 | 4 502 | 83 / 74 | 1 988 / 1 861 | 50 | depth v=0.335 |

**Закрывает оси:** (ii) тонкий/умеренный · **сердце** · nucleus · **самые сбалансированные группы во всём списке (25v29)** при ≥13v13.

**Envelope:** σ=0.35 — **все 11**. σ=0.5 — **все 11**. Потолок `min(A,B) = 25` → в тир σ=0.7 (≥23) проходят стратумы с A=25/B=29 и A=23/B=29. Второй по устойчивости датасет после SEA-AD.

**Риски:**
- Ткань «ex vivo right atrial», забор при кардиохирургии — контроль это не здоровое сердце, а пациент без AF, идущий на операцию. Разница «AF vs не-AF» тоньше, чем «болезнь vs здоровый», что и делает это хорошим тонким контрастом, но одновременно означает, что ожидаемый эффект может быть близок к нулю.
- Сиблинг `8f4f8502` (CAREBANK, та же коллекция, 9–15 vs 12–31) — **не** независимый датасет по D2. Если брать, то как внутриколлекционный контроль.

---

### #7 — `d18736c3-6292-4379-919a-d6d973204c87` · Binvignat 2024 RA blood

**Название/публикация:** «Single-cell RNA-seq analysis reveals cell subsets and gene signatures associated with Rheumatoid Arthritis Disease Activity» · Binvignat et al. (2024) JCI Insight · `10.1172/jci.insight.178499`
**Ассей:** 10x 3′ v3 · **suspension:** cell · **ткань:** blood · 108 717 клеток
**Контраст:** `rheumatoid arthritis` vs `normal`

**Рекомендуемые стратумы:**

| cell_type | A | B | n_cells | cpd A/B | **cnts A/B** | rdf | Флаги |
|---|---|---|---|---|---|---|---|
| classical monocyte | 18 | 18 | 15 391 | 340 / 392 | **548 / 520** | 34 | depth v=0.136 |
| central memory CD4-positive T cell | 18 | 18 | 21 531 | 556 / 574 | **369 / 321** | 34 | depth v=0.272 |
| CD8-positive, alpha-beta memory T cell | 18 | 18 | 2 292 | 57 / 56 | **426 / 670** | 34 | **флагов нет вообще** |
| natural killer cell | 18 | 18 | 9 698 | 175 / 324 | **318 / 327** | 34 | depth v=0.136 |
| naive B cell | 18 | 18 | 6 760 | 104 / 208 | **417 / 346** | 34 | depth v=0.236 |

**Закрывает оси:** (ii) тонкий · кровь · **экстремум глубины снизу: 300–800 counts/cell**, на порядок ниже медианы по кандидатам, и на два порядка ниже нейронных стратумов #1/#10. Это делает датасет прямым тестом на то, не является ли наблюдаемая инфляция артефактом глубины (риск 2 / D1). Плюс — **идеально сбалансированные 18v18**.

**Envelope:** σ=0.35 — 14 из 15. σ=0.5 — **13 из 15**. Потолок `min(A,B) = 18`. σ=0.7 (≥23) — не проходит.

**Риски:**
- Очень низкая глубина сама по себе — риск для `frozen_universe_size` (§1 п.5, ≥ 200 генов). При 300 counts/cell псевдобалк на 57 клеток донора (CD8 memory) даёт ≈ 24k counts на псевдосэмпл; универсум может не набраться. **Этот датасет — наиболее вероятный кандидат на вылет по `frozen_universe_size`**, и это надо проверить до заморозки, а не после.
- Стратум `CD8-positive, alpha-beta memory T cell` не имеет ни одного конфаунд-флага — редкость в этом манифесте, стоит держать его как «чистый» референс.

---

### #8 — `a12ccb9b-4fbe-457d-8590-ac78053259ef` · KPMP kidney snRNA-seq v1.5

**Название:** «Single-nucleus RNA-seq of the Adult Human Kidney (Version 1.5)»
**Коллекция:** An atlas of healthy and injured cell states and niches in the human kidney (Version 1.5) · **DOI в Discover: `null`** — публикация коллекцией не проставлена, *не установлено* из API (для протокола: это KPMP-консорциум, но подтвердить по Discover нельзя)
**Ассей:** 10x 3′ v3 · **suspension:** nucleus · **ткань:** cortex of kidney / kidney / renal medulla / renal papilla · 304 989 клеток
**Контрасты:** `acute kidney failure` vs `normal` **и** `chronic kidney disease` vs `normal`

**Рекомендуемые стратумы:**

| Болезнь | cell_type | A | B | n_cells | cpd A/B | cnts A/B | rdf | Флаги |
|---|---|---|---|---|---|---|---|---|
| AKI | epithelial cell of proximal tubule | 16 | 24 | 53 374 | 1 098 / 998 | 3 355 / 2 646 | 38 | depth v=0.278 |
| AKI | kidney loop of Henle thick ascending limb | 16 | 23 | 29 067 | 312 / 581 | 5 078 / 5 092 | 37 | depth v=0.266 |
| AKI | endothelial cell | 16 | 24 | 12 836 | 238 / 264 | 2 431 / 3 122 | 38 | depth v=0.26 |
| CKD | epithelial cell of proximal tubule | 38 | 24 | 68 118 | 690 / 998 | 2 074 / 2 646 | 60 | depth v=0.208 |
| CKD | kidney interstitial fibroblast | 37 | 24 | 16 911 | 213 / 118 | 2 228 / 3 047 | 59 | depth v=0.167 |
| CKD | podocyte | 32 | 17 | 5 102 | 80 / 141 | 2 824 / 3 552 | 47 | depth v=0.265 |

**Закрывает оси:** **почка** (шестая система органов) · nucleus · **острое vs хроническое повреждение к одной norm-группе** — ещё одна внутридатасетная ось «сила эффекта» · донор-богатый.

**Envelope:** σ=0.35 — 29 из 37. σ=0.5 — **25 из 37**. Потолок `min(A,B) = 24`. σ=0.7 (≥23) — проходят CKD-стратумы с B=24 (`epithelial cell of proximal tubule`, `endothelial cell`, `kidney interstitial fibroblast`). Третий по устойчивости после SEA-AD и PERIHEART.

**Риски:**
- Регион почки (cortex/medulla/papilla) внутри одного `dataset_id` — как и с Rexach, `tissue_general` этого не ловит, регион остаётся неучтённым батчем.
- Сиблинг `dea717d4` («Single-cell RNA-seq of the Adult Human Kidney», **suspension = cell**, 10x 3′ v3, 43 стратума, 26 ≥8v8) из **той же коллекции** — соблазнительная пара «cell vs nucleus при одной ткани и одной когорте». Но по D2 это не второй независимый датасет. Если брать — только как внутриколлекционный контроль суспензии, явно помеченный.
- DOI коллекции отсутствует в Discover → провенанс публикации из API **не установлен**.

---

### #9 — `19e46756-9100-4e01-8b0e-23b557558a4c` · Heimlich 2024 CHIP PBMC

**Название:** «Single Cell Sequencing of Human PBMCs in Clonal Hematopoeisis of Indeterminant Potential»
**Коллекция:** Multiomic Profiling of Human Clonal Hematopoiesis Reveals Genotype and Cell-Specific Inflammatory Pathway Activation · Heimlich et al. (2024) Blood Advances · `10.1182/bloodadvances.2023011445`
**Ассей:** 10x 3′ v3 · **suspension:** cell · **ткань:** blood · 66 985 клеток
**Контраст:** `clonal hematopoiesis` vs `normal`

**Рекомендуемые стратумы:**

| cell_type | A | B | n_cells | cpd A/B | cnts A/B | rdf | perm | Флаги |
|---|---|---|---|---|---|---|---|---|
| CD14-positive monocyte | 17 | 6 | 18 604 | 685 / 568 | 3 332 / 3 571 | 21 | 1.01e5 | depth v=0.214 |
| CD4-positive, alpha-beta T cell | 17 | 7 | 19 046 | 667 / 598 | 1 735 / 2 058 | 22 | 3.46e5 | depth v=0.13 |
| natural killer cell | 17 | 7 | 9 014 | 358 / 330 | 1 603 / 1 771 | 22 | 3.46e5 | depth v=0.259 |
| CD16-positive, CD56-dim natural killer cell | 15 | 6 | 913 | 46 / 28 | 4 060 / 3 422 | 19 | 5.43e4 | **флагов нет** |

**Закрывает оси:** **(ii) самый тонкий ожидаемый эффект во всём предложении.** CHIP — это доклональное состояние без клинического фенотипа; сама публикация говорит об активации воспалительных путей, то есть об эффекте на уровне отдельных программ, а не глобального транскриптомного сдвига. Если псевдобалк и здесь даёт ноль, а naive — сотни генов, это самая чистая демонстрация floor'а на реальных данных.

**Envelope:** σ=0.2 — все 7. σ=0.35 — **0**. Потолок `min(A,B) = 7` (контрольная группа 5–7 доноров). Датасет выпадает уже при σ_donor = 0.35.

**Риски:**
- Косые 17v6/17v7. `perm` большой, разрешающая способность — нет.
- Ожидаемый эффект настолько мал, что стратум может оказаться де-факто «no-effect стратумом» — что по §4/A1 как раз ценно («clean truth=0 guarantee lives in oracle (b) and in **no-effect strata**»), но должно быть заявлено ЗАРАНЕЕ, а не после просмотра результата. **Предлагаю пререгистрировать этот датасет именно как ожидаемый no-effect якорь.**

---

### #10 — `c893ddc3-f25b-45e2-8c9e-155918b4261c` · Phan 2024, opioid use disorder, striatum

**Название/публикация:** «Transcriptional responses of the human dorsal striatum in opioid use disorder implicates cell type-specific programs» · Phan et al. (2024) Nat Commun · `10.1038/s41467-024-45165-7`
**Ассей:** 10x 3′ v3 · **suspension:** nucleus · **ткань:** caudate nucleus / putamen · 98 848 клеток
**Контраст:** `opiate dependence` vs `normal`

**Рекомендуемые стратумы:**

| cell_type | A | B | n_cells | cpd A/B | **cnts A/B** | rdf | perm | Флаги |
|---|---|---|---|---|---|---|---|---|
| oligodendrocyte | 6 | 6 | 62 982 | **5 872 / 3 000** | 6 916 / 7 999 | 10 | 9.24e2 | depth v=0.408 |
| direct pathway medium spiny neuron | 6 | 6 | 5 385 | 520 / 280 | **45 516 / 56 842** | 10 | 9.24e2 | depth v=0.408 |
| indirect pathway medium spiny neuron | 6 | 6 | 6 172 | 602 / 290 | **41 745 / 43 479** | 10 | 9.24e2 | depth v=0.408 |
| astrocyte | 6 | 6 | 8 463 | 684 / 624 | 10 920 / 12 819 | 10 | 9.24e2 | depth v=0.408 |
| microglial cell | 6 | 6 | 6 630 | 463 / 328 | 4 710 / 4 914 | 10 | 9.24e2 | depth v=0.408 |

**Закрывает оси:** (ii) тонкий · мозг · nucleus · **экстремум counts/cell сверху: до 57 000** на клетку у MSN — в 150 раз выше, чем у RA (#7). Пара #7 ↔ #10 задаёт диапазон глубины почти в два с половиной порядка при том, что обе — 10x 3′ v3. Плюс `min(A,B) = 6` при 5 872 клетках на донора — редкая комбинация «мало доноров, очень много клеток», то есть именно тот угол, где псевдорепликация должна бить сильнее всего.

**Envelope:** σ=0.2 — 9 из 10. σ=0.35 — **0**. Потолок `min(A,B) = 6`.

**Риски:**
- `sequencing_depth_bin` V = 0.408 на всех пяти рекомендованных стратумах — глубина частично разделяет условия, а глубина здесь экстремальная. При rdf = 10 ковариату добавить можно (C4: df ≥ 3), но запас невелик.
- 6v6 → `perm = 924` (C(12,6) = 924), то есть **полное перечисление обязательно** (§4: «Small D → enumerate all»), и разрешение нуля грубое. Это надо заявить заранее как ограничение, а не открыть постфактум.
- Opiate dependence — посмертная когорта; токсикология, PMI, причина смерти не сматчены и в obs не выражены.

---

### #11 — `8e47ed12-c658-4252-b126-381df8d52a3d` · Elmentaite 2020 paediatric gut, Crohn

**Название:** «Paediatric Human Gut (4-14y)»
**Коллекция:** Single-Cell Sequencing of Developing Human Gut Reveals Transcriptional Links to Childhood Crohn's Disease · Elmentaite et al. (2020) Developmental Cell · `10.1016/j.devcel.2020.11.010`
**Ассей:** **10x 3′ v2** · **suspension:** cell · **ткань:** ileal mucosa · 22 502 клетки
**Контраст:** `Crohn disease` vs `normal`

**Рекомендуемые стратумы:**

| cell_type | A | B | n_cells | cpd A/B | cnts A/B | rdf | perm | Флаги |
|---|---|---|---|---|---|---|---|---|
| intestinal crypt stem cell | 7 | 8 | 1 190 | 88 / 47 | 10 266 / 10 380 | 13 | 6.44e3 | depth v=0.189 |
| enterocyte | 6 | 8 | 5 060 | 190 / 443 | 5 636 / 7 265 | 12 | 3.00e3 | depth v=0.285 |
| intestine goblet cell | 7 | 8 | 1 049 | 68 / 42 | 5 540 / 5 032 | 13 | 6.44e3 | depth v=0.5 |
| CD4-positive, alpha-beta T cell | 5 | 5 | 1 254 | 110 / 73 | 3 628 / 3 235 | 8 | 2.52e2 | depth v=0.447 |
| memory B cell | 5 | 5 | 1 032 | 79 / 108 | 4 072 / 3 473 | 8 | 2.52e2 | depth v=0.408 |

**Закрывает оси:** **кишечник** · **чистый 10x 3′ v2** (единственный в основном списке датасет на одном только v2) · cell · **низкие cells-per-donor (18–190)** — нижний конец оси D1 · малые группы 5v5–7v8.

**Envelope:** σ=0.2 — 15 из 18. σ=0.35 — **0**. Потолок `min(A,B) = 7`.

**Риски:**
- Педиатрическая когорта 4–14 лет: возраст сильно варьирует и в дизайн не входит.
- `perm` от 252 до 6 440 — для стратумов 5v5 полное перечисление (C(10,5)=252) обязательно, разрешение нуля грубое.
- Флаги глубины до V = 0.655 (`fibroblast`) и `near_confound_v=0.851` у двух малых стратумов — эти два (`activated CD4 T`, `conventional dendritic cell`) я в рекомендацию **не** включаю.
- Всего 22 502 клетки в датасете — самый маленький в списке; риск не набрать `frozen_universe_size` вторичен после #7.

---

### #12 — `4b6af54a-4a21-46e0-bc8d-673c0561a836` · Emphysema Cell Atlas (non-immune) — якорь 3v3

**Название:** «non-immune cells»
**Коллекция:** Emphysema Cell Atlas · Wang et al. (2023) Immunity · `10.1016/j.immuni.2023.01.032`
**Ассей:** 10x 3′ v3 · **suspension:** cell · **ткань:** alveolus of lung · 18 386 клеток
**Контраст:** `pulmonary emphysema` vs `normal` — **все 8 стратумов ровно 3v3**

**Рекомендуемые стратумы:**

| cell_type | A | B | n_cells | cpd A/B | cnts A/B | rdf | perm | Флаги |
|---|---|---|---|---|---|---|---|---|
| fibroblast | 3 | 3 | 5 415 | 1 080 / 627 | 9 573 / 7 730 | 4 | **20** | **флагов нет** |
| endothelial cell | 3 | 3 | 5 874 | 915 / 1 499 | 5 436 / 4 726 | 4 | **20** | **флагов нет** |
| ciliated cell | 3 | 3 | 117 | 18 / 19 | 15 361 / 14 230 | 4 | **20** | **флагов нет** |
| epithelial cell of lower respiratory tract | 3 | 3 | 710 | 33 / 185 | 16 425 / 12 844 | 4 | **20** | **флагов нет** |
| pulmonary alveolar type 2 cell | 3 | 3 | 2 891 | 359 / 716 | 20 316 / 22 559 | 4 | 20 | depth **near_confound v=0.816** |

**Сиблинг `1e5bd3b8` («immune cells», та же коллекция, 9 стратумов, все 3v3)** — стоит взять 2–3 стратума как внутриколлекционный контроль, **не** считая за независимый датасет:

| cell_type | A | B | n_cells | cpd A/B | cnts A/B | rdf | Флаги |
|---|---|---|---|---|---|---|---|
| monocyte | 3 | 3 | 8 823 | 1 113 / 1 576 | 9 883 / 9 299 | 4 | depth near_confound v=0.816 |
| T cell | 3 | 3 | 12 953 | 3 106 / 871 | 3 873 / 3 539 | 4 | depth near_confound v=0.816 |
| alveolar macrophage | 3 | 3 | 680 | 91 / 110 | 1 522 / 6 788 | 4 | **флагов нет** |

**Закрывает оси:** **обязательное требование §1 «some exactly 3v3»** · лёгкое · cell (в пару к nucleus-лёгкому #3 Melms — та же ткань, другая суспензия, другой датасет) · cells-per-donor от 18 до 3 106 внутри пары.

**Envelope:** **вне envelope при любом σ_donor.** Даже при σ = 0.2 требуется ≥ 4 донора на группу; здесь ровно 3. По Amendment 3 **ни один результат отсюда не может быть отчитан как pbcheck-измерение мощности.** Датасет включён потому, что §1 явно требует 3v3, и потому, что по §4/A1 такие стратумы всё равно нужны — но исключительно как иллюстрация того, что при 3v3 sharp-null не ортогонален истине («at 3v3 no balanced permutation is orthogonal to the true grouping»).

**Риски:**
- `perm = 20` — полное перечисление, из них после снятия identity+complement остаётся 10 (§4). Разрешение нулевого распределения предельно грубое; это надо писать в отчёте рядом с каждым числом.
- `residual_df = 4` — на грани C4/C5.
- `near_confound_v = 0.816` по глубине у половины стратумов: при 3v3 добавить ковариату глубины **невозможно** без обнуления df. Эти стратумы придётся отчитывать без коррекции, с явной пометкой.
- Эмфизема против «normal» лёгкого — эффект скорее умеренный, но при 3v3 это уже не важно: датасет здесь не ради эффекта, а ради нижней границы дизайна.

---

## 5. Mathys 2019 / оракул §8(d) — критическая находка

### 5.1 Что требует спека

`PHASE0_SPEC.md:212–213`:

> **(d) Real anchor — Mathys 2019 AD snRNA-seq** (via Census if present, else Synapse syn18485175).
> Run the full stratified pipeline. PASS = qualitatively reproduces Murphy & Skene 2023: naive per-cell DE grossly exceeds pseudobulk and the permutation-null floor accounts for most naive calls (high `lambda_naive`), consistent with pseudoreplication dominance. **This is the BINDING check if the simulators are optimistic.**

Amendment 3, «What this does NOT settle»:

> **The real-data anchor is still untouched.** Oracle (d), Mathys 2019 (§8(d)), has not been run.

### 5.2 Что я установил

**Mathys 2019 отсутствует в CELLxGENE Discover целиком, а следовательно и в Census 2025-01-30.**

Метод проверки (воспроизводим, см. Приложение A.3): выгружен полный индекс `GET /curation/v1/datasets` — **2216 датасетов**; по объединённому полю `title + collection_name + collection_doi + collection_doi_label + citation` выполнен регистронезависимый поиск подстрок `mathys`, `rosmap`, `religious order`, `memory and aging`, `1195-2`, `s41586-019-1195`. **Совпадений — ноль.**

Контрольная выборка: все коллекции, где `disease` содержит Alzheimer или это слово есть в названии — **10 коллекций**, ни одна не является Mathys/ROSMAP:

| collection_id | Коллекция | Публикация |
|---|---|---|
| `1ca90a2d` | SEA-AD: Seattle Alzheimer's Disease Brain Cell Atlas | Gabitto 2024 Nat Neurosci |
| `c53573b2` | Cross-dementia human brain snRNA-seq | Rexach 2024 Cell |
| `0d35c0fd` | Molecular Signatures of Resilience to AD in Neocortical Layer 4 Neurons | Dharshini 2026 Nat Commun |
| `b953c942` | Single-soma transcriptomics of tangle-bearing neurons in AD | Otero-Garcia 2022 Neuron |
| `180bff9c` | Molecular characterization of selectively vulnerable neurons in AD | Leng 2021 Nat Neurosci |
| `84ce6837` | Population-scale cross-disorder atlas of the human prefrontal cortex | Lee 2024 medRxiv |
| `433700dc` | Brain vascular single-cell multi-omics | — |
| `7c4552fd` | Deciphering glial contributions to CSF1R-related disorder | Pan 2024 |
| `fcb3d1c1` | Live Human Microglia Single-cell RNA-seq | Olah 2020 Nat Commun |
| `8b35aa1f` | Gateway atlas of living human olfactory epithelium | Zhu 2026 bioRxiv |

Дополнительно: в полном манифесте (все 2190 строк, включая отсеянные) фигурирует **73 датасета**, все 73 присутствуют в текущем индексе Discover — то есть версия Census не «потеряла» датасет, которого сейчас нет.

**Атрибуция датасетов, названных в постановке задачи:**

| ID | Что это на самом деле | Mathys/ROSMAP? |
|---|---|---|
| `cff99df2` (27–31 vs 9–15) | Dharshini et al. (2026) Nat Commun, «Molecular Signatures of Resilience to AD in Neocortical Layer 4 Neurons»; PFC + precuneus + V1; snRNA + spatial | **Нет.** Провенанс доноров по Discover **не установлен** — описание коллекции источник тканей не называет. |
| `ac0c6561` | Rexach et al. (2024) Cell, cross-dementia (AD/PSP/Pick), BA4 + insula + V1 | **Нет** |
| `85c60876` (8v8) | Otero-Garcia et al. (2022) Neuron, «Single-soma transcriptomics of tangle-bearing neurons in AD — **Excitatory**» | **Нет** |
| `9813a1d4` (8v8) | Otero-Garcia et al. (2022) Neuron, то же — **Inhibitory** | **Нет** |
| `9f222629` | Sikkema et al. (2023) Nat Med, Human Lung Cell Atlas — **лёгкое, не мозг** | **Нет** |

### 5.3 Что это означает для проекта — и что я предлагаю

**Оракул §8(d) не может быть выполнен через Census-путь.** Спека это предусмотрела («else Synapse syn18485175»), но следствие архитектурное и его надо назвать вслух:

1. `io_counts.load_stratum(dataset_id, cell_type)` (§9 п.2) спроектирован под Census. Для Mathys понадобится **второй путь загрузки** — из Synapse syn18485175 (по публичным описаниям: 80 660 ядер, PFC/BA10, 48 доноров ROSMAP, 24 с AD-патологией vs 24 без; отфильтрованная матрица — syn18681734). Это работа, которой нет в плане §9 и которая должна быть оценена **до** заморозки списка, а не после.
2. Данные ROSMAP через Synapse требуют **отдельного соглашения о доступе (DUA)**. Это не техническая, а юридическая зависимость, и она лежит на критическом пути binding-проверки. Срок получения доступа — риск расписания, который стоит зафиксировать в этом же документе.
3. Целевой результат зафиксирован независимо: Murphy & Skene (2023), eLife 90214 — переанализ Mathys 2019 псевдобалком дал **в 549 раз меньше DEG при FDR 0.05**. Это и есть «качественное воспроизведение», которого требует §8(d).

**Особая роль якорного датасета — предложение.**
Mathys **не** входит в стратум-лист как строка списка и **не** должен в него входить: это оракул, а не измеряемый стратум. Предлагаю зафиксировать его отдельным разделом пререгистрации со статусом «BINDING, внешний по отношению к Census, зависит от DUA», и **не** ставить заморозку стратум-листа в зависимость от него — иначе весь список блокируется юридической процедурой.

**Ближайшие структурные суррогаты внутри Census** (на случай, если Synapse-путь окажется недоступен — как *дополнение*, никогда не *замена* binding-оракула):

| Датасет | Дизайн | Насколько близко к Mathys | Почему не замена |
|---|---|---|---|
| `85c60876` / `9813a1d4` Otero-Garcia 2022 | **ровно 8v8**, prefrontal cortex, Braak VI AD vs age-matched healthy | Ближайший: та же область, тот же тип контраста, донор-баланс 8v8 | Не snRNA-seq: FACS-отсортированные **отдельные сомы** (`suspension_type = cell`). Плюс `assay: v=0.577` и `depth: v=0.52` — сильные частичные конфаунды. Единственный стратум в `85c60876`: `neuron`, 96 129 клеток, cpd 6 562/5 624, cnts 4 591/**7 015**, rdf=14, perm=1.29e4 |
| `cff99df2` Dharshini 2026 | 31v15 (glutamatergic / GABAergic), PFC + precuneus + V1 | Донор-богатый AD, в тире σ=0.5 на двух стратумах | Три ассея в одном датасете (3′ v2 + 3′ v3 + Drop-seq), `assay: v=0.11–0.196`; провенанс доноров не установлен |
| `6c600df6` / `2727d83a` Leng 2021 | 7v3, SFG и энторинальная кора | Классический AD snRNA-seq | 7v3 → вне envelope при любом σ; `perfect_separation` по глубине на двух стратумах в `2727d83a` |

---

## 6. Envelope-разбор

### 6.1 Норма (Amendment 3 Change 1)

| σ_donor | min доноров/группу | Источник |
|---|---|---|
| 0.2 | 4 | Amendment 1 frontier (derived, \|err\| < 0.033); в сетке отсутствует |
| **0.35** | **8** | Amendment 1 frontier; сетка: `ebayes` power 0.793 при 8v8 → n\* ≤ 8 |
| **0.5** | **13** | Amendment 1 frontier; сетка: `ebayes` power 0.486 при 12v12 → n\* > 12 |
| 0.7 | 23 | Amendment 1 frontier; сетка: `ebayes` power 0.003 при 8v8 |

Формулировка амендмента, которую нельзя терять при чтении таблицы ниже:

> **Power 0.60 at `sigma_donor` = 0.5 with 8 donors per group remains UNMET and is NOT claimed.**
> …
> the honest one-line summary of the instrument's status after this amendment is *"valid within a stated envelope that most real strata may well fall outside"*, not *"valid"*.

### 6.2 Выживаемость предложенного списка по тирам

| # | Датасет | Стратумов-кандидатов | ≥4v4 (σ=0.2) | ≥8v8 (**σ=0.35**) | ≥13v13 (**σ=0.5**) | ≥23v23 (σ=0.7) | Потолок `min(A,B)` |
|---|---|---|---|---|---|---|---|
| 1 | `6f7fd0f1` SEA-AD DLPFC | 18 | 18 | **18** | **18** | **18** | 39 |
| 2 | `ac0c6561` Rexach | 27 | 24 | **23** | 0 | 0 | 11 |
| 3 | `d8da613f` Melms | 28 | 25 | **0** | 0 | 0 | 7 |
| 4 | `2a498ace` Yoshida | 47 | 46 | **34** | **27** | 0 | 20 |
| 5 | `ebc2e1ff` COMBAT | 25 | 23 | **21** | 0 | 0 | 10 |
| 6 | `f1606894` PERIHEART | 11 | 11 | **11** | **11** | **6** | 25 |
| 7 | `d18736c3` RA | 15 | 15 | **14** | **13** | 0 | 18 |
| 8 | `a12ccb9b` KPMP kidney | 37 | 34 | **29** | **25** | **6** | 24 |
| 9 | `19e46756` CHIP | 7 | 7 | **0** | 0 | 0 | 7 |
| 10 | `c893ddc3` opioid | 10 | 9 | **0** | 0 | 0 | 6 |
| 11 | `8e47ed12` Crohn gut | 18 | 15 | **0** | 0 | 0 | 7 |
| 12 | `4b6af54a` Emphysema | 8 | **0** | 0 | 0 | 0 | 3 |
| | **Итого датасетов** | | **11 / 12** | **7 / 12** | **5 / 12** | **3 / 12** | |

### 6.3 Что из этого следует для GO/NO-GO

Решающее правило §1 требует «majority of **independent datasets**» (D2, кластеризация по `dataset_id`). Значит:

- **При σ_donor ≈ 0.2:** 11 из 12 датасетов внутри envelope → большинство достижимо с запасом. Список работает.
- **При σ_donor ≈ 0.35:** внутри — 7 из 12 (SEA-AD, Rexach, Yoshida, COMBAT, PERIHEART, RA, KPMP). Большинство от 7 = 4 датасета. **Работоспособно, но запаса нет: выпадение любых двух по `integer_check` / `frozen_universe_size` ставит правило на грань.**
- **При σ_donor ≈ 0.5:** внутри — **5 из 12** (SEA-AD, Yoshida, PERIHEART, RA, KPMP). Из них RA — наиболее вероятный кандидат на вылет по размеру универсума (см. #7). **Фактически 4–5 датасетов на всё исследование.** Это ниже нижней границы §1 («8–12 datasets»), то есть **при σ_donor = 0.5 предложенный список формально не удовлетворяет §1**, и никакой другой список из этого манифеста удовлетворить не сможет — донор-богатых датасетов в кандидатах всего 5 сверх уже выбранных (см. §7 запасные).
- **При σ_donor ≈ 0.7:** внутри — 3 датасета (SEA-AD, PERIHEART, KPMP), из них два — nucleus/сердце-почка. Исследование в предзаявленной форме **невозможно**.

Это не аргумент за понижение планки. Это ровно тот исход, который Amendment 3 назвал заранее:

> **Whether the real sweep is feasible at all.** If real strata cluster at `sigma_donor` ≈ 0.5–0.7, the envelope admits them only at ≥ 13–23 donors per group. Whether enough CELLxGENE strata clear that is an open empirical question, **and a negative answer is a live outcome of this study, not a failure mode to be designed around.**

**Практическое предложение:** заморозить список в двух слоях —
**слой A (12 датасетов, выше)** как заявленный первый проход, и
**слой B — явную предзаявленную усечённую версию:** «если σ_donor-якорь даст ≥ 0.5, первый проход схлопывается до датасетов #1, #4, #6, #7, #8, и это отчитывается как ограничение области применимости, а не как выбор». Записав слой B заранее, мы лишаем себя возможности выбрать удобное подмножество постфактум.

### 6.4 Замечание о `permutation_count` как ложном ориентире

В манифесте `permutation_count` = C(D, n_A) от **суммарного** числа доноров. У COMBAT COVID-стратумов он равен 4.69e13 при контрольной группе в 10 доноров; у Melms — 8.88e5 при контроле в 7. Число велико, информации мало. Предлагаю в отчёте всегда печатать рядом `min(n_A, n_B)` и не использовать `permutation_count` как индикатор адекватности — иначе он систематически льстит косым дизайнам.

---

## 7. Литературная оценка «сильный vs тонкий» — и её источники

**Это суждение, а не измерение.** Ни один эффект здесь не измерен нами и не может быть измерен до прогона. Классификация нужна ровно для одного: обеспечить, что список покрывает обе стороны оси, как требует §1 (i)/(ii). Источник для каждого — собственная публикация датасета (DOI взят из Discover, не из памяти).

| # | Датасет | Оценка | Обоснование | Источник |
|---|---|---|---|---|
| 3 | Melms lethal COVID lung | **Сильный (максимум)** | Аутопсийное лёгкое при летальном COVID-19: диффузное альвеолярное повреждение, массивная миелоидная инфильтрация, фиброзная перестройка | Melms 2021 Nature, `10.1038/s41586-021-03569-1` |
| 1 | SEA-AD DLPFC | **Сильный–умеренный** | Нейродегенерация коры при AD-континууме: потеря уязвимых нейронных подтипов + глиальная реакция. Ослаблено тем, что термин Census — клинический `dementia` | Gabitto 2024 Nat Neurosci, `10.1038/s41593-024-01774-5` |
| 2 | Rexach cross-dementia | **Сильный** | Тауопатии (AD/PSP/Pick) с выраженным глиальным и нейронным ответом; дизайн специально построен на сравнении трёх диагнозов | Rexach 2024 Cell, `10.1016/j.cell.2024.08.019` |
| 5 | COMBAT COVID-плечо | **Сильный** | Острая тяжёлая вирусная инфекция в крови: интерферон-ответ, эмердженс-популяции моноцитов, плазмобласты | Ahern 2022 Cell, `10.1016/j.cell.2022.01.012` |
| 4 | Yoshida COVID-плечо | **Сильный** | То же — острая SARS-CoV-2 в PBMC | Yoshida 2022 Nature, `10.1038/s41586-021-04345-x` |
| 8 | KPMP AKI-плечо | **Умеренный–сильный** | Острое повреждение почки: выраженные «injured»-состояния проксимального канальца — центральный результат атласа | коллекция KPMP v1.5; DOI в Discover отсутствует — **не установлено** |
| 8 | KPMP CKD-плечо | **Умеренный** | Хроническая перестройка: фиброз, потеря подоцитов — медленнее и мягче, чем AKI | там же |
| 6 | PERIHEART AF | **Тонкий–умеренный** | Ремоделирование предсердия при AF; контроль — не здоровое сердце, а кардиохирургический пациент без AF, что дополнительно сужает контраст | Linna-Kuosmanen 2024 Cell Rep Med, `10.1016/j.xcrm.2024.101556` |
| 7 | Binvignat RA | **Тонкий** | Периферическая кровь при RA: сигнатуры на уровне субпопуляций и активности болезни, а не глобального сдвига (это и заявлено в названии работы) | Binvignat 2024 JCI Insight, `10.1172/jci.insight.178499` |
| 11 | Elmentaite Crohn (педиатр.) | **Тонкий–умеренный** | Ткань подвздошной кишки при детской болезни Крона; работа позиционирована как «транскрипционные связи», а не как описание массивного сдвига | Elmentaite 2020 Dev Cell, `10.1016/j.devcel.2020.11.010` |
| 10 | Phan opioid striatum | **Тонкий** | Опиоидная зависимость в дорсальном стриатуме: клеточно-специфичные программы, а не глобальная нейродегенерация | Phan 2024 Nat Commun, `10.1038/s41467-024-45165-7` |
| 4 | Yoshida post-COVID-плечо | **Тонкий** | Post-COVID-состояние в PBMC: клинически гетерогенный синдром без структурного повреждения ткани | Yoshida 2022 Nature, там же |
| 9 | Heimlich CHIP | **Тонкий (минимум)** | Клональный гемопоэз неопределённого потенциала — доклональное состояние без клинического фенотипа; заявленный эффект — активация воспалительных путей в отдельных типах клеток | Heimlich 2024 Blood Advances, `10.1182/bloodadvances.2023011445` |
| 12 | Wang emphysema | **Умеренный (нерелевантно)** | Включён ради дизайна 3v3, не ради эффекта | Wang 2023 Immunity, `10.1016/j.immuni.2023.01.032` |

**Методологический контекст (не датасетные источники, но опора всей рамки):**
- Murphy & Skene (2023), eLife 90214 — переанализ Mathys 2019: 549-кратное сокращение числа DEG при переходе на псевдобалк. Это целевой качественный результат §8(d).

**Как этим НЕ пользоваться.** Эта таблица не является предсказанием λ_naive и не должна сравниваться с результатами постфактум как «проверка гипотезы». Единственное её назначение — доказать, что список не подобран под win.

---

## 8. Чего это предложение НЕ решает

Раздел написан в жанре «What this does NOT settle» из `AMENDMENTS.md`, намеренно.

1. **Допуск в sweep закрыт, и этот документ его не открывает.** `admitted_to_sweep = False` у всех 2190 строк; блокеры — `integer_check`, `frozen_universe_size`, `sigma_donor_estimate`, `envelope_membership`. Первые два вычисляются на загрузке X и в этом артефакте отсутствуют; последние два упираются в **σ_donor-якорь, который Amendment 3 оставил OPEN третий амендмент подряд**. Никакой стратум отсюда не допущен, и никакая часть этого документа не является допуском.

2. **Членство в envelope не установлено ни для одного стратума.** Колонки «σ=0.35 / 0.5» в §6 — это **арифметика по числу доноров при гипотетическом σ**, а не измерение. Пока `sqrt(s0^2) · ln 2 → donor_sigma` не выведено и не провалидировано против симулятора (Amendment 3 явно называет это «required work that this amendment does not do»), реальный σ каждого стратума неизвестен, и таблица §6.2 читается как сценарный анализ, не как результат.

3. **Пререгистрация — отдельный человеческий акт.** §1: «Pre-register the stratum list before computing any metric». Этот файл — предложение к заморозке. Заморозка = решение Александра + фиксация списка в репозитории с датой и хешем **до** любого вычисления метрик по этим стратумам. Пока этого не произошло, список можно менять; после — нельзя, и изменения идут через `AMENDMENTS.md`.

4. **Критерий «сильный/тонкий» — литературное суждение.** См. §7. Он может оказаться неверным для конкретного датасета, и это не будет сбоем: назначение оси — покрытие, а не предсказание.

5. **`pooled` не разрешён ни для одного стратума.** Флаг `pooled: unresolved — no pool/library id in obs (D3 lower bound)` стоит на всех 1197 кандидатах. По §1/D3 это значит: донор-псевдобалк на **любом** из предложенных стратумов остаётся *нижней оценкой* правильной единицы репликации, и заявление «донор-псевдобалк калиброван» (золотой стандарт) ни на одном из них сделать нельзя. Выбор датасетов этого не чинит — это свойство того, какие колонки Census отдаёт в этом пине.

6. **Bins для cells-per-donor в спеке не определены.** §1 требует список, покрывающий «the pre-registered bins (D1)», §7 п.3 требует «pre-register the cells-per-donor bins», решающее правило (строка 13) ссылается на «the floor of the pre-registered bins» — но **нигде в `PHASE0_SPEC.md` или `AMENDMENTS.md` численные границы bins не заданы** (проверено grep'ом, см. Приложение A.5). Следовательно **строго говоря, требование §1 «spanning the pre-registered bins» сейчас невыполнимо**: не существует объекта, который можно перекрыть. Я показал перекрытие диапазона (18 … 6 071 клеток на донора), но это не то же самое. **Bins должны быть пререгистрированы тем же актом, что и стратум-лист, иначе оба пререгистрации неполны.**

7. **`cell_type_ontology_depth` = `pending` у всех строк.** D5 требует не пулить headline через стратумы разной гранулярности. Конфликт виден невооружённым глазом: SEA-AD даёт `L2/3-6 IT`, `chandelier pvalb`, `sncg GABAergic`; Rexach даёт `astrocyte`, `glutamatergic neuron`. **Оба датасета в списке, и пулить их headline нельзя.** Гармонизирующий уровень должен быть выбран до заморозки.

8. **Регион ткани внутри датасета не проверен.** У Rexach (BA4/insula/V1), KPMP (cortex/medulla/papilla), Dharshini (PFC/precuneus/V1) регион — это батч, который `tissue_general` не различает. Скрин его не поймал не потому, что его нет, а потому, что он не в тех колонках. До заморозки надо решить, входит ли регион в определение стратума или в ковариаты.

9. **Баланс по `sex`, `development_stage`, `self_reported_ethnicity` не проверен.** Эти колонки запрошены в obs (`obs_columns_requested`), но конфаунд-скрин §1 их не покрывает — он проверяет только {assay, suspension_type, tissue_general, depth bin, pool}. Для Yoshida (дети + взрослые) и Elmentaite (4–14 лет) возраст — очевидный кандидат в конфаундеры. **Это дыра в §1, а не в манифесте**, и её стоит закрыть отдельным пунктом до заморозки.

10. **Ни один из 12 датасетов не проверен на целочисленность counts и на размер универсума.** Оба — блокеры допуска, оба вычисляются только на загрузке X. Наиболее вероятный кандидат на вылет по универсуму — #7 (RA, 300–800 counts/cell). Список надо считать предложением *до* этой проверки, и после неё он может сократиться.

11. **§8(d) остаётся невыполненным и внешним.** См. §5.3: Census-путь закрыт, Synapse-путь требует DUA. Это binding-проверка, и она не должна блокировать заморозку списка, но её статус обязан быть записан в тот же документ.

---

## 9. Запасные датасеты (5)

По одной строке — почему второй.

| dataset_id | Что это | Env | Почему второй |
|---|---|---|---|
| `1c739a3e` | Kuppe 2022 Nature, инфаркт миокарда, snRNA, heart left ventricle, 16v4 | 0.2 | Сильный эффект и хорошая ткань, но контроль зафиксирован на 4 донорах → выпадает уже при σ=0.35, а сердце уже закрыто #6 при 25v29. |
| `cff99df2` | Dharshini 2026 Nat Commun, AD, PFC/precuneus/V1, 31v15 | 0.5 (2 стратума) | Донор-богатый AD и почти второй SEA-AD, но три ассея в одном датасете (3′ v2 + 3′ v3 + Drop-seq) и сильный перекос групп; берётся, если SEA-AD отвалится на `integer_check`. |
| `3c75a463` | Rodríguez-Ubreva 2022 Nat Commun, CVID, PBMC, **10x 5′ v1**, 10v8 | 0.35 (20 стратумов) | Чистейшее 5′ 8v8+ с ровными группами, но это **активированные** PBMC (in-vitro стимуляция) — лишний слой между болезнью и транскриптомом; берётся, если понадобится третье 5′-плечо. |
| `c7775e88` | «Single-cell multi-omics analysis of the immune response in COVID-19»; Discover атрибутирует коллекции *CITIID-NIHR COVID-19 BioResource Collaboration et al. (2021) Nat Med*, `10.1038/s41591-021-01329-2`; кровь, 10x 3′ transcription profiling, до 86v29 | 0.5 (21 стратум) | Донор-богатый и в тире σ=0.5, но полностью дублирует ось «COVID в PBMC», уже занятую #4 и #5; берётся как 3′-замена одному из 5′, если понадобится сместить баланс ассеев. |
| `9f222629` | Sikkema 2023 Nat Med, Human Lung Cell Atlas (full), 191 стратум, 8 ассеев, 76 стратумов ≥13v13 | 0.5 (76 стратумов) | Крупнейший пул кандидатов и единственный источник Drop-seq/Seq-Well, но `assay` V ≈ 0.34–0.86 и `suspension_type` V ≈ 0.27–0.67 почти на всех донор-богатых стратумах: мета-атлас, где ассей частично разделяет условия. Берётся **только** как сознательный «трудный конфаундированный дизайн» с явной ковариатой, никогда как рядовая строка. |

---

## Приложение A. Как проверить каждое утверждение

Все команды воспроизводимы с нуля. `$S` = каталог со скрэтчпадом.

### A.1 Пересчёт сводки манифеста (§1)

```bash
python - <<'PY'
import json, collections
d = json.load(open(r"$S\census_candidates_full.json", encoding="utf-8"))
rows = d["rows"]
print("rows", len(rows))
print("gate_status", collections.Counter(r["gate_status"] for r in rows))
print("admitted", sum(r["admitted_to_sweep"] for r in rows))
print("blockers", collections.Counter(b for r in rows for b in r["admission_blockers"]))
print("datasets", len({r["dataset_id"] for r in rows if r["gate_status"]=="candidate"}))
PY
```
Ожидается: `rows 2190`; `candidate 1197`, `excluded_inclusion_gate 981`, `excluded_confound 12`; `admitted 0`; каждый из четырёх блокеров = 2190; `datasets 68`.

### A.2 Тиры envelope (§1.4, §6.2)

```bash
python - <<'PY'
import json, collections
d = json.load(open(r"$S\census_candidates_full.json", encoding="utf-8"))
cand = [r for r in d["rows"] if r["gate_status"] == "candidate"]
t = collections.Counter()
for r in cand:
    n = min(r["n_donors_A"], r["n_donors_B"])
    t["3" if n < 4 else ">=4" if n < 8 else ">=8" if n < 13 else ">=13"] += 1
print(t)                      # ожидается 3:180, >=4:463, >=8:243, >=13:311
for ds in ["6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3"]:      # подставить любой id
    rs = [r for r in cand if r["dataset_id"] == ds]
    print(ds, len(rs),
          sum(min(r["n_donors_A"], r["n_donors_B"]) >= 8  for r in rs),
          sum(min(r["n_donors_A"], r["n_donors_B"]) >= 13 for r in rs),
          max(min(r["n_donors_A"], r["n_donors_B"]) for r in rs))
PY
```

### A.3 Отсутствие Mathys/ROSMAP в Discover (§5.2) — ключевая проверка

```bash
curl -s --compressed "https://api.cellxgene.cziscience.com/curation/v1/datasets" -o "$S/cxg_datasets_index.json"

python - <<'PY'
import json
idx = json.load(open(r"$S\cxg_datasets_index.json", encoding="utf-8"))
print("datasets in Discover:", len(idx))          # на 2026-08-16: 2216
pats = ["mathys", "rosmap", "religious order", "memory and aging", "1195-2", "s41586-019-1195"]
hits = [e for e in idx
        if any(p in " ".join(str(e.get(k)) for k in
              ("title","collection_name","collection_doi","collection_doi_label","citation")).lower()
              for p in pats)]
print("hits:", len(hits))                          # ожидается 0
for e in hits: print(e["dataset_id"], e["title"])
PY
```
Утверждение §5.2 верно ⟺ `hits: 0`. Если однажды станет ≠ 0 — Mathys появился в Discover и §8(d) можно вернуть на Census-путь.

### A.4 Метаданные любого предложенного датасета (assay / suspension / tissue / публикация)

```bash
python - <<'PY'
import json
idx = {x["dataset_id"]: x for x in json.load(open(r"$S\cxg_datasets_index.json", encoding="utf-8"))}
e = idx["a12ccb9b-4fbe-457d-8590-ac78053259ef"]        # подставить любой id из §3
lab = lambda v: sorted({x["label"] if isinstance(x, dict) else str(x) for x in v})
for k in ("title", "collection_name", "collection_doi", "collection_doi_label", "cell_count"):
    print(k, "=", e.get(k))
print("assay      =", lab(e["assay"]))
print("suspension =", lab(e["suspension_type"]))
print("tissue     =", lab(e["tissue"]))
PY
```
Альтернатива через коллекцию (даёт описание/провенанс):
`curl -s "https://api.cellxgene.cziscience.com/curation/v1/collections/{collection_id}"`
Замечание: `GET /curation/v1/datasets/{dataset_id}` **возвращает 404** — датасет по одному id этим маршрутом не отдаётся; работает только индекс или маршрут через коллекцию.

### A.5 Отсутствие численных bins для cells-per-donor (§8 п.6)

```bash
cd <REPO>
grep -n -iE "bin(s)?\b" docs/PHASE0_SPEC.md
grep -n -iE "cells.per.donor bin|pre-register.*bin" docs/AMENDMENTS.md
```
Ожидается: в спеке 4 попадания (строки 13, 71, 75, 192), из них строка 71 — про `sequencing-depth bin` (другое), остальные три ссылаются на «pre-registered bins», **не задавая их численно**; в `AMENDMENTS.md` — ноль попаданий.

### A.6 Все стратумы конкретного датасета с полями манифеста (карточки §4)

```bash
python - <<'PY'
import json
ds = "ac0c6561-7a48-4185-af6f-af799f699172"      # подставить любой id
d = json.load(open(r"$S\census_candidates_full.json", encoding="utf-8"))
for r in d["rows"]:
    if r["dataset_id"] != ds or r["gate_status"] != "candidate":
        continue
    A, B = r["cells_per_donor_by_group"]["A"], r["cells_per_donor_by_group"]["B"]
    print(f'{r["disease"][:26]:28s} {r["cell_type"][:36]:38s} '
          f'A={r["n_donors_A"]:3d} B={r["n_donors_B"]:3d} n={r["n_cells"]:7d} '
          f'cpd={A["median"]:.0f}/{B["median"]:.0f} '
          f'cnts={r["median_counts_per_cell_by_group"]["A"]:.0f}/'
          f'{r["median_counts_per_cell_by_group"]["B"]:.0f} '
          f'rdf={r["residual_df"]} perm={r["permutation_count"]:.2e} '
          f'| {[f for f in r["confound_flags"] if not f.startswith("pooled")]}')
PY
```
Каждая строка любой таблицы в §4 воспроизводится этой командой дословно.

### A.7 Цитаты §1 и Amendment 3

```bash
sed -n '56,78p'   <REPO>/docs/PHASE0_SPEC.md   # §1 целиком
sed -n '201,214p' <REPO>/docs/PHASE0_SPEC.md   # §8 оракулы, (d) на 212–213
sed -n '598,711p' <REPO>/docs/AMENDMENTS.md    # Amendment 3 Change 1 + σ_donor OPEN
sed -n '745,771p' <REPO>/docs/AMENDMENTS.md    # Amendment 3 "What this does NOT settle"
```

### A.8 Флаг `pooled` на 100% кандидатов (§1.5, §8 п.5)

```bash
python - <<'PY'
import json
d = json.load(open(r"$S\census_candidates_full.json", encoding="utf-8"))
cand = [r for r in d["rows"] if r["gate_status"] == "candidate"]
print(sum(any(f.startswith("pooled") for f in r["confound_flags"]) for r in cand), "/", len(cand))
print({r["pooled_flag"] for r in cand})
PY
```
Ожидается `1197 / 1197` и `{'unresolved'}`.

---

*Составлено по манифесту от 2026-08-15T22:18:37Z и индексу CELLxGENE Discover от 2026-08-16. Ни одно число в этом документе не выведено из памяти: всё либо пересчитано из манифеста, либо взято из ответа Discover API, либо процитировано из `docs/`. Там, где источник не ответил — сказано «не установлено». Литературные оценки силы эффекта помечены как суждения и снабжены DOI.*
