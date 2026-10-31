# Piano d'attacco per il paper — MathBehaviours

*Nota di lavoro, 31/10/2026. Non è il paper: è la ricognizione della situazione e il piano per scriverlo.*

> **Placeholder verificati (31/10/2026):** completata la figura del disegno con sequenze T0/T1 e schermate originali delle slide (`figures/design.py` la rigenera). Per i self-report è stata aggiunta una fonte pertinente: Quintero et al. (2022), doi:10.1111/bjep.12482; tolta l'affermazione non documentata che fossero la maggioranza delle misure. Verificato il testo completo di Wang et al. (2025), doi:10.1016/j.lindif.2025.102681: la discussione considera campione adulto, AMAS e scelte dopo errore; anche quel compito forniva hint al computer. Nei dati deduplicati le 363 scelte mancanti hanno tutti e tre gli indicatori NA e accuratezza già pari a 0; la causa resta da confermare con i log. I 547 ID sono presenti in entrambe le sessioni: l'export non documenta il flusso dei bambini non inclusi. Restano aperti età, reclutamento/esclusioni, superamento della stopping rule, regione, protocollo etico, link OSF/deposito, funding, conflitti e ruoli, oltre alle verifiche bibliografiche specifiche degli strumenti.

> **Aggiornamento dopo le analisi finali (31/10/2026):** `paper.qmd` usa ora esclusivamente le sintesi finali in `analysis/results/` (N = 547), con tabelle generate dal codice e figura delle probabilità predette. Le ipotesi di lavoro e le questioni analitiche sotto restano come ricognizione storica; lo stato attuale è il seguente.
>
> - Screening ordinale con indici scalati e nessuna covarianza residua: avoidance escluso (CFI = .923); buoyancy ammesso. M1 aggiunge quindi solo buoyancy (ΔELPD = 0.28, SE = 0.70), senza superare M0. L'effetto specifico per categoria di worry è mantenuto (21.81, SE = 7.24); quello di competenza percepita no (−0.24, SE = 1.91).
> - Il multinomiale resta preferito al cumulativo completo (ΔELPD del cumulativo = −718.83, SE = 35.87) e all'adjacent-category con effetti specifici (−847.31, SE = 39.21). Sono riportati tutti e tre i contrasti. Worry, self vs help: −0.36 [−0.57, −0.14].
> - L'esclusione di avoidance dipende dallo screening, non da LOO; la sua inclusione nei follow-up è dichiarata esplorativa. Nei binari completi l'associazione con self vs non-self è incerta (p = .076), mentre quelle con skip vs help e skip vs non-skip rimangono presenti. Buoyancy mostra il pattern di soppressione già ipotizzato.
> - LOO lascia fuori trial di bambini già osservati: non valuta nuovi bambini e non dimostra l'assenza di associazioni tra bambini. Il vantaggio del multinomiale riguarda le specificazioni confrontate. Il vantaggio di accuratezza dopo help è osservazionale.
> - Gli script non salvano più fit o catene; producono automaticamente solo le sintesi necessarie. Workspace R disattivato nel progetto. Le catene precedenti sono state eliminate dopo la sintesi.
>
> Restano da recuperare: età (colonna vuota), flusso di reclutamento/esclusioni e stopping rule, metadati istituzionali/OSF, log di timeout/scelte mancanti/prompt, item delle altre scale e verifiche bibliografiche già indicate. Il dataset fornito contiene solo gli item di buoyancy e avoidance; l'affidabilità delle altre scale non è stimabile da questo export.

> **Aggiornamento 31/10/2026 (pomeriggio):** il repository è stato riorganizzato (vedi `README.md`). I percorsi citati sotto sono quelli *precedenti*: `Analisys_truedata/` → dati in `data/`, resto in `archive/first_real_data_run/`; `Script/`, `Data/`, `outputs/`, `index.*`, PDF della prereg → `Preregistration/`. Gli script finali sono in `analysis/` (§9 punti 2–5). `vo4a16` è risultato un export duplicato (righe identiche): N = 547.
*Fonti consultate: `Slide/7_Pagan_AIRIPA2026.pdf`, `OSF final preregistration text.pdf`, `index.Qmd` (note di progetto, nov. 2025), `Analisys_truedata/` (script, dati, oggetti `.rds`), `Script/`, `Data/`, `outputs/`, storia git, più alcune descrittive ricalcolate al volo sui dati reali.*

---

## 1. Stato attuale: cosa c'è e dove

| Cartella / file | Cosa contiene davvero | Note |
|---|---|---|
| `Data/SimulatedDatasetBehaviours.csv` | Dataset **simulato** della preregistrazione | — |
| `Script/` | Materiali della **preregistrazione**: simulazione dati, power analysis (ordinale CLMM e categoriale "fast"), script pre-pianificati (`Planned_analysis.R`, `New_PrePlannedAnalysis.R`, `Completed_PrePlannedAnalysis_from_modelselection.R` = la pipeline preregistrata) | `Completed_…` legge ancora il dataset simulato |
| `outputs/` | LOO/WAIC categoriale vs ordinale su **dati simulati** (9000 righe; Δelpd = −728) | Non sono i risultati reali: possono generare confusione |
| `Analisys_truedata/DatasetBehaviours.csv` | **Dati reali** (formato long, T0 + T1) | |
| `Analisys_truedata/Analysis_DEF.R` | Analisi che *in parte* sta dietro alle slide (LOO categoriale vs ordinale; CFA; modello di mediazione categoriale) | Vedi problemi §4 |
| `Analisys_truedata/*.rds`, `results_*.csv` | Output della **pipeline preregistrata (ordinale) eseguita sui dati reali** | Vedi §2.2 — risultati mai entrati nelle slide |
| `index.Qmd/html` | Note di progetto nov. 2025: RQ1 = ordinale vs categoriale, RQ2 = predittori | È questa la struttura che le slide seguono, non la prereg finale |
| `Template_OSF_Behaviours_020426.pdf` | Bozza precedente della prereg | |
| `CheckDatiFeedback.R` | Script di un altro progetto (gruppo "CG", feedback) | Fuori tema, probabilmente da rimuovere |
| `results5000_*_v2.RData` in root | Residui della power analysis | Da spostare in `Script/` o eliminare |

Il codice che produce **la tabella e i grafici delle slide** (coefficienti multinomiali, probabilità predette, regressione lineare di robustezza, descrittive per genere/classe, matrice di correlazione) **non è nel repo** (in parte sembra fatto in Python: la heatmap ha lo stile matplotlib). È la prima cosa da recuperare da Riccardo.

---

## 2. Le due "storie" a confronto

### 2.1 Preregistrazione (OSF, finale)

- **RQ1 (prereg)**: atteggiamenti di tratto e di stato (worry, competenza percepita di stato, buoyancy, evitamento) predicono la scelta (skip / help / self), controllando abilità T0 e genere? Gli atteggiamenti di stato mediano i corrispondenti di tratto (ansia matematica → worry; self-competence → competenza percepita di stato)?
- **RQ2 (prereg)**: all'interno della rappresentazione **ordinale** (skip < help < self), alcuni predittori richiedono effetti category-specific (violazione degli odds proporzionali)?
- Assunzione di lavoro: la scelta è **ordinale**. Il modello categoriale è previsto **solo** come follow-up "se la struttura ordinale si rivela chiaramente inadeguata", da riportare come analisi **supplementare non preregistrata** (p. 3 della prereg). C'è però una power analysis anche per il modello categoriale.
- Pipeline: CFA (buoyancy, evitamento; CFI > .95, RMSEA < .08; alpha ordinale) → M0 (solo stato + covariate, = mediazione totale) → M1 (+ buoyancy, evitamento) → M2 (+ diretti di tratto: mediazione parziale) con criterio |Δelpd| > 2 SE → estensioni `cs()` → test direzionali one-tailed con FDR, in un quadro "frequentista pragmatico" (MCMC solo come strumento di stima) → effetti indiretti a×b.
- Predittori standardizzati (z).

### 2.2 Cosa ha dato la pipeline preregistrata sui dati reali (dagli `.rds` in `Analisys_truedata/`)

- Righe analizzate: 17 171 (= 17 536 − 365 scelte mancanti).
- **M0 vs M1**: M1 (con buoyancy ed evitamento) è *peggiore* di M0, Δelpd = −3.4 (SE 1.1) → si tiene **M0**; M2 non stimato; mediazione "totale" per default.
- **Stadio 2**: `cs(worry)` preferito agli odds proporzionali, Δelpd = 21.4 (SE 7.3), ≈ 2.9 SE. È già un segnale di violazione dell'ordinalità *dentro* la pipeline preregistrata.
- Effetti: worry ≈ 0 sulla soglia skip|help, **negativo sulla soglia help|self** (−0.24 [−0.38, −0.09]); competenza percepita di stato ≈ 0. Indiretto ansia → worry → scelta credibile solo sulla seconda soglia.
- Conseguenza importante: **seguendo la prereg alla lettera, l'evitamento — il "predittore più consistente" delle slide — sarebbe stato escluso dal modello finale.** Questo va gestito apertamente (§3).

### 2.3 Slide AIRIPA 2026 (la storia da raccontare)

- Ri-etichettatura delle RQ: **RQ1 (slide) = natura della misura** (ordinale vs categoriale), **RQ2 (slide) = predittori**. È l'inverso della prereg e ricalca `index.Qmd`.
- N = 547 (297 F); 4ª classe 378, 5ª 169. Distribuzione scelte: self 59.2%, help 30.3%, skip 8.3%, NA 2.1%.
- **Categoriale ≫ ordinale**: Δelpd = −721.1, SE 35.9 (≈ 20 SE).
- Modello multinomiale (riferimento = skip), tutti i predittori di stato e tratto + abilità + genere, intercetta random soggetto:
  - *self vs skip*: abilità T0 + ; **evitamento −**; **buoyancy −** (direzione opposta all'ipotesi); worry, competenza percepita (stato e tratto), ansia di tratto ≈ 0.
  - *help vs skip*: **worry (stato) +** (0.13 [0.05, 0.21]); evitamento − (marginale, PP = .95); il resto ≈ 0.
- Robustezza: regressione lineare sul punteggio per soggetto (skip = 0, help = 1, self = 2, z): abilità .51, worry −.13, buoyancy −.11, evitamento −.10; R² = .33.
- Conclusioni: le scelte si comportano come categorie, non come un continuum di engagement; l'evitamento di tratto è il predittore più consistente; il worry aumenta specificamente la richiesta di aiuto; risultato inatteso sulla buoyancy; "aiuto come evitamento mascherato".

---

## 3. Come conciliare prereg e slide: la cornice proposta

La deviazione è **legittima e in parte anticipata dalla prereg stessa** (che prevede il follow-up categoriale e ha una power analysis categoriale), ma va dichiarata con molta trasparenza. Proposta:

1. **Il risultato centrale del paper diventa un risultato sulla misura**: la sequenza skip < help < self, assunta come continuum di engagement nella prereg, *non regge*. Evidenze convergenti: (a) LOO enorme a favore del categoriale; (b) anche dentro la pipeline preregistrata gli odds proporzionali cadono (`cs(worry)`); (c) a livello descrittivo help e skip hanno profili di correlazione quasi identici con i predittori (vedi §5), cioè l'aiuto non sta "in mezzo".
2. **Riportare la pipeline preregistrata per intero** (non in appendice nascosta): una sezione "Preregistered analyses" con CFA, M0/M1/M2, `cs()`, test FDR e indiretti, onestamente, inclusa l'esclusione di buoyancy/evitamento in M1. Poi "Deviation from preregistration and follow-up analyses", che giustifica il passaggio al multinomiale con l'evidenza del punto 1.
3. **Il modello multinomiale è presentato come analisi principale *esplorativa / non preregistrata*** (oppure: "prespecified as contingency, but not part of the preregistered decision pipeline"). Le inferenze vanno formulate di conseguenza (stime + intervalli, niente linguaggio confermativo).
4. **Tabella di deviazione dalla prereg** (formato raccomandato, es. Willroth & Atherton, 2024): ogni deviazione, motivo, impatto.
5. Le RQ del paper: usare la numerazione delle slide (RQ1 = struttura della misura; RQ2 = predittori), ma con una nota esplicita che nella prereg l'ordine era inverso e la domanda sulla struttura era formulata come domanda sugli odds proporzionali. In alternativa, tenere la numerazione della prereg per minimizzare attriti con i revisori. **Da decidere con Riccardo/Sara** — io propenderei per la numerazione della prereg + una RQ "strutturale" dichiarata come derivata dalla RQ2 preregistrata.
6. **Mediazione**: le slide la abbandonano. Nel paper va comunque riportata (era preregistrata), ma ridimensionata: il disegno non la sostiene bene (vedi §4, ordine temporale). Proposta: risultati della mediazione preregistrata in breve + discussione dei limiti; nessuna enfasi.

---

## 4. Problemi analitici da risolvere prima di scrivere

In ordine di priorità.

### Alta priorità (possono cambiare i risultati)

1. **Predittori non standardizzati in `Analysis_DEF.R`.** La prereg prevede z-scores; i grafici delle slide dicono "(z)", ma lo script fitta sulle scale grezze (le β in tabella sembrano per punto grezzo: es. abilità 0.14 per item corretto ≈ 1.1 per DS). Rifittare con predittori standardizzati, e chiarire nelle tabelle.
2. **Effetto "inatteso" della buoyancy = probabile soppressione.** A livello bivariato la buoyancy correla *positivamente* con la proporzione di self (r = .12) e negativamente con skip (−.09) e help (−.10). Il segno negativo appare solo condizionando su worry, competenza di stato, ansia, abilità (con cui la buoyancy correla .30–.46). Prima di interpretarlo come "più buoyancy → meno autonomia" servono: modelli un predittore alla volta (+ covariate), VIF, e una discussione esplicita del fenomeno di soppressione. Probabilmente nel paper diventa "nessun effetto unico robusto" più che "effetto opposto".
3. **Struttura di clustering incompleta.** Solo `(1|id)`. Mancano:
   - **classe** (37 classi, ricavabili dal prefisso dell'id, es. `dv4a`; la colonna `Class` è tutta NA) → `(1|class)` almeno come sensitivity;
   - **item**: gli item sono campionati da un pool di 100, `item` è nel dataset → `(1|item)` incrociato. La difficoltà dell'item è probabilmente un determinante forte della scelta.
4. **Un ID duplicato**: `vo4a16` ha 64 righe nella sessione 2 (32 in T0). Probabilmente due bambini con lo stesso codice o un doppio export. Da chiarire con chi ha raccolto i dati e correggere a monte.
5. **La robustezza con regressione lineare usa il punteggio skip = 0 / help = 1 / self = 2**, cioè proprio l'ordinamento che il paper sostiene essere sbagliato. Va sostituita (o affiancata) da analisi coerenti con la conclusione categoriale, ad es.:
   - modelli binari separati (self vs non-self; skip vs non-skip; help vs non-help), o proporzioni per soggetto con modelli beta/binomiali;
   - modello multinomiale frequentista (es. `mclogit`/`nnet` + cluster SE) come check di stima.
6. **Riportare tutti e tre i contrasti** (anche *self vs help*, che manca nelle slide ed è il più informativo sull'"aiuto come evitamento mascherato") e, soprattutto, **effetti sulla scala delle probabilità** (average marginal effects per 1 DS su P(skip), P(help), P(self)), che rendono il multinomiale leggibile e confrontabile con le ipotesi.

### Media priorità (solidità / interpretazione)

7. **Accuratezza molto bassa**: tra le risposte self l'accuratezza è 26.7%, tra le help 27.5% (abilità T0 media 11.8/32). Va capito: limite dei 25 s? codifica delle risposte in ritardo? Ha implicazioni per la validità del compito (scegliere "self" non significa riuscirci) e per l'interpretazione della scelta come "engagement". Analisi secondaria utile: l'hint migliora l'accuratezza a parità di abilità?
8. **Scelte mancanti (365, 2.1%)**: cosa sono? timeout prima della scelta? Descriverle; eventualmente sensitivity che le tratta come quarta categoria o come skip.
9. **Casi estremi**: 57 bambini scelgono sempre self, 3 sempre help, nessuno sempre skip; 215 non saltano mai. Descriverli (la prereg dice di non escluderli) e verificare la stabilità della varianza random.
10. **Ordine temporale delle misure** (rilevante per mediazione e per l'interpretazione causale):
    - T0: 16 item → worry e competenza di stato → 16 item → ansia di tratto e buoyancy (quindi lo stato è misurato *prima* del tratto);
    - T1: compito di scelta → **evitamento e self-competence di tratto misurati *dopo* il compito**. Possibile influenza dell'esperienza del compito sulle risposte; va dichiarato tra i limiti.
11. **Dinamiche nel compito**: effetti di posizione/blocco (fatica, pausa), e i prompt automatici ogni 5 skip/help consecutivi (che per costruzione alterano le scelte successive). Almeno descrittivi; eventualmente `block` come covariata.
12. **Confronto ordinale vs categoriale: sondare alternative intermedie**, per non contrapporre solo due estremi. Proposta minimale: (a) cumulativo con `cs()` su tutti i predittori; (b) modello sequenziale (prima "provo da solo vs no", poi "aiuto vs salto"), che traduce direttamente l'idea di processo decisionale a due stadi; (c) ordinale con ordine alternativo. Se il categoriale resta superiore, la conclusione è molto più forte; se il sequenziale va bene quanto il categoriale, offre un'interpretazione psicologica più ricca.
13. **Quadro inferenziale da uniformare.** Prereg: "frequentista pragmatico", p one-tailed + FDR. Slide: posterior probability bayesiane. Scegliere una lingua sola: proposta = stime posteriori con CrI al 95% come principale; per gli effetti preregistrati riportare anche p one-tailed e FDR come da prereg.
14. **Validità psicometrica**: riportare CFA (buoyancy, evitamento — come da prereg) e affidabilità per tutte le scale (ordinal alpha/omega), oltre alle correlazioni.
15. **Genere ed età**: `Age` è tutto NA → recuperare l'età o usare la classe (4ª/5ª) come covariata/descrittiva. Differenze di genere marcate su ansia e competenza percepita (descrittive nelle slide) → menzionare; eventuali moderazioni solo come esplorative (la prereg le prevede come tali).

### Bassa priorità / bug di codice in `Analysis_DEF.R`

- `loo(..., resp = "choicec")` funziona solo nel ramo categoriale (nel ramo ordinale la risposta è `choiceord`).
- `write.excel()` non esiste in R base; `h_diff` esiste solo nel ramo categoriale.
- Il modello di confronto RQ1 (`fit_cat`/`fit_ord`) include tutti i predittori; il modello di mediazione ne include un sottoinsieme: chiarire quale modello genera la tabella delle slide.
- I file `results_*.csv` e gli `.rds` provengono da un run diverso (pipeline ordinale) rispetto alle slide: rinominarli per evitare confusione.

---

## 5. Evidenze già emerse dalle verifiche rapide (da confermare nella pipeline definitiva)

Correlazioni bivariate per soggetto (N = 547) tra proporzione di scelte e predittori grezzi:

| | buoyancy | evitamento | ansia (tratto) | SC tratto | worry | SC stato | abilità T0 |
|---|---|---|---|---|---|---|---|
| P(skip) | −.09 | .21 | .16 | −.15 | .12 | −.21 | −.40 |
| P(help) | −.10 | .14 | .19 | −.13 | .26 | −.22 | −.41 |
| P(self) | .12 | −.20 | −.22 | .17 | −.27 | .27 | .50 |

Lettura: **help e skip si comportano in modo quasi identico** rispetto a tutti i predittori, e il contrasto vero è *self vs non-self*. È un argomento descrittivo semplice e forte a favore della non-ordinalità (e dell'idea di "help come evitamento mascherato"), ed è coerente con un modello sequenziale "provo da solo? → se no, aiuto o salto". Il worry è l'unica variabile che sembra distinguere un po' help da skip. Buoyancy bivariata nella direzione attesa (vedi §4.2).

---

## 6. Messaggio del paper (bozza di tesi)

> Un paradigma comportamentale di scelta item per item (skip / help / self) mostra che, nella scuola primaria, chiedere aiuto e saltare non sono gradi diversi di uno stesso continuum di engagement: le scelte si comportano come categorie qualitativamente distinte, e la distinzione principale è tra affrontare il problema da soli e non farlo. L'abilità matematica è di gran lunga il determinante più forte; tra gli atteggiamenti, l'evitamento di tratto riduce la scelta autonoma e il worry di stato è associato in particolare alla richiesta di aiuto. Gli effetti attesi di competenza percepita e buoyancy non emergono come effetti unici.

Contributi da rivendicare: (1) metodo — una misura comportamentale oggettiva in età primaria e un test empirico della sua struttura, con una lezione generale (non assumere l'ordinalità di categorie comportamentali); (2) sostanza — quali atteggiamenti contano per quali scelte; (3) trasparenza — prereg eseguita e riportata, deviazione giustificata dai dati.

Formulazioni da evitare (presenti nelle slide): "evitamento e buoyancy non sono estremi di uno stesso continuum" (confonde la struttura della **misura comportamentale** con i costrutti di tratto); "buoyancy riduce l'autonomia" (vedi soppressione); linguaggio causale e mediazionale forte.

---

## 7. Struttura proposta del manoscritto (in inglese)

1. **Introduction** — atteggiamenti/emozioni → comportamenti in matematica; limiti di self-report e outcome distali (Choe et al., 2019; Cuder et al., 2024; Pizzie & Kraemer, 2017); paradigma di Wang et al. (2025) e Hilz et al. (2025); distinzione tratto/stato (Mammarella et al., 2023); help-seeking come comportamento ambivalente (adattivo vs evitante: letteratura su *instrumental vs executive help-seeking*, Newman, Karabenick — da verificare); obiettivi e RQ.
2. **Method** — partecipanti (N, classi, scuole, consenso, esclusioni, stopping rule vs N raggiunto: 547 > 400 previsto); procedura T0/T1 (figura del disegno come nelle slide); compito (pool di 100 item, 25 s, hint, dwell time 15 s, prompt ogni 5 skip/help); misure con affidabilità; **preregistrazione e deviazioni**; analisi.
3. **Results**
   3.1 Descrittive, psicometria (CFA, affidabilità), correlazioni; distribuzione delle scelte, accuratezza per scelta, casi estremi.
   3.2 Analisi preregistrate (pipeline ordinale completa, sintetica; tabella).
   3.3 Struttura della scelta: ordinale vs `cs()` vs sequenziale vs categoriale (LOO).
   3.4 Modello multinomiale: tutti i contrasti + effetti marginali sulla scala delle probabilità (figura stile slide, ma con predittori in z e bande di incertezza).
   3.5 Robustezza: clustering classe/item, modelli un predittore alla volta (soppressione), modelli binari, trattamento dei NA, multinomiale frequentista.
4. **Discussion** — la misura non è un continuum; help come evitamento mascherato vs supporto adattivo; ruolo di abilità, evitamento, worry; mancati effetti di competenza percepita e buoyancy; limiti (ordine temporale, misure post-compito, accuratezza bassa, singolo compito aritmetico, contesto italiano); implicazioni per la ricerca (misure comportamentali, non assumere ordinalità) e per la didattica (lettura della richiesta d'aiuto).
5. **Supplementary** — codice e dati (OSF), tabelle complete, diagnostica MCMC, tabella delle deviazioni.

Figure chiave (4–5): (1) disegno T0/T1 + schermate del compito; (2) distribuzione scelte e accuratezza per scelta; (3) confronto modelli (Δelpd con SE); (4) effetti marginali per predittore × scelta (un pannello unico a "forest plot" sulla scala delle probabilità); (5) eventualmente probabilità predette per abilità ed evitamento (come le aree impilate delle slide).

Rivista target (da discutere): *Learning and Individual Differences*, *Contemporary Educational Psychology*, *British Journal of Educational Psychology*, *Journal of Experimental Child Psychology*. Il taglio metodologico (struttura della misura) potrebbe anche suggerire *Collabra: Psychology* o *Advances in Methods and Practices in Psychological Science* se si spinge di più su quella parte.

---

## 8. Riorganizzazione del repository (prima di rifare le analisi)

Proposta (da concordare, **non ancora eseguita**):

```
Data/
  raw/DatasetBehaviours.csv          # dati reali (da Analisys_truedata/)
  simulated/SimulatedDatasetBehaviours.csv
Preregistration/
  OSF final preregistration text.pdf
  Template_OSF_Behaviours_020426.pdf
  Script/ ...                        # attuale Script/ (power, simulazione, pipeline preregistrata)
Script/
  01_data_prep.R                     # pulizia, ID duplicato, classe, z-score, codifica scelta
  02_psychometrics.R                 # CFA, affidabilità, descrittive, correlazioni
  03_preregistered_pipeline.R        # pipeline OSF eseguita sui dati reali
  04_structure_models.R              # ordinale / cs / sequenziale / categoriale, LOO
  05_multinomial_main.R              # modello principale + effetti marginali + figure
  06_robustness.R
Results/                             # .rds, tabelle, figure (con nomi parlanti)
Report/
  note.md                            # questo file
  paper.qmd                          # manoscritto (Quarto) — quando si inizia a scrivere
Slide/
```

Altre pulizie: spostare `results5000_*_v2.RData` dalla root; rimuovere/archiviare `CheckDatiFeedback.R` (altro progetto) e `outputs/` (risultati su dati simulati); aggiornare `README.md`; i modelli `.rds` da ~9 MB forse fuori da git (o Git LFS). Valutare se i dati reali possono stare in un repo pubblico (pseudonimizzazione, consenso) prima del deposito su OSF/Zenodo.

---

## 9. Piano operativo (sequenza)

1. **Recuperare da Riccardo** il codice di slide/figure/robustezza (e la heatmap) e chiarire: ID duplicato, NA nelle scelte, codifica accuratezza/timeout, età, classe/scuola, ordine dei questionari.
2. **Riorganizzare il repo** (§8) e scrivere `01_data_prep.R` pulito.
3. **Rieseguire** la pipeline preregistrata sui dati puliti (riproducibilità dei `.rds` esistenti) e il confronto di struttura esteso (§4.12).
4. **Modello principale multinomiale** con predittori standardizzati, `(1|id)` + `(1|class)` + `(1|item)`, tre contrasti, effetti marginali.
5. **Robustezza** (§4.2, 4.3, 4.5, 4.8) e analisi secondarie (accuratezza × scelta, dinamiche nel compito).
6. **Decisioni di framing** con i coautori: numerazione RQ, peso della mediazione, rivista target, primo autore e ordine autori (slide: Pagan, Toffalini, Masi, Caviola).
7. **Scrittura**: prima Method + Results (tabelle e figure definitive), poi Introduction e Discussion, abstract per ultimo. Tabella deviazioni dalla prereg sin dall'inizio.
8. Pacchetto di condivisione (dati + codice + README) su OSF, collegato alla prereg.

### Domande aperte per i coautori

- Le slide sono la versione "approvata" della storia, o c'è margine per rivedere l'interpretazione della buoyancy alla luce della soppressione?
- Va data enfasi metodologica (struttura della misura) o sostantiva (atteggiamenti → scelte)? Influenza rivista e introduzione.
- La mediazione si riporta solo perché preregistrata, o qualcuno ci tiene?
- Esistono dati dello stesso campione usati in altri lavori (es. ansia da test, citata nella prereg)? Va dichiarato per evitare problemi di "salami slicing".
