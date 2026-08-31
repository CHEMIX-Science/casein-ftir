# CHEMIX · Casein FTIR

**Explorer des spectres. Comparer des bandes. Documenter les limites.**

[English](README.en.md) · [Méthode et limites](docs/SCIENTIFIC_SCOPE.md) · [Données](docs/DATA_PROVENANCE.md) · [Installation](#installation)

**CONFIRMÉ** — Programme Python développé au sein de CHEMIX pour le projet
Galalithe, à partir du programme interne `casein_ftir`.
Cette distribution fournit une version autonome du logiciel sous licence MIT.

**PROPOSITION** — Un outil de recherche exploratoire et de médiation : il aide
à prétraiter, comparer et décomposer des spectres FTIR de caséine et de matériaux
associés. Il ne remplace pas l'interprétation scientifique des mesures.

> **À VALIDER** — Les attributions, seuils et modèles ne constituent pas une
> méthode analytique validée. Une similarité élevée ne prouve ni pureté, ni
> innocuité, ni réticulation. Les fractions d'aires ajustées ne sont pas des
> populations structurales démontrées. `residual_formol` n'est pas un dosage du
> formaldéhyde et `crosslink_index` n'est pas un taux de conversion.

## Fonctions disponibles

**CONFIRMÉ — fonctions présentes dans le code :**

| Fonction | Accès |
|---|---|
| Lecture CSV, TSV, TXT, DAT, ASC, DPT | Bibliothèque et CLI |
| Lecture JCAMP-DX et Bruker OPUS | Extensions optionnelles ; compatibilité instrumentale à vérifier |
| Ligne de base ALS, polynomiale ou rubberband | `analyze` |
| Lissage Savitzky–Golay et normalisation | `analyze` |
| Détection de bandes, aires et ratios | `analyze` |
| Ajustement gaussien de la région amide I | `deconv` ou `analyze --deconv` |
| Comparaison, référence et différences | `compare`, `analyze --reference` |
| Spectres synthétiques reproductibles | `simulate` |
| Conversion de sorties Gaussian en spectres | `gaussian` |
| Rapport texte ou HTML local | `analyze --report` |
| Étalonnage et relation Beer–Lambert | API Python `casein_ftir.quantification` |

**CONFIRMÉ** — La lecture Thermo `.spa` et la commande CLI `calibrate` ne sont
pas implémentées. Exporter les données instrumentales en CSV et utiliser l'API
Python pour l'étalonnage. Le logiciel ne lance pas Gaussian.

## Dépôt et téléchargements

[Code source](https://github.com/CHEMIX-Science/casein-ftir) · [Versions](https://github.com/CHEMIX-Science/casein-ftir/releases) · [Tests automatiques](https://github.com/CHEMIX-Science/casein-ftir/actions) · [Signaler un problème](https://github.com/CHEMIX-Science/casein-ftir/issues)

Pour récupérer le projet avec Git :

```bash
git clone https://github.com/CHEMIX-Science/casein-ftir.git
cd casein-ftir
```

## Installation

**CONFIRMÉ — configuration de cette distribution :** Python 3.10 ou plus récent.
Depuis le dossier du dépôt téléchargé :

```bash
python -m venv .venv
```

Windows PowerShell :

```powershell
.\.venv\Scripts\python.exe -m pip install "."
.\.venv\Scripts\casein-ftir.exe --help
```

Linux / macOS :

```bash
source .venv/bin/activate
python -m pip install "."
casein-ftir --help
```

Les commandes ci-dessous supposent l'environnement activé. Sous PowerShell,
on peut remplacer `python` par `.\.venv\Scripts\python.exe`, sans modifier la
politique d'exécution de Windows. `python -m casein_ftir` équivaut à `casein-ftir`.

Extensions facultatives :

```bash
python -m pip install ".[deconv,baseline]"
python -m pip install ".[jcamp,opus]"
```

**CONFIRMÉ** — SciPy fournit les méthodes de repli pour ALS et l'ajustement amide I.
`lmfit` et `pybaselines` sont facultatifs ; leur présence peut modifier les résultats
numériques. Noter les versions et paramètres lors de toute analyse.

## Démonstration sans données expérimentales

**CONFIRMÉ** — Les CSV `synthetic_*` fournis sont générés par le logiciel et ne
sont pas des mesures CHEMIX. Pour les régénérer :

```bash
python examples/generate_demo.py
```

Analyse avec ajustement et rapport local :

```bash
python -m casein_ftir analyze examples/synthetic_casein.csv --deconv --deconv-export deconv.csv --report rapport.html
python examples/plot_amide_deconv.py deconv.csv --output figure.png
```

Comparaison numérique :

```bash
python -m casein_ftir compare examples/synthetic_casein.csv examples/synthetic_galalithe.csv --diff-out difference.csv
python -m casein_ftir analyze examples/synthetic_galalithe.csv --galalithe --reference examples/synthetic_casein.csv --report comparaison.html
```

**PROPOSITION** — Ces exemples illustrent le fonctionnement ; ils ne valident
pas les hypothèses chimiques qui ont servi à générer les spectres.
Les rapports restent sur votre ordinateur. Relire les noms d'échantillons et les
fichiers exportés avant tout partage. Les fichiers de sortie existants peuvent
être remplacés : utiliser des noms distincts pour chaque analyse.

## Format des mesures

**CONFIRMÉ — contrat d'entrée :** deux colonnes, nombre d'onde en cm⁻¹ puis
absorbance ; point décimal ; virgule, tabulation, point-virgule ou espaces comme
séparateur. Un en-tête textuel est admis. Pas de `%T` dans un CSV : convertir
explicitement en absorbance avant importation.

```csv
wavenumber_cm-1,absorbance
1700,0.10
1698,0.12
1696,0.15
```

Cet extrait illustre seulement le format : il est trop court pour l'analyse
amide I. Les points sont triés par nombre d'onde décroissant ; doublons, valeurs
non finies et lignes de mesure malformées sont refusés. Le lissage demande un axe
régulier et au moins 11 points par défaut. L'ajustement à six composantes demande
au moins 19 points dans la fenêtre sélectionnée.

## Références personnelles et calculs Gaussian

**CONFIRMÉ** — La distribution publique ne contient pas le CSV DFT interne.
Vous pouvez comparer un spectre à une référence dont vous avez les droits,
ou convertir votre propre sortie Gaussian.

```bash
python -m casein_ftir analyze mon_spectre.csv --reference ma_reference.csv --report comparaison.html
python -m casein_ftir gaussian convert mon_calcul.log --scale 1.0 --output reference.csv
python -m casein_ftir analyze mon_spectre.csv --reference reference.csv --no-preprocess-ref --report analyse.html
```

**PROPOSITION** — `--scale 1.0` signifie « sans correction » et sert d'exemple,
ce n'est pas un facteur recommandé. Vérifier le facteur adapté à la méthode et
à la base employées. La conversion refuse des intensités IR manquantes.
L'ancienne option `--use-default-galalithe-ref` explique l'absence du fichier ;
utiliser `--reference` pour fournir votre référence autorisée.

## Développement et tests

```bash
python -m pip install -e ".[test]"
python -m pytest
python -m pip install build
python -m build
```

**CONFIRMÉ** — 201 tests passent localement sur la distribution sans données internes. Les contrôles et leurs limites sont consignés dans [VALIDATION.md](docs/VALIDATION.md).

**CONFIRMÉ** — La configuration CI inclut Windows, Linux et macOS, et vérifie
l’installation d’une roue Python. Consulter les résultats des exécutions dans
[GitHub Actions](https://github.com/CHEMIX-Science/casein-ftir/actions).
Les tests synthétiques contrôlent le logiciel ; ils ne prouvent pas sa validité
sur des mesures réelles.

## Attribution, licence et contributions

**CONFIRMÉ** — Projet développé chez CHEMIX et distribué sous [licence MIT](LICENSE),
validée par la présidence le 31 août 2026. Attribution collective : CHEMIX contributors.
Le périmètre et les licences des dépendances sont précisés dans [LICENSING.md](LICENSING.md).
Les limites scientifiques restent applicables indépendamment de la licence.

[Contribuer](CONTRIBUTING.md) · [Sécurité](SECURITY.md) · [Historique](CHANGELOG.md) · [CHEMIX](https://chemix-paris.com)
