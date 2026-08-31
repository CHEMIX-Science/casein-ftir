# CHEMIX · Casein FTIR

**Un outil Python pour explorer les spectres infrarouges de la caséine.**

Développé par CHEMIX pour le projet Galalithe, Casein FTIR accompagne l'analyse
FTIR, de l'import du spectre à la comparaison des bandes et à la création d'un
rapport. Il s'adresse aux étudiants, enseignants et chercheurs qui souhaitent
examiner des spectres de caséine et de matériaux associés avec des paramètres
explicites et des résultats exportables.

[English](README.en.md) · [Télécharger](https://github.com/CHEMIX-Science/casein-ftir/releases) · [Méthode et limites](docs/SCIENTIFIC_SCOPE.md)

## Ce que vous pouvez faire

| Besoin | Fonction du programme |
|---|---|
| Préparer un spectre | Corriger la ligne de base, lisser le signal et normaliser son intensité. |
| Examiner les bandes | Repérer des pics dans des fenêtres prédéfinies, intégrer leurs aires et calculer des ratios. |
| Explorer la région amide I | Ajuster six composantes gaussiennes et visualiser leur contribution au signal. |
| Comparer des spectres | Calculer des similarités, repérer des écarts de bandes et exporter une différence numérique. |
| Utiliser une référence théorique | Charger le spectre DFT fourni pour un fragment moléculaire du modèle Galalithe, ou votre propre référence. |
| Exploiter un calcul Gaussian | Convertir une sortie de calcul de fréquences en spectre IR. |
| Présenter une analyse | Produire un rapport texte ou HTML et un graphique de décomposition amide I. |
| Découvrir l'outil | Utiliser les exemples synthétiques fournis ou générer de nouveaux spectres de démonstration. |

Le traitement s'effectue localement, en ligne de commande ou depuis Python.
Les rapports sont enregistrés sur votre ordinateur.

> Casein FTIR est un outil exploratoire : ses résultats aident à examiner un
> spectre, mais ne suffisent pas à identifier un matériau ou à établir sa pureté,
> son innocuité ou son degré de réticulation. Les composantes ajustées et la
> référence théorique doivent être interprétées avec les conditions de mesure
> et d'autres éléments scientifiques.

## Installation

Il faut **Python 3.10 ou plus récent**. Téléchargez et décompressez le code depuis
[les versions publiées](https://github.com/CHEMIX-Science/casein-ftir/releases),
ou utilisez Git :

```bash
git clone https://github.com/CHEMIX-Science/casein-ftir.git
cd casein-ftir
```

Dans le dossier du projet, créez un environnement :

```bash
python -m venv .venv
```

Sous Windows PowerShell :

```powershell
.\.venv\Scripts\python.exe -m pip install "."
.\.venv\Scripts\python.exe -m casein_ftir --help
```

Sous Linux ou macOS :

```bash
source .venv/bin/activate
python -m pip install "."
python -m casein_ftir --help
```

Les commandes suivantes supposent l'environnement activé. Sous PowerShell,
vous pouvez remplacer `python` par `.\.venv\Scripts\python.exe`.
La commande `casein-ftir` est également disponible après installation.

## Première analyse

Cet exemple lit un spectre synthétique de caséine, applique les prétraitements
par défaut, ajuste la région amide I et enregistre un rapport et les courbes ajustées :

```bash
python -m casein_ftir analyze examples/synthetic_casein.csv --deconv --deconv-export deconv.csv --report rapport.html
python examples/plot_amide_deconv.py deconv.csv --output figure.png
```

Ouvrez `rapport.html` dans votre navigateur et `figure.png` pour examiner le
résultat. Remplacez ensuite le chemin du CSV par celui de votre spectre.
Choisissez des noms de sortie distincts : un fichier existant peut être remplacé.

Pour comparer deux spectres :

```bash
python -m casein_ftir compare examples/synthetic_casein.csv examples/synthetic_galalithe.csv --diff-out difference.csv
```

Les exemples `synthetic_*` sont des simulations, pas des mesures expérimentales.
Vous pouvez les régénérer avec `python examples/generate_demo.py`.

## Comparaison avec une référence

Une référence DFT est fournie pour explorer les signatures d'un petit fragment
moléculaire associé au modèle Galalithe. Elle représente un modèle théorique,
pas le spectre complet de la caséine ni un étalon certifié de galalithe.

```bash
python -m casein_ftir analyze examples/synthetic_galalithe.csv --galalithe --use-default-galalithe-ref --no-preprocess-ref --report comparaison_dft.html
```

L'option `--no-preprocess-ref` conserve la référence théorique telle quelle.
Pour consulter sa description, utilisez `python -m casein_ftir info`.
Les paramètres du calcul et leurs limites sont décrits dans
[la documentation des données](docs/DATA_PROVENANCE.md).

Vous pouvez également fournir votre propre référence :

```bash
python -m casein_ftir analyze mon_spectre.csv --reference ma_reference.csv --report comparaison.html
```

Pour convertir votre sortie Gaussian en un CSV utilisable comme référence :

```bash
python -m casein_ftir gaussian convert mon_calcul.log --scale 1.0 --output reference.csv
```

Le programme lit la sortie du calcul ; il ne lance pas Gaussian. `--scale 1.0`
signifie sans correction d'échelle : choisissez un facteur adapté à votre
méthode et à votre base. Les intensités IR doivent être présentes dans le fichier.

## Formats et options

Le format CSV comporte **deux colonnes numériques** : nombre d'onde en cm⁻¹,
puis absorbance. Utilisez un point décimal. Les séparateurs virgule, tabulation,
point-virgule et espaces sont acceptés, ainsi qu'un en-tête textuel.

```csv
wavenumber_cm-1,absorbance
1700,0.10
1698,0.12
1696,0.15
```

Cet extrait illustre le format ; une analyse demande un spectre couvrant les
bandes étudiées. Les doublons et valeurs non finies sont refusés. Convertissez
explicitement la transmittance en absorbance avant d'importer un CSV.

Les fichiers TSV, TXT, DAT, ASC et DPT sont également acceptés.
Pour activer les lecteurs JCAMP-DX et Bruker OPUS :

```bash
python -m pip install ".[jcamp,opus]"
```

Leur compatibilité doit être vérifiée avec les fichiers de votre instrument.
Le format Thermo `.spa` n'est pas pris en charge ; exportez-le en CSV.

Les moteurs optionnels `lmfit` et `pybaselines` s'installent avec :

```bash
python -m pip install ".[deconv,baseline]"
```

SciPy fournit les méthodes de repli. Les résultats peuvent varier selon le
moteur et les paramètres : conservez ces informations avec votre analyse.
L'étalonnage et la relation Beer-Lambert sont accessibles dans l'API Python
`casein_ftir.quantification`, sans commande CLI `calibrate`.

## Documentation et contributions

- [Méthodes, indicateurs et limites scientifiques](docs/SCIENTIFIC_SCOPE.md)
- [Exemples et référence DFT](docs/DATA_PROVENANCE.md)
- [Tests et compatibilité](docs/VALIDATION.md)
- [Références bibliographiques](docs/REFERENCES.md)
- [Contribuer](CONTRIBUTING.md) · [Signaler un problème](https://github.com/CHEMIX-Science/casein-ftir/issues) · [Sécurité](SECURITY.md)

Pour lancer les tests depuis le dossier du projet :

```bash
python -m pip install -e ".[test]"
python -m pytest
```

## Licence

Le code, la documentation, les exemples et la référence DFT fournie sont
accessibles sous [licence MIT](LICENSE). Voir [les attributions](LICENSING.md).

Développé par **[CHEMIX](https://chemix-paris.com)** pour le projet **Galalithe**.
