# Exemples et référence théorique

## Spectres de démonstration

`examples/synthetic_casein.csv` et `examples/synthetic_galalithe.csv` sont générés
par `examples/generate_demo.py`. Ils servent à découvrir les commandes et à
examiner le comportement du logiciel. Ils ne représentent pas des mesures
expérimentales ni une validation des hypothèses chimiques du simulateur.

## Référence DFT Galalithe

Le fichier `casein_ftir/data/galalithe_dft_reference.csv` fournit un spectre
théorique pour explorer les signatures d'un fragment moléculaire associé au
modèle Galalithe. Il contient 1 801 points de 4 000 à 400 cm⁻¹, espacés de
2 cm⁻¹, avec une intensité normalisée à un maximum de 1.

Cette référence est fournie par CHEMIX et distribuée sous [licence MIT](../LICENSE).
Ses valeurs numériques sont conservées sans modification.

### Paramètres associés au fichier

| Paramètre | Valeur consignée |
|---|---|
| Fragment | N,N'-diéthylméthanediamine, C5H14N2 |
| Modèle | Pont aminal minimal associé à la réticulation caséine-formaldéhyde |
| Méthode et base | B3LYP/6-311+G(2d,p) |
| Solvant | Eau, IEFPCM |
| Nombre d'atomes / modes | 21 / 57 |
| Facteur d'échelle | 0,9679 |

Ces paramètres sont ceux associés au spectre fourni. Le calcul n'a pas été
reproduit indépendamment pour cette distribution ; la traçabilité complète des
paramètres et du facteur d'échelle reste à établir. Le CSV seul ne permet pas
de vérifier la géométrie, la convergence ou le traitement complet des fréquences.

Cette référence n'est ni le spectre complet de la caséine ni un étalon certifié
de galalithe. Une coïncidence de bandes ne démontre pas l'identité d'une liaison
ou la réticulation du matériau.

### Utilisation

```python
from casein_ftir.defaults import get_galalithe_dft_reference

reference = get_galalithe_dft_reference()
```

En ligne de commande, utilisez `analyze --use-default-galalithe-ref` et
`--no-preprocess-ref` pour conserver le spectre théorique sans prétraitement.
Le marquage `dft_packaged_reference` adapte les messages de comparaison à une
référence de fragment. Vous pouvez aussi charger une référence avec
`--reference` ou `load_reference()`.

## Partager vos résultats

Les rapports et exports CSV usuels n'ajoutent pas le chemin absolu du fichier
source. Les CSV conservent le type synthétique ou DFT reconnu. Relisez toutefois
les noms d'échantillons et les données exportées avant de les partager ; ne
redistribuez que les fichiers pour lesquels vous disposez des droits nécessaires.
