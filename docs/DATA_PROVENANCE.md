# Provenance et séparation des données

| Élément distribué | Provenance | Statut et portée |
|---|---|---|
| `examples/synthetic_casein.csv` | `examples/generate_demo.py`, simulateur fourni | **CONFIRMÉ** : données synthétiques, pas de mesure expérimentale |
| `examples/synthetic_galalithe.csv` | Même générateur, paramètres explicites | **PROPOSITION** : scénario illustratif, pas de résultat matériel |
| Tests Python | Jeux synthétiques et courts blocs Gaussian artificiels | **CONFIRMÉ** : validation logicielle seulement |

**CONFIRMÉ — données exclues :** la référence DFT interne, les résultats
expérimentaux, rapports d'échantillons, graphiques historiques, anciens CSV
non qualifiés, journaux réels de calcul, entrées Gaussian historiques, scripts
de campagne, archives d'origine et environnements Python ne sont pas distribués.
L'historique Git de préparation contenant le CSV DFT n'est pas importé.

**À VALIDER** — Toute publication ultérieure de la référence DFT interne exige
une autorisation explicite de redistribution de ce fichier. La validation MIT
du code ne vaut pas validation de ses droits ni de sa provenance scientifique.

**CONFIRMÉ** — L'import de références autorisées reste disponible avec
`--reference` ou `load_reference()`. La conversion Gaussian reste disponible
pour les fichiers fournis par l'utilisateur. Les anciens accesseurs à une
référence intégrée expliquent son absence ; aucun faux spectre DFT n'est substitué.

**CONFIRMÉ** — Les CSV conservent le type synthétique/DFT reconnu. Les chemins
absolus ne sont plus exportés par défaut dans les rapports et CSV usuels.
Les noms de fichiers et d'échantillons doivent toujours être relus avant partage.
