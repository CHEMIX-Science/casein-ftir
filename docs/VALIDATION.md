# Validation de la distribution publique

**CONFIRMÉ** — La suite utilise des données synthétiques. Aucun test ne dépend
de la référence DFT interne ; un contrôle vérifie explicitement son absence.
L'import/export du type DFT est testé à l'aide d'un objet synthétique portant
une métadonnée de test, sans simuler un résultat réel de calcul quantique.

**CONFIRMÉ** — 201 tests réussis sous Windows/Python 3.12.3 avec lmfit et
pybaselines le 31 août 2026 (18,82 s). La roue et l’archive source sont construites
et leur contenu est contrôlé : licence MIT présente, CSV DFT interne absent.

**CONFIRMÉ** — La configuration GitHub Actions couvre Linux/Python 3.10,
Windows/Python 3.12, macOS/Python 3.13 et Linux/Python 3.12 avec les moteurs
optionnels. Les résultats de chaque exécution sont disponibles dans
[GitHub Actions](https://github.com/CHEMIX-Science/casein-ftir/actions).

**À VALIDER** — Compatibilité des lecteurs JCAMP/OPUS avec les instruments réels
et validité analytique sur un jeu de données indépendant.
Les tests ne prouvent pas la validité analytique sur des mesures expérimentales.
