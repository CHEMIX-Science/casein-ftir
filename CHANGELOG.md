# Historique

## 0.2.0 — 31 août 2026

**CONFIRMÉ — changements par rapport au dossier interne fourni :**

- Distribution Python moderne et commande `python -m casein_ftir`.
- NumPy minimal 2.0, cohérent avec l'usage de `numpy.trapezoid`.
- Guide français, introduction anglaise, provenance, périmètre scientifique et règles de contribution.
- Exemples synthétiques régénérables ; résultats expérimentaux et environnements locaux exclus.
- Validation des entrées numériques, tri de l'axe, refus des doublons, lissage irrégulier refusé et comparaison limitée au domaine commun.
- Échecs CLI explicites, refus des options incohérentes et des conversions Gaussian sans intensités valides.
- Échappement HTML des textes et suppression des chemins absolus des exports usuels.
- Préservation du type synthétique/DFT lors d'un export CSV.
- Suppression des conclusions de pureté, aptitude à l'usage et confirmation de réticulation dans les rapports.
- Initialisation des aires lmfit corrigée et solveur least_squares borné à 5 000 évaluations.
- Graphique : suppression de la double inversion de l’axe partagé et rappel des limites des fractions ajustées.
- Conservation des principales clés de l’API ; référence numérique interne exclue, avec message explicite des anciens accesseurs et import de références externes maintenu.

**CONFIRMÉ** — Licence MIT validée le 31 août 2026 et intégrée à la distribution ;
attribution collective CHEMIX contributors, sans noms individuels inventés.
**CONFIRMÉ** — Organisation GitHub CHEMIX-Science créée ; distribution destinée au dépôt `CHEMIX-Science/casein-ftir`.

## 0.1.0 — Version interne fournie

**CONFIRMÉ** — Code et tests présents dans le dossier d'origine ; aucun historique
Git antérieur importé ni résultat de test historique repris sans vérification.

**CONFIRMÉ** — La publication GitHub exclut le CSV DFT interne et tout historique Git qui le contient. Une autorisation explicite de redistribution sera nécessaire pour l’ajouter.
