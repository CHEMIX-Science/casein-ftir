# Portée scientifique et niveau de preuve

## Ce que le code calcule

- Sommes de bandes paramétriques pour les démonstrations, avec graines fixes de bruit.
- Correction de ligne de base, lissage et normalisation ; résultats dépendants des paramètres.
- Aires intégrées dans des fenêtres fixes. Une ligne de base locale linéaire est soustraite par défaut, puis les valeurs négatives sont tronquées à zéro : ce choix modifie les aires.
- Ratios entre ces aires ; un dénominateur quasi nul produit `NaN`.
- Ajustement amide I à six gaussiennes contraintes (lmfit : least_squares, au plus 5 000 évaluations ; repli SciPy : curve_fit, maxfev 20 000) ; pourcentages calculés sur les aires des gaussiennes complètes (et non sur leur seule portion dans la fenêtre).
- Indicateurs de similarité : Pearson, cosinus, RMSE et `1 - aire(|A-B|)/aire((|A|+|B|)/2)`. Le dernier indice peut être négatif ; ce n'est pas une probabilité.
- Catégories `excellent/good/fair/poor` dérivées de seuils Pearson 0,95/0,85/0,70 ; les clés historiques sont conservées pour compatibilité, sans jugement de qualité du matériau.

## Hypothèses du modèle

Les centres, largeurs, intensités relatives, fenêtres et décalages entre espèces
sont des paramètres de modèle. Les six catégories amide I sont
des étiquettes de modèle. Les noms historiques `crosslink_degree`,
`crosslink_index` et `residual_formol` ne sont pas des grandeurs expérimentales
étalonnées. Leur conservation évite de casser les scripts existants.

La référence DFT fournie modélise un fragment moléculaire ; elle n’est pas un
étalon certifié du matériau. La même prudence s’applique à vos propres références. Une coïncidence de bandes doit
être interprétée comme un indice à examiner, jamais comme une preuve exclusive.

## Limites et validations restantes

- Attribution et spécificité des bandes sur chaque échantillon, notamment les fenêtres de marqueurs Galalithe.
- Répétabilité, incertitudes, effets de l'eau, préparation, géométrie ATR/transmission, normalisation et réglages instrumentaux.
- Conversion des fractions d'aires en populations structurales ; ambiguïté des ajustements à plusieurs bandes et dépendance aux contraintes.
- Choix d'une référence appropriée, de son facteur d'échelle et de sa largeur de convolution.
- Comparabilité des fenêtres disponibles : ne pas interpréter une fenêtre non mesurée comme une absence de bande. Utiliser des spectres couvrant intégralement les fenêtres analysées.
- Lissage et ALS supposent un échantillonnage suffisamment régulier ; aucun rééchantillonnage automatique n'est appliqué à l'entrée.
- Intensités et facteur d'absorption adaptés à l'étalonnage Beer–Lambert. Une aire intégrée nécessite un coefficient étalonné pour cette aire ; ne pas lui appliquer sans adaptation un coefficient de hauteur de pic. Ne pas normaliser l'amplitude avant un dosage absolu.
- Le lecteur Gaussian rassemble les blocs présents et ne certifie ni convergence, ni minimum stable, ni identité chimique. Fournir un fichier ne contenant qu'un calcul de fréquences pertinent et contrôler le journal original.
- Aucun test externe sur jeu expérimental indépendant n'est revendiqué. Les tests numériques ne constituent pas une validation analytique.

## Bibliographie

Les références sont conservées dans `database.py` pour la
traçabilité des modèles. L'existence d'un article ne prouve pas que chaque
nombre du modèle en est extrait. Voir le [registre](REFERENCES.md) pour les
contrôles et les éléments non vérifiés.
