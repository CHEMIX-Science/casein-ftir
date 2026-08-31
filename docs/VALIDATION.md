# Tests et compatibilité

Les tests vérifient les imports, les prétraitements, les calculs numériques,
les exports et les commandes sur des jeux synthétiques. Des contrôles dédiés
vérifient également l'intégrité du CSV DFT fourni, son chargement et le maintien
des avertissements propres à une référence de fragment.

La configuration GitHub Actions couvre Linux/Python 3.10, Windows/Python 3.12,
macOS/Python 3.13 et Linux/Python 3.12 avec les moteurs optionnels `lmfit` et
`pybaselines`. Chaque configuration construit le paquet puis l'installe dans
un environnement séparé pour exercer les commandes et la référence intégrée.

Consultez [les résultats des tests automatiques](https://github.com/CHEMIX-Science/casein-ftir/actions)
pour la version qui vous intéresse. Pour exécuter la suite localement :

```bash
python -m pip install -e ".[test]"
python -m pytest
```

Ces tests contrôlent le logiciel, pas la validité d'une méthode analytique sur
des mesures réelles. Les lecteurs JCAMP/OPUS doivent être vérifiés avec les
fichiers de votre instrument. Voir [les limites scientifiques](SCIENTIFIC_SCOPE.md).
