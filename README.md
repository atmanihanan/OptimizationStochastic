# Projet d'implémentation de l'algorithme de Newton avec région de confiance en Scilab

## Description du projet

Ce projet consiste à implémenter l'algorithme de Newton avec région de confiance (TR Nwt) pour résoudre des problèmes d'optimisation unidimensionnels. L'algorithme est appliqué à des fonctions qui sont au moins deux fois dérivables. L'objectif principal est de trouver les zéros de la dérivée d'une fonction donnée, en affichant les valeurs pertinentes à chaque itération.

## Prérequis

Avant de commencer, assurez-vous que le **toolbox diffcode** est installé dans Scilab pour utiliser les dérivées des fonctions. 

## Fichiers inclus

- `Fcts et Drvs.sci`: Ce script fournit les dérivées des fonctions en utilisant le toolbox diffcode.
- `TR Nwt.sci`: Ce script doit être complété avec le code de l'algorithme de Newton avec région de confiance.
- `TestTP1.sce`: Un script de test qui affiche la fonction et son approximation quadratique (polynôme de Taylor).
- `Fct Exemple.sci`: Un exemple de fonction à optimiser.

## Algorithme

L'algorithme de Newton avec région de confiance est implémenté selon les étapes suivantes :

1. Initialisation des données : valeur de départ `x`, taille de région de confiance `∆`, critère d'arrêt `ϵ`, et fonction `f`.
2. Calcul de la forme quadratique `q(d)` à partir de la fonction et de ses dérivées.
3. Exécution de la boucle principale :
   - Détermination de la direction `dR` basée sur la comparaison des valeurs de `q(∆)` et `q(−∆)`.
   - Calcul de la direction de Newton `dN`.
   - Mise à jour de la direction de recherche `dR` si les conditions sont satisfaites.
   - Calcul des valeurs `ared` (gain réel) et `pred` (gain prédit).
   - Mise à jour de la taille de la région de confiance `∆` basée sur le rapport `r`.
4. Répétition jusqu'à ce que le critère d'arrêt soit atteint (|f'(x)| < ϵ).

## Résultats

Chaque itération affiche les valeurs suivantes :
- `xk`: la valeur courante de x
- `f'(xk)`: la dérivée de f en xk
- `∆k`: la taille de la région de confiance
- `ared`: le gain réel
- `pred`: le gain prédit

Les résultats peuvent être comparés aux valeurs de référence fournies dans le tableau ci-dessous :

| Iter | x         | |f'(x)|   | ∆   | predared       |
|------|-----------|---------|------|-----------------|
| 0    | 2.0      | 0.08284 | 1    | -0.142467       |
| 1    | 1.0      | 0.25    | 2    | -0.25          |
| 2    | 1.0      | 0.25    | 1    | -0.1875        |
| 3    | 1.0      | 0.25    | 0.5  | -0.109375      |
| 4    | 0.5      | 0.09467 | 0.5  | -0.0050441     |
| 5    | 0.606557 | 0.0064  | 1    | -0.0000210     |
| 6    | 0.599998 | 0.0000017| 2   | -1.517 x 10^-12 |
| 7    | 0.6      | 0       | 4    | -0.1730769     |

## Exécution

Pour exécuter le projet, suivez les étapes ci-dessous :
1. Complétez le fichier `TR Nwt.sci` avec le code de l'algorithme.
2. Exécutez le fichier `TestTP1.sce` pour afficher la fonction et son approximation.
3. Comparez les résultats avec ceux obtenus en utilisant la commande `optim` de Scilab.

## Acknowledgments

Ce projet s'inspire des travaux pratiques de l'algorithme de Newton et de l'optimisation, et utilise les ressources de Scilab pour l'implémentation.

## License

Ce projet est sous licence MIT - voir le fichier [LICENSE](LICENSE) pour les détails.

