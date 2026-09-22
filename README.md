# Projet MS14 — Partie 2

## Compression d'image par triangulation de Delaunay

Cette partie porte sur la compression d'une image en niveaux de gris à l'aide de la triangulation de Delaunay et de l'interpolation.

Le projet comprend :
- L'implémentation d'un noyau de Delaunay avec la méthode **Bowyer-Watson**.
- La localisation des points, la création de cavités et la re-triangulation.
- L'application de la triangulation à la compression d'images.
- La reconstruction de l'image par interpolation et l'évaluation avec le **PSNR**.
- Une comparaison entre une méthode de compression naïve et une méthode basée sur un seuil d'erreur d'interpolation.
