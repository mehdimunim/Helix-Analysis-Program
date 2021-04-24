*Adapted from simulaid.html*

## Perform helix analyses with a selected set of helices for each helix 
* calculate its axis
* find and plot the angles between the axis and the laboratory frame coordinate axes
* calculate the and plot radius if circle fitted to the alpha carbons (a bend indicator) and provide a novel shape estimator.
* calculate and plot the length of the helix
* calculate and plot the average turn angle per residue.
*when more than one helix is selected for analysis, the angles between each pair and the distance between the centers of each helix will also be calculated*
## When scanning a trajectory, for each frame
  * calculate and plot the angle between the axis in the current and in the reference structure (read in at the start of the run)
  * calculate and plot the angle the helix was rotated around the axis (again, compared to the reference structure). The overall rotation of the protein can be applied to the    reference helix before the comparison.
  * calculate and plot the displacement fo the helix center compared to the reference structure. The displacement of the protein's COM can be subtracted.
  * calculate the angle between the normals of the plane fitted to the reference structure and the current structure, rotate both normals so that the reference planes normal is parallel to the Z axids and plot the progress of the X-Y projection of the current (rotated) normal.

*At the end of the trajectory scan the correlation coefficients among all pairs of calculated properties will be calculated. For properties involving two angles the Fisher-Lee circular correlation coefficient will also be calculated: For properties involving one angle only the Mardia linear-circular correlation coefficient will also be calculated:*

C(a,x) = (r122+r122-2*r12*r13*r23)/(1-r232)

where

r12 = corr(sin(a),x);   r13 = corr(cos(a),x);   r23 = corr(sin(a),cos(a))
