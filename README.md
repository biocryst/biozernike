# BioZernike
Protein structure descriptors and alignment based on 3D Zernike moments.

See it in action: http://shape.rcsb.org/


This library implements 3D Zernike moment calculation and normalization as introduced
in [Canterakis 1996](https://link.springer.com/chapter/10.1007/978-3-642-80294-2_36) and
 [Canterakis 1999](https://lmb.informatik.uni-freiburg.de/people/canterakis/publications/cant_scia99.pdf).
Routines are provided for calculation of:
* Trivial rotational invariants (norms of the vectors Ωnl), 
[commonly referred](http://www.eurekaselect.com/88710/article) as 3D Zernike Descriptors. 
Calculation of these descriptors is based on the 3D Zernike Moments library by Marcin Novotni (see [Novotni and Klein 2003](https://cg.cs.uni-bonn.de/aigaion2root/attachments/novotni-2003-3d.pdf)).
The implementation here fixes a bug that causes the invariants of the same order to be cumulative.
* Complete rotational invariants (Canterakis norms), not available in the Novotni library.
* Alignments, based on the complete rotational invariants. 

See the [preprint](https://www.biorxiv.org/content/biorxiv/early/2019/11/16/845123.full.pdf) describing this work.

## License
The `zernike` package is derived from the ["3D Zernike Moments" library](http://www.cg.cs.uni-bonn.de/project-pages/3dsearch/) by Marcin
Novotni  and is distributed under the terms of LGPL v2.0.

The `volume` package is derived from the ["gmconvert" program](https://pdbj.org/gmfit/) by Takeshi Kawabata
 and is distributed under the terms of LGPL v3.0.

The `complex` package is derived from the [code by John B. Matthews](https://sites.google.com/site/drjohnbmatthews/polyroots) and is distributed 
under the terms of LGPL v3.0.

The BioZernike library as a whole is distributed under the terms of the MIT License.
