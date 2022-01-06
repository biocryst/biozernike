# BioZernike
Protein structure descriptors and alignment based on 3D Zernike moments.

See it in action:

- https://www.rcsb.org : assembly and chain search integrated with other types of searches (text, sequence etc)

- http://shape.rcsb.org : standalone frontend application that performs assembly and chain search and calculates alignments
 on the fly, displaying them with NGL  


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

The test directory contains tests that demonstrate how to read PDB-deposited protein structures (with the help of [BioJava](https://github.com/biojava/biojava))
and perform Zernike moment invariant calculation and alignment.

See the publication describing this work: [Real time structural search of the Protein Data Bank. Guzenko D, Burley SK, Duarte JM. PLoS Computational Biology 2020](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007970).

## Use it from your project
We publish jar artifacts to maven central. In a maven project, you can use this library by adding this dependency:
```xml
    <dependencies>
      <dependency>
        <groupId>org.rcsb</groupId>
        <artifactId>biozernike</artifactId>
        <version>1.0.0-alpha10</version>
      </dependency>
    </dependencies>
```

## License
The `zernike` package is derived from the ["3D Zernike Moments" library](http://www.cg.cs.uni-bonn.de/project-pages/3dsearch/) by Marcin
Novotni  and is distributed under the terms of LGPL v2.0.

The `volume` package is derived from the ["gmconvert" program](https://pdbj.org/gmfit/) by Takeshi Kawabata
 and is distributed under the terms of LGPL v3.0.

The `complex` package is derived from the [code by John B. Matthews](https://sites.google.com/site/drjohnbmatthews/polyroots) and is distributed 
under the terms of LGPL v3.0.

The BioZernike library as a whole is distributed under the terms of the MIT License.
