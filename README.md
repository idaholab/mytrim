MyTRIM [![Build Status](https://travis-ci.org/idaholab/mytrim.svg?branch=master)](https://travis-ci.org/idaholab/mytrim)
======

Three dimensional binary collision Monte Carlo library for sampling ion
collision cascades in materials.

The math behind this is based on the book
*"The Stopping and Range of Ions in Solids"* New York: Pergamon Press.
By J. F. Ziegler, J. P. Biersack, and U. Littmark (1985 (new edition in 1996)).

[Wikipedia article](https://en.wikipedia.org/wiki/Stopping_and_Range_of_Ions_in_Matter)

![TRIM Schema](http://idaholab.github.io/img/mytrim/trims_hor.png)

*MyTRIM* is fully tree dimensional and extendable to arbitrary sample geometries.

Examples
-----

### Energy deposition
![TRIM Schema](http://idaholab.github.io/img/mytrim/trim.png)

Displacement events calculated by *MyTRIM*, colored by inelastic energy loss. The cube is 100nm on the side.

### Gas resolution

The video linked below shows ion cascades in a ceramic material knocking gas atoms out of a bubble inclusion.

[MyTRIM cascade video](http://idaholab.github.io/img/mytrim/3dtrim.mp4)

JSONCPP
----
The `runmytrim` executable needs the jsoncpp library compile it in teh jsoncpp directory with

```
mkdir build && cd build
cmake -DBUILD_STATIC_LIBS=ON -DBUILD_SHARED_LIBS=OFF  -G "Unix Makefiles" ..
sudo make install
```


About
----

The MyTRIM source code is licensed under the LGPL 2.1 license.
The data files supplied in data/ are prepared by James. F. Ziegler as part of the
[SRIM](http://www.srim.org) distribution.

This software is provided "AS IS" without any express or implied warranty.
