#Assemblathon Codebase

This is a set of evaluation scripts developed by [Dent Earl](https://github.com/dentearl/) and [Benedict Paten](https://github.com/benedictpaten/) to assess the [Assemblathon](http://assemblathon.org/) 1 dataset.

##Dependencies

* [jobTree](https://github.com/benedictpaten/jobTree)
* [sonLib](https://github.com/benedictpaten/sonLib)
* [cactus](https://github.com/benedictpaten/cactus)
* [cactusTools](https://github.com/benedictpaten/cactusTools)
* [assemblaLib](https://github.com/benedictpaten/assemblaLib)

##Running Tests
Download the data set [here](http://compbio.soe.ucsc.edu/assemblathon1/data.tar.gz), expand it and place the data/ directory in the same directory as this README. Then type

    make testLittle

and off you go. There are larger and longer running tests inside of the Makefile if you'd like to justify spending a couple hours doing something else.
