#Installing assemblathon.

Tip: place all the downloaded projects in the same root directory, this way no makefile variables need be edited and the python path needs to be pointed at this base directory.

1. Download and install [assemblaLib](https://github.com/benedictpaten/assemblaLib) (this will involve installing [cactus](https://github.com/benedictpaten/cactus) and [sonLib](https://github.com/benedictpaten/sonLib)).
2. Download and install [jobTree](https://github.com/benedictpaten/jobTree)
3. Place the directory containing assemblathon on your python path environmental variable, i.e.

    <code>PYTHONPATH=${PYTHONPATH}:FOO</code>

where <code>FOO/assemblaLib</code> is the path to the parent directory of assemblaLib. 
4. Compile the C code:
Modify the <code>rootPath</code> variable in the src/Makefile to point at where you installed cactus.
In the base directory type '<code>make all</code>' 

See https://github.com/benedictpaten/ for different projects.
