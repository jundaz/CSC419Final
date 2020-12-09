# CSC419 libigl project-Bihamonic
A libigl style implementation of Biharmonic distance method by Y.lipman et.al <br />
link to paper: https://www.cs.princeton.edu/~funk/biharmonic.pdf

biharmonic_distance

biharmonic_distance(V: array, F: array, approach: int, k: int)

Constructs a matrix of biharmonic distances of every pair of vertices given by V and F.
We implement the exact approach according to the 3.2 Discrete construction section in the paper
and the approximate approach according to the 3.3 Practical computation section in the paper.
The output can be visualize by passing the returned matrix D to libigl viewer and there will be
color change appearing on the input surface: purple means the closest distance and yellow means
the farthest distance. According to our results, the approximate approach is usually 10 times faster
than the exact approach when taking the first 10 eigen vector of generalized eigenvalue problem L.phi = A.lambda.phi
where lambda is the eigen value and phi is the eigen vector

Remark: when using approximate approach, some large input surfaces might
have undesirable performance due to the eigen decomposition not converging.

we have included some triangle mesh file in the same directory with executable file
after compiling you can simply call <br />./biharmonic <name_of_mesh> <approach you want to use(exact 0/approx1)> <number of eigen vectors used to compute approx distance><br />
to run the function.

example:<br />
.\biharmonic.exe .\bunny.off 0 10

dependency:<br />
requires libigl, eigen and spectra library

How to compile:<br />
The procedure is the same with all other assignments we did

Parameters:<br />
    V: #V by 3 list of mesh vertex positions<br />
    F: #F by 3 list of mesh face indices into V<br />
    approach: 0 for exact approach and 1 for approximate approach<br />
    k: first k eigen vectors to keep<br />

Returns:<br />
    D: #V by #V matrix of biharmonic distances where an entry (i,j) is the distance between vertex i and j

Example Output Speed:<br />
    bunny.off: exact approach is 1.56345s, approximate approach with k = 10 is 0.27575s<br />
    cactus.obj: exact approach is 3.64s, approximate approach with k = 10 is 0.32s<br />
