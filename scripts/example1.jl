using MatrixPolynomials

Jcol = colBlock([1,2,3])
Jrow = rowBlock([1,2,3])
Jinf = infBlock([1,2,3])

J2 = eigBlock(2, [1,2])
J3 = eigBlock(3, [2,2])

J = blkdiag([Jcol, Jrow, Jinf, J2, J3])
