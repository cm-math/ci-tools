Q := [[1,1,1,0,0,-3], [0,0,0,1,1,1], [0,1,2,1,2,0]];
mu := [0,3,0];
SearchPrimeBinomial(Q, mu, 1 : T := [3]);
/* output:
true [
    [ 0, 0, 0, 0, 3, 0 ],
    [ 0, 4, 2, 1, 0, 2 ]
]
*/

SearchPrimeBinomial(Q, mu, 6 : T := [3]);
// output: false
