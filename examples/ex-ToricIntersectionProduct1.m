Q := [ [1,1,1,1,1,1,1],
       [0,0,0,0,1,1,1],
       [0,0,1,1,0,0,1] ];
mu := [ [2], [2], [2] ];
u := [1];
K := [1];
D := mu cat [K, K, K];
ToricIntersectionProduct(Q, u, D : T := [2,2]);
// output: 2
