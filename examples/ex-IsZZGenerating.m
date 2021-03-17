u1 := [1,0,0];
u2 := [0,1,0];
u3 := [3,3,3];
u4 := [4,4,4];
IsZZGenerating([u1,u2,u3,u4]);
// output: true

w1 := [1,2,1];
w2 := [1,-1,1];
IsZZGenerating([w1,w2] : T := [2]);
// output : false
