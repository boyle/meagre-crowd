A=zeros(5,5);

A(1,1) =   1.0;
A(3,1) =   2.0;
A(5,1) =   3.0;
A(1,2) =  -4.0;
A(4,2) =   5.0;
A(2,3) =  -6.0;
A(5,3) =  -7.0;
A(1,4) =  -8.0;
A(4,4) =  -9.0;
A(2,5) =  10.0;
A(5,5) =  11.0;

A
save('unsym4.mat', 'A', '-v4');
save('unsym6.mat', 'A', '-v6');
save('unsym7.mat', 'A', '-v7');
save('unsym73.mat', 'A', '-v7.3');

A=sparse(A);
save('unsym4sp.mat', 'A', '-v4');
save('unsym6sp.mat', 'A', '-v6');
save('unsym7sp.mat', 'A', '-v7');
save('unsym73sp.mat', 'A', '-v7.3');
