function [ matrix ] = skewSym( vector )
%SKEWSYM.m

x = vector(1);
y = vector(2);
z = vector(3);

matrix = [0 -z y;
    z 0 -x;
    -y x 0];

end
