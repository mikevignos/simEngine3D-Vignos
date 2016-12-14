function yp = unforced1(t,y)
c = 2; m = 1; k = 100;
yp = [y(2);(-((c/m)*y(2))-((k/m)*y(1)))]; 

