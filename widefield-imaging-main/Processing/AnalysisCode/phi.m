function x = phi(x)

id = find(x<0);
x(id) = 0;