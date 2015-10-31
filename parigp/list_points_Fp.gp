{
   \\ 2-component vector [a_4,a_6]: Weierstrass equation Y^2 = X^3 + a_4 X + a_6,
   E = ellinit([0,1], 101);
   \\P = random(E);
   P = [75, 10];
   \\ order = ellorder(E, P);
   print("Point order: ", order);
   \\ for(k=1, order+1,
   for(k=1, 20,
      Q=ellmul(E, P, k);
      print("k=" k ", " Q);
   );
   quit();
}

