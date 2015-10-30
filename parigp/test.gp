{
   \\ 2-component vector [a_4,a_6]: Weierstrass equation Y^2 = X^3 + a_4 X + a_6,
   E = ellinit([0,1]);
   r = ellgroup(E, 43);
   print(r);
   r = ellgroup(E, 101);
   print(r);
   quit();
}
