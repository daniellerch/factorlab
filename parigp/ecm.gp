ecm(N, B = 1000!, nb = 100)=
{
  for(a = 1, nb,
    iferr(ellmul(ellinit([a,1]*Mod(1,N)), [0,1]*Mod(1,N), B),
      E, return(gcd(lift(component(E,2)),N)),
      errname(E)=="e_INV" && type(component(E,2)) == "t_INTMOD"))
}

{
   f=ecm(2^101-1);
   print(f);
   quit();
}
