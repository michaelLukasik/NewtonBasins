# NewtonBasins
Quick Little program for generating newton basins according to newtons method. The plotting is in a python script to folow. 

The basic idea is to drive each step via the ratio between the fuinction value and its derivative, hopefully driving each point towards its root (or Ia supopose a saddle point). The point can then be colored according to the root it belongs to, or according to how many iterations it took to "settle" in its final position

A nicer version of this is coming in the future but this is all still preliminary/in the works

![-0 785397_0 785397_-0 785397_0 785397_ _30_3000_magma_1e-15](https://github.com/michaelLukasik/NewtonBasins/assets/138163589/6869a0da-64d1-4022-8c37-bdfe8bda0538)
*Newton Basin using p(z) = $j_0(z)$ [Bessel function of order 0]*


![-3 141590_3 141590_-3 141590_3 141590_ _300_1000_copper_1e-12](https://github.com/michaelLukasik/NewtonBasins/assets/138163589/61ac5259-b96e-408e-b381-8f0ac5bcf235)
*Newton Basin using p(z) = $z^3 + 2z^2 + 3z^3 -1$*

![-3 141590_3 141590_-3 141590_3 141590_ _300_1000_gray_r_1e-12](https://github.com/michaelLukasik/NewtonBasins/assets/138163589/38cac2e1-e0ea-44f2-aa41-d68739733639)
*Newton Basin using p(z) = $sinh(z)$*
