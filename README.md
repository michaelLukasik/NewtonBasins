# NewtonBasins
Quick Little program for generating newton basins according to newtons method. The plotting is in a python script to folow. 

The basic idea is to drive each step via the ratio between the fuinction value and its derivative, hopefully driving each point towards its root (or Ia supopose a saddle point). The point can then be colored according to the root it belongs to, or according to how many iterations it took to "settle" in its final position

A nicer version of this is coming in the future but this is all still preliminary/in the works
| Driving Function | Basin after 100 Iterations |
| :-: | :-: |
| $p(z) = J_0(z)$ <br /> <br />  $p'(z) = -J_1(z)$ |  <img src="https://github.com/michaelLukasik/NewtonBasins/assets/138163589/1759249c-6645-4e07-9d09-d55daa98de06" alt="alt text" width="1000" height="700"> 
| $p(z) = z^3 + 2z^2 + 3z^3 -1$ <br /> <br />  $p'(z) = 3z^2 +4z +9$ | <img src="https://github.com/michaelLukasik/NewtonBasins/assets/138163589/9717e975-2795-40cf-8829-af61bf7475c5" alt="alt text" width="1000" height="700">   |
| $p(z) = cos(z)(z^3 -a) e^{-cos(z)-ib}$  <br /> <br />  $p'(z) =e^{-cos(z) - ib} (a-z^3)sin(z)$  <br /> $+$ $cos(z)[(z^3-a)sin(z) +3z^2]$  <br /> <br /><br /> $a = (0, -\pi)$, <br /> $b= (0, 1/\pi)$ | <img src="https://github.com/michaelLukasik/NewtonBasins/assets/138163589/3f33eed0-7465-4a8a-942f-05c4bb7ca46c" alt="alt text" width="1000" height="700">  |



