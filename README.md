# NewtonBasins
Quick Little program for generating newton basins according to newtons method As well as higher order householder methods. Right now, HHM up to n=2 are supported (aka up to Halley's Method), but more to follow. 

The basic idea is to drive each step via the ratio between the fuinction value and its derivative, hopefully driving each point towards its root (or, I supopose, a saddle point). The point can then be colored according to the root it belongs to telling us about the over all structure of the function in the complex plane. Areas of high instability (i.e. where the domains swap quickly in regards to pixel value) represent areas of equal-ish "influence" from the given functions roots. 


A nicer version of this is coming in the future but this is all still preliminary/in the works
| Driving Function  | Basin after 100 Iterations |
| :-: | :-: |
| $p(z) = J_0(z)$ <br /> <br />  $p'(z) = -J_1(z)$ |  <img src="https://github.com/michaelLukasik/NewtonBasins/assets/138163589/1759249c-6645-4e07-9d09-d55daa98de06" alt="alt text" width="1000" height="700"> 
| $p(z) = z^3 + 2z^2 + 3z -1$ <br /> <br />  $p'(z) = 3z^2 +4z + 3$ | <img src="https://github.com/michaelLukasik/NewtonBasins/assets/138163589/9717e975-2795-40cf-8829-af61bf7475c5" alt="alt text" width="1000" height="700">   |
| $p(z) = cos(z)(z^3 -a) e^{-cos(z)-b}$  <br /> <br />  $p'(z) =e^{-cos(z) - b} (a-z^3)sin(z)$  <br /> $+$ $cos(z)[(z^3-a)sin(z) +3z^2]$  <br /> <br /><br /> $a = (0, -\pi)$, <br /> $b= (0, 1/\pi)$ | <img src="https://github.com/michaelLukasik/NewtonBasins/assets/138163589/3f33eed0-7465-4a8a-942f-05c4bb7ca46c" alt="alt text" width="1000" height="700">  |

Varrying the input offset parameter c (a complex value) effects the basin's shape and behvaior depending on the input function. For an example, we look again at the case where the driving function is the zero order Bessel function (sinc function) with a variable c whose value is represented by the red circle traveling along c's path in the complex plane on the top left. The corresponding Newton Basin after 100 iterations is represeted on the right. These types of curves and paths can be evaluated for all values of c, some of which may lead to erratic behavior or improper focusing to the roots of the equation. 


<picture>
<img src= "https://github.com/michaelLukasik/NewtonBasins/blob/master/Examples/BesselFinal.gif" alt="Newton Basin driven by Bessel Function" width="1000" height="700">
</picture>




