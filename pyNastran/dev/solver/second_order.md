Second Order System
-------------------
A problem of the form will be solved:

$$ m \ddot x + b \dot x + k x = F(t) $$

Putting it in standard form (divide by m):

$$ \ddot x + b/m \dot x + k/m x = F(t)/m $$

Let:

 - $\omega_n = \sqrt{k/m}$
 - $\zeta = b / \sqrt{2mk}$

So:

$$ \ddot x + 2 \zeta \omega_n \dot x + \omega_n^2 x = F(t)/m $$

The characteristic equation can be found with the Laplace:

$$ s^2 + 2 \zeta \omega_n s + \omega_n^2 = 0 $$

With roots of:

$$ s = -\zeta \omega_n +/- \omega_n sqrt{z^2 - 1} $$
$$ s = -\zeta \omega_n +/- \omega_d $$

The Laplace Transform assuming 0 initial conditions is:

$$ TF(s) = \frac{ \omega_n^2 }{s^2 + 2 \zeta \omega_n s + \omega_n^2} $$

Transforms
----------
$$ L(0) = 0 $$
$$ L(x) = X(s) $$
$$ L(\dot x) = s X(s) - x(0) $$
$$ L(\ddot x) = s^2 X(s) - s x(0) - \dot x(0) $$

$$ L(sin(\omega t + \phi) ) = \frac{s sin(\phi) + \omega cos(\phi)}{s^2 + \omega^2} $$

$$ L(e^{-\omega t) sin(\omega t + \phi) ) = ? $$

$$ L(e^{-a t) = \frac{1}{s+a} $$
$$ L(e^{a t) = \frac{1}{s-a} $$
$$ L(e^{-a t) sin(\omega t) ) = \frac{\omega}{(s + a)^2 + \omega^2} $$

More
----
The Laplace Transform with initial conditions is:


$$ \ddot x + 2 \zeta \omega_n \dot x + \omega_n^2 x = 0 $$
$$ \ddot x + b \dot x + k x = 0 $$


$$ (s^2 X(s) - s x_0 - \dot x_0) + b (X(s) - x_0) + k X(s) = 0 $$
$$ s^2 X(s) - s x_0 + b X(s) + k X(s) = \dot x_0 + b x_0 $$
$$ (s^2  + b X(s) + k X(s)) - s x_0 = \dot x_0 + b x_0 $$
$$ (s^2  + b + k) X(s) = \dot x_0 + b x_0 + s x_0 $$

$$ X(s) = \frac{\dot x_0 + b x_0 + s x_0}{s^2  + b + k} $$

$$ X_1(s) = \frac{\dot x_0}{s^2 + b + k} $$
$$ X_2(s) = \frac{b x_0   }{s^2 + b + k} $$
$$ X_3(s) = \frac{s x_0   }{s^2 + b + k} $$

$$ Y(s) = \frac{1}{s^2 + 2 \zeta \omega_n s + \omega_n^2} $$

$$ s = -\zeta \omega_n + \omega_d $$
$$ s = -\zeta \omega_n - \omega_d $$

s^2 + 2 \zeta \omega_n s + \omega_n^2 = $$ 
$$ = (s -\zeta \omega_n + \omega_d)(s -\zeta \omega_n - \omega_d) $$
