Equation of Motion

$$ [M]{\ddot x} + [B]{\dot x} + [K]{x} = F(t) $$

Free Vibration
--------------

$$ [M]{\ddot x} + [K]{x} = 0 $$

$$ x(t) = A sin(\omega t) + B cos(\omega t) $$

$$ \dot x(t) = A \omega cos(\omega t) - B \omega sin(\omega t) $$

$$ \ddot x(t) = -A \omega^2 sin(\omega t) - B \omega^2 cos(\omega t) = -\omega^2 x $$

$$ [M]{\ddot x} + [K]{x} = 0 $$

$$ (-\omega^2 [M] + [K]){x} = 0 $$

For the non-trivial solution, this must be 0:

$$ -\omega^2 [M] + [K] = 0 $$

$$ \omega^2 = [M]^{-1}[K] $$

so:

$$ \omega_{natural} = \omega = \sqrt (\frac{K}{M} ) $$

For:

$$ m \ddot x + b \dot x + k x = F $$
$$ \ddot x + \frac{b}{m} \dot x + \frac{k}{m} x = \frac{F}{m} $$
$$ \ddot x + 2 \zeta \omega_n \dot x + \omega_n^2 x = \frac{F}{m} $$

$$ \frac{k}{m} = \omega_n^2 $$
$$ \zeta = \frac{b}{b_{cr}} = \frac{G}{2} $$
$$ \frac{b}{m} = 2 \zeta \omega_n $$
$$ b = 2 \zeta \sqrt{km} $$

For zeta=1

$$ b = b_{cr} $$
$$ b_{cr} = 2 \sqrt{km} = 2 m \omega_n $$

The steady-state solution is:

$$ u(t) = p/k \frac{sin(\omega t + \phi)} {\sqrt{(1-\omega^2 / \omega_n^2)^2 + (2 \zeta \omega / \omega_n)^2}} $$

$$ \phi = -tan^{-1}(\frac{2 \zeta \omega / \omega_n}{1 - \omega^2 / \omega_n^2}) $$

where for unit mass:

$$ k = m \omega_n^2 = \omega_n^2 $$

Frequency Response
------------------
Let:

$$ {x} = U(\omega) e^{j \omega t}$$

So:

$$ {\dot x} = j \omega   U(\omega) e^{j \omega t} $$

$$ {\ddot x} = -\omega^2 U(\omega) e^{j \omega t} $$

Cancelling  $ e^{j \omega t} $ on both sides:

$$ (-\omega^2 [M] U(\omega) + j \omega [B] U(\omega) + [K] U(\omega)) = {F} $$

Solving for $U(\omega)$:

$$ (-\omega^2 [M] + j \omega [B] + [K]) U(\omega) = {F} $$

$$ U(\omega) = (-\omega^2 [M] + j \omega [B] + [K])^{-1} {F} $$

This is the direct response method.


For mass normalized mode shapes and dropping [B]:

$$ [I] = [\hat M] = [\Phi]^T [M] [\Phi] $$

$$ [\omega_n^2] = [\hat K] = [\Phi]^T [K] [\Phi] $$

$$ -\omega^2 [I] + [\omega_n^2] = [\Phi]^T {F} $$

For Free Vibration:

$$ -\omega^2 [\hat M] + [\hat K] = 0 $$

$$ \omega^2 = [\hat M]^{-1} [\hat K] = ([\Phi]^T [M] [\Phi])^{-1} [\Phi]^T [K] [\Phi] $$

$$ (AB)^{-1} = B^{-1} A^{-1} $$

$$ ([\Phi]^T [M] [\Phi])^{-1} = (A [\Phi])^{-1} = [\Phi]^{-1} ([\Phi]^T [M])^{-1} = [\Phi]^{-1} [M]^{-1} [\Phi] $$

$$ \omega^2 = [\Phi]^T [M]^{-1} [K] [\Phi] = [M]^{-1} [K] $$


Modal Dynamics
--------------
Let:

$$ {x} = U(\omega) e^{j \omega t} = [\Phi] {q} e^{j \omega t} $$

$$ {\dot x} = j \omega [\Phi] {\dot q}   e^{j \omega t} $$

$$ {\ddot x} = -\omega^2 [\Phi] {\ddot q} e^{j \omega t} $$

From:

$$ -\omega^2 [M] U(\omega) + j \omega [B] U(\omega) + [K] U(\omega) = {F} $$

Plug in q:

$$ -\omega^2 [M] [\Phi] {q} + j \omega [B] [\Phi] {q} + [K] [\Phi] {q} = {F} $$

Pre-multiply:

$$ -\omega^2 [\Phi]^T [M] [\Phi] {q} + j \omega [\Phi]^T [B] [\Phi] {q} + [\Phi]^T [K] [\Phi] {q} = [\Phi]^T {F} $$

$$ -\omega^2 [\hat M] {q} + j \omega [\hat B] {q} + [\hat K] {q} = [\Phi]^T {F} $$

$$ (-\omega^2 [\hat M] + j \omega [\hat B] + [\hat K]) {q} = [\Phi]^T {F} $$

Given \omega and assuming B=0, solve:

$$ ([\hat K] - \omega^2 [I]) {q} = {\hat F} $$

$$ (\omega_n^2 - omega^2) {q} = {\hat F} $$

Reintroducing B:

$$ (-\omega^2 [\hat M] + j \omega [\hat B] + [\hat K]) {q} = [\Phi]^T {F} $$
$$ (-[\hat M] + j \frac{1}{\omega} [\hat B] + \frac{1}{\omega^2}[\hat K]) {q} = \frac{1}{\omega^2} [\Phi]^T {F} $$

$$ \\ddot q + [2 \zeta \omega_n] \dot q + [omega_n^2] {q} = {P} $$



Modal Transient
---------------
From the definition of q

$$ [M] [\Phi] \ddot q + [B] [\Phi] \dot q + [K] [Phi] q = {F} $$

$$ [\Phi]^T [M] [\Phi] \ddot q + [\Phi]^T [B] [\Phi] \dot q + [\Phi]^T [K] [Phi] q = [\Phi]^T {F} $$

$$ [\hat M] \ddot q + [\hat B] \dot q + [\hat K] q = [\Phi]^T {F} $$


