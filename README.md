# rownania-rozniczkowe

# Rozwiązanie układu równań różniczkowych

Równanie: $y''' = 6y'' − 11y' + 6y$

Warunki początkowe: $y(0) = 1, y'(0) = 0, y''(0) = 0$

Rozwiązanie analityczne: $y(t) = 3e^t − 3e^{2t} + e^{3t}$ 

Krok całkowania wynosi 0.1, toleranca błędu dla adaptacyjnuch metod (RKF oraz RKF-CK) wynosi 1e-6, przedział całkowania jest równy $t \in [0, 1]$.
<img width="1678" height="855" alt="obraz" src="https://github.com/user-attachments/assets/8383cfa1-d5b8-45b5-b494-a2480c02e64e" />

# Równanie (układ Lorenza):

Układ Lorenza ma postać: 

$\frac{dx}{dt}  = \sigma (y-x)$,

$\frac{dy}{dt} = x(\rho−z)−y$,

$\frac{dz}{dt} = xy−\beta z​$,

gdzie $\sigma = 10$, $\rho = 28$, $\beta = 8/3$.

Warunki początkowe: $x(0) = 1, y(0) = 1, z(0) = 1$

Rozwiązanie metodą rk4 z krokiem całkowania $h = 0.1$ na przedziale $t \in [0, 100]$ jest przedstawione na poniższym rysunku.

<img width="1078" height="824" alt="obraz" src="https://github.com/user-attachments/assets/44a2ff4e-be57-437f-9a67-768010c3df4b" />


a tutaj rozwiązanie metodą Rungego-Kutty-Fehlberga z początkowym krokiem całkowania $h = 0.1$ oraz tolerancją błędu równą $1e-6$ na przedziale $t \in [0, 100]$.

<img width="861" height="639" alt="obraz" src="https://github.com/user-attachments/assets/904c3b07-7609-486a-aeb1-f8d5e3db6ed0" />



Wybrane Źródła:

Kiusalaas J., Numerical methods in engineering with Python, Cambridge university press, 2010

Press W.H., Teukolsky S.A., Vetterling W.T., Flannery B.P., Numerical Recipes in Fortran 77: The Art of Scientific Computing, vol 1, 1992.













