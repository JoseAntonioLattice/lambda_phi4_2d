# 2d $\lambda \phi^4$ theory

This code simulates the 1-component $\lambda \phi^4$ theory.

## Requisites

- A Fortran compiler: gfortran
- Makefile installed
- Use Linux

# Theory
Let $\phi(x)\in\mathbb{R}$ be a field. The Lagrangian density in Euclidian time reads
 $$ \mathcal{L} = {1 \over 2}\partial_{\mu}\phi(x)\partial_{\mu}\phi(x) +  {1 \over 2}m^2\phi(x)^2 +  {1 \over 4}\lambda \phi(x)^4,$$
 where $\lambda \ge 0$, for the potential to be bounded from below. The continuum action reads

$$S[\phi] = \int dx^d \ \mathcal{L}(\phi,\partial_{\mu}\phi).$$

With the lattice regularization
$$
\begin{aligned}
\phi(x) & \to  \phi_x \\
\partial_{\mu}\phi(x) & \to  \frac{\phi_{x+a\hat{\mu}} - \phi_x}{a}\\
\int dx^d & \to  \sum_x a^d
\end{aligned}
$$

The unit vectors are defined as
$$ \hat{1}=  \begin{pmatrix}
0 \\
1 \\
\end{pmatrix}, \ \hat{2}=  \begin{pmatrix}
1 \\
0 \\
\end{pmatrix}
$$

On the lattice the action takes the form
$$ S[\phi] = \sum_x a^d \left\{ \frac{1}{2}\sum_{\mu}\left( \frac{\phi_{x+a\hat{\mu}} - \phi_x}{a}\right)^2 + \frac{1}{2}m^2\phi_x^2 + \frac{1}{4}\lambda \phi_x^4 \right\}.$$

In lattice units we set $a = 1$.

Expanding the action through the lattice
$$
\begin{aligned}
S[\phi] & =   \frac{1}{2}\sum_{\mu=1}^2\left(\phi_{(1,1)+\hat{\mu}} - \phi_{(1,1)}\right)^2 + \frac{1}{2}m^2\phi_{(1,1)}^2 + \frac{1}{4}\lambda \phi_{(1,1)}^4 \\
& +  \frac{1}{2}\sum_{\mu=1}^2\left(\phi_{(1,2)+\hat{\mu}} - \phi_{(1,2)}\right)^2 + \frac{1}{2}m^2\phi_{(1,2)}^2 + \frac{1}{4}\lambda \phi_{(1,2)}^4 \\
& \vdots  \\
&  +  \frac{1}{2}\sum_{\mu=1}^2\left(\phi_{(i,j-1)+\hat{\mu}} - \phi_{(i,j-1)}\right)^2 + \frac{1}{2}m^2\phi_{(i,j-1)}^2 + \frac{1}{4}\lambda \phi_{(i,j-1)}^4 \\
& \vdots \\
&  +  \frac{1}{2}\sum_{\mu=1}^2\left(\phi_{(i-1,j)+\hat{\mu}} - \phi_{(i-1,j)}\right)^2 + \frac{1}{2}m^2\phi_{(i-1,j)}^2 + \frac{1}{4}\lambda \phi_{(i-1,j)}^4 \\
&  +  \frac{1}{2}\sum_{\mu=1}^2\left(\phi_{(i,j)+\hat{\mu}} - \phi_{(i,j)}\right)^2 + \frac{1}{2}m^2\phi_{(i,j)}^2 + \frac{1}{4}\lambda \phi_{(i,j)}^4 \\
& \vdots \\
& = \cdots \\
& + \frac{1}{2}\left[\left( \phi_{(i+1,j-1)} - \phi_{(i,j-1)}\right)^2 + \left( \phi_{(i,j)} - \phi_{(i,j-1)}\right)^2 \right] + \frac{1}{2}m^2\phi_{(i,j-1)}^2 + \frac{1}{4}\lambda \phi_{(i,j-1)}^4 \\
& \vdots \\
& + \frac{1}{2}\left[\left( \phi_{(i,j)} - \phi_{(i-1,j)}\right)^2 + \left( \phi_{(i-1,j+1)} - \phi_{(i-1,j)}\right)^2 \right] + \frac{1}{2}m^2\phi_{(i-1,j)}^2 + \frac{1}{4}\lambda \phi_{(i-1,j)}^4 \\
& + \frac{1}{2}\left[\left( \phi_{(i+1,j)} - \phi_{(i,j)}\right)^2 + \left( \phi_{(i,j+1)} - \phi_{(i,j)}\right)^2 \right] + \frac{1}{2}m^2\phi_{(i,j)}^2 + \frac{1}{4}\lambda \phi_{(i,j)}^4 \\
& \vdots
\end{aligned}
$$

If we update the field $\phi$ at site $x = (i,j)$, that is, $\phi_x \to \phi'_x$, the action becomes
$$
\begin{aligned}
S[\phi'] & = \cdots \\
& + \frac{1}{2}\left[\left( \phi_{(i+1,j-1)} - \phi_{(i,j-1)}\right)^2 + \left( \phi'_{(i,j)} - \phi_{(i,j-1)}\right)^2 \right] + \frac{1}{2}m^2\phi_{(i,j-1)}^2 + \frac{1}{4}\lambda \phi_{(i,j-1)}^4 \\
& \vdots \\
& + \frac{1}{2}\left[\left( \phi'_{(i,j)} - \phi_{(i-1,j)}\right)^2 + \left( \phi_{(i-1,j+1)} - \phi_{(i-1,j)}\right)^2 \right] + \frac{1}{2}m^2\phi_{(i-1,j)}^2 + \frac{1}{4}\lambda \phi_{(i-1,j)}^4 \\
& + \frac{1}{2}\left[\left( \phi_{(i+1,j)} - \phi'_{(i,j)}\right)^2 + \left( \phi_{(i,j+1)} - \phi'_{(i,j)}\right)^2 \right] + \frac{1}{2}m^2\phi_{(i,j)}'^2 + \frac{1}{4}\lambda \phi_{(i,j)}'^4 \\
& \vdots
\end{aligned}
$$

Therefore the change in the action $\Delta S = S[\phi'] - S[\phi]$ is

$$
\begin{aligned}
\Delta S & = \frac{1}{2}\left[  \left( \phi'_{(i,j)} - \phi_{(i,j-1)}\right)^2 - \left( \phi_{(i,j)} - \phi_{(i,j-1)}\right)^2\right] \\
& + \frac{1}{2}\left[ \left( \phi'_{(i,j)} - \phi_{(i-1,j)}\right)^2 - \left( \phi_{(i,j)} - \phi_{(i-1,j)}\right)^2\right] \\
 & + \frac{1}{2}\left[\left( \phi_{(i+1,j)} - \phi'_{(i,j)}\right)^2 - \left( \phi_{(i+1,j)} - \phi_{(i,j)}\right)^2\right] \\
 & + \frac{1}{2}\left[\left( \phi_{(i,j+1)} - \phi'_{(i,j)}\right)^2 - \left( \phi_{(i,j+1)} - \phi_{(i,j)}\right)^2\right] \\
& + \frac{1}{2}m^2\left(\phi_{(i,j)}'^2 - \phi_{(i,j)}^2\right) + \frac{1}{4}\lambda\left( \phi_{(i,j)}'^4 - \phi_{(i,j)}^4\right) \\
& = \frac{1}{2}\left[ \left(\phi_{(i,j)}'^2 - \phi_{(i,j)}^2\right) - 2\phi_{(i,j-1)}\left( \phi'_{(i,j)} - \phi_{(i,j)}\right) \right] \\
& + \frac{1}{2}\left[ \left(\phi_{(i,j)}'^2 - \phi_{(i,j)}^2\right) - 2\phi_{(i-1,j)}\left( \phi'_{(i,j)} - \phi_{(i,j)}\right) \right] \\
& + \frac{1}{2}\left[ \left(\phi_{(i,j)}'^2 - \phi_{(i,j)}^2\right) - 2\phi_{(i+1,j)}\left( \phi'_{(i,j)} - \phi_{(i,j)}\right) \right] \\
& + \frac{1}{2}\left[ \left(\phi_{(i,j)}'^2 - \phi_{(i,j)}^2\right) - 2\phi_{(i,j+1)}\left( \phi'_{(i,j)} - \phi_{(i,j)}\right) \right] \\
& + \frac{1}{2}m^2\left(\phi_{(i,j)}'^2 - \phi_{(i,j)}^2\right) + \frac{1}{4}\lambda\left( \phi_{(i,j)}'^4 - \phi_{(i,j)}^4\right) \\
& = 2  \left(\phi_{(i,j)}'^2 - \phi_{(i,j)}^2\right)\\
& - 2\left( \phi'_{(i,j)} - \phi_{(i,j)}\right) \left(\phi_{(i-1,j)} + \phi_{(i,j-1)} + \phi_{(i+1,j)} + \phi_{(i,j-1)} \right)\\
& + \frac{1}{2}m^2\left(\phi_{(i,j)}'^2 - \phi_{(i,j)}^2\right) + \frac{1}{4}\lambda\left( \phi_{(i,j)}'^4 - \phi_{(i,j)}^4\right) 
\end{aligned}
$$


## Metropolis Algorithm

We generate a change in the field, $\phi_x \to \phi'_x$ and accept it with probability

$$ p = \min(1,\exp\left(-\Delta S\right)).$$