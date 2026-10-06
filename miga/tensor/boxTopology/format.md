## Description of topological parts of a axis aligned box

$\Omega=\times_{i=1}^n [a_i,b_i)$ then to all $v\in\{-1,0,1\}^d$ we associate

$$P_{v}=\times_{i=1}^n T_i,\qquad T_i=\begin{cases}a_i & v_i=-1\\ [a_i,b_i)& v_i=0\\ b_i & v_i=1\end{cases}$$

In the program the vector $v$ must be a column vector and lists of sides are
represented as matrices in which each side is column. 
It is suggested to write explicit sides as row vectors and transpose to improve
readability.