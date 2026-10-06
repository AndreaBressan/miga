
# Basic data types are 

When not specified all data is assumed to be stored in a $\mathtt{double}$ array.

1. A point in $p\in \mathbb R^n$ is represented by the vector of its coordinates $$\begin{pmatrix}x_1(p) \\ x_2(p)\\\vdots\\ x_n(p)\end{pmatrix}.$$ By extension a set of points $P=\{p_1,\dots,p_m\} \subset \mathbb{R}^n$ is described by a matrix $$\begin{pmatrix}x_1(p_1) & x_1(p_2) & \dots &x_1(p_m)\\\vdots & \vdots&  & \vdots\\ x_n(p_1) & x_n(p_2)& \dots &x_n(p_m) \end{pmatrix}.$$ A $m$-dimensional array is also accepted and it has the same meaning as the matrix obtained by flattening the last $m-1$ dimensions.
2. An axis aligned cuboids in $\mathbb{R}^n$, the cuboid $$[a_1,b_1]\times\dots\times [a_n,b_n]$$ is represented by the matrix array of size $n\times 2$ $$\begin{pmatrix}a_1 & b_1\\\vdots & \vdots \\a_n & b_n \end{pmatrix}$$ i.e. by the pair of opposing vertices having respectively minimum and maximum coordinates' values.
3. The representation of axis aligned Cartesian grids $G$ can be compressed by using a cell array containing $$\{x_1(G),\dots,x_n(G)\}$$ where $x_j(G)=\{t:\exists p\in G:x_j(p)=t \}$. 
4. A partial derivative of a function is represented by the associated multi-index as column vector $$\begin{pmatrix}o_1\\\vdots\\o_n\end{pmatrix}\rightarrow\frac{\partial^{o_1}}{\partial x_1^{o_1}} \dots \frac{\partial^{o_n}}{\partial x_n^{o_n}}.$$ Only integer values of $o_1,\dots,o_n$ are allowed. A list of partial derivatives is represented by a matrix or a multidimensiona array in which each column represents a partial derivative, i.e., the same format as for points, but only integer values are accepted.
