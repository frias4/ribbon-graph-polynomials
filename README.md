## Presentation

In the present project some tools are implemented in Python to draw and compute classical polynomials for graphs embedded in orientable surfaces. The implementations could 
be useful for mathematicians and scientists from other research areas interested in computing polynomials associated to graphs embedded in surfaces. 
<br>
<br>
The project is divided into three main modules *embedded_graph_canvas.py*, *ribbon_graph_polynomial.py* and *Penrose_polynomial.py*. The module *embedded_graph_canvas.py* 
allows the user to draw an embedded graph in an orientable surface, and the main goal of modules *ribbon_graph_polynomial.py* and *Penrose_polynomial.py* is to compute 
the ribbon graph polynomial and Penrose polynomial, respectively, of a given embedded graph.
<br>
<br>

## Embedded Graphs and Ribbon Graphs

A surface is a topological manifold of dimension 2. In the case of compact connected surfaces there is a complete classification theorem: 
every compact connected surface $F$ is homeomorphic to the $2$-sphere, to the connected sum of $g\geq 1$ 
torus (if $F$ is orientable) or to the connected sum of $g\geq 1$ projective planes (in the non-orientable case).
<br>
<br>
A compact connected surface without boundary can be represented by a polygon with an even number of edges which are pairwise identified and inducing the quotient topology. 
In Figure 1, a torus is decomposed into an octagon after *cutting* along four curves having a common point. Each one of the  cutting curves is represented 
by two edges in the octagon with the same color and type of arrow; the arrows indicate the orientation of the identification to recover the surface. In general, given a
compact orientable surface of genus $g$ (homeomorphic to the connected sum of $g$ tori), it can be represented by a $4g$-gon with pairs of edges identified.   

<br>
<figure>
  <img src="/images/identify.JPG" width="600" >
  <figcaption>Figure 1. A compact surface represented by a polygon with edges pairwise identified. </figcaption>
</figure>
<br>
<br>

On the other hand, a graph is a fundamental combinatorial object defined by two sets: a set of vertices (V), and a set of edges (E). In recent years, the study of graphs 
and the spaces in which they can be embedded is taking relevance, particularly, the understanding of  embedding of graphs in surfaces is an active research topic.


**Definition.** A graph $G=(V,E)$ is *cellularly embedded* in a compact connected surface $F$ if $G\subset F$ is such that the edges of $G$ are simple arcs
which intersect at their endpoints in vertices, and such that each component in $F\setminus G$ is homeomorphic to a disc. Two graphs $G_1$  and $G_2$ cellularly
embedded in  surfaces $F_1$ and $F_2$ are *isomorphic* as embedded graphs is there is a homeomorphism $\varphi: F_1\rightarrow F_2$ such that $\varphi|_{G_1}: G_1\rightarrow G_2$ 
is an isomorphism of abstract graphs.      
<br>

A ribbon graph is a topological and combinatorial object that has a duality as a graph and a as a surface, namely:
<br>

**Definition.** A *ribbon graph* $G$ is a surface that has a handle decomposition: a set of $0$-handles or *fat vertices* $V$, and a set of disjoint $1$-handles or *bands*, 
$E$, attached to the $0$-handles along disjoint arcs on their boundaries.    

There is a close relationship between cellularly embedded graphs in compact surfaces and ribbon graphs. Given a cellularly embedded graph $G\subset F$ in the surface $F$, 
after taking a regular neighborhood of $G$ in $F$, we get a ribbon graph $\tilde{G}$. On the other hand, if we start with a ribbon graph $\tilde{G}$, after capping off all 
boundary components of $\tilde{G}$ with discs and inducing quotient topology, we get a compact surface $F$ without boundary, in which the core $G$ of the initial ribbon graph 
$\tilde{G}$ is a cellularly embedded graph in $F$. In the Figure 2, a cellularly embedded graph in a torus is presented, as well as the ribbon graph associated to it.   

<br>
<figure>
 <img src="/images/embed.JPG" width="600" >
  <figcaption>Figure 2. A cellularly embedded graph and its associated ribbon graph. </figcaption>
</figure>
<br>
<br>
 
In the module *embedded_graph_canvas.py*, the function *embed_graph_canvas(g,v)* can be used to draw  the embedded graph of Figure 2. The parameter
$g$ is the genus of the containing surface, in this case $g=1$ since the surface is a torus, while $v$ is the number of vertices of the embedded graph, $v=3$ in our
case.

```python
aristas = embed_graph_canvas(1,3)
print(aristas)
```
<br>

A canvas is displayed with a polygon (square, since $g=1$) representing the surface and a set of blue discs representing the vertices of the embedded graph (left-hand
image in Figure 3). The user can draw the edges of the embedded graph using the mouse or a pointer to mark the vertices of each arc, and pressing *enter* when an arc is complete. When the embedded graph is complete, the user should press *enter* button. In the right side image in Figure 3, the embedded graph of Figure 2 was drawn.

<figure>
 <img src="/images/ejemplo.jpg" width="700" >
  <figcaption>Figure 3. The Ribbon Graph Canvas to draw an embedded graph. </figcaption>
</figure>
<br>
<br>

The output of the *embed_graph_canvas(g,v)* function is an embedded graph represented by an array named *aristas* of size $(2e,2)$, where $e$ is the number of edges of the embedded graph and the values 
$aristas[i,:]$ and $aristas[i+e,:]$ represent the endpoints of the same edge of the embedded graph. In our example we get the array:

$$ [[1 ,1],[0, 0],[0, 1],[2, 2],[2, 1],[1, 3],[1, 2],[2, 0],[0, 2],[1, 0]] $$

Note that the entrances  $[0,0]$ and  $[1,2]$ in the array represent the two endpoint of the same edge of the graph, and indicate that the $0$-th point around 
the $0$-th vertex of the graph is connected to the second point around the vertex with label $1$ by an edge of the embedded graph.    
<br>

## Ribbon Graph Polynomial

The ribbon graph polynomial is an invariant of ribbon graphs introduced  by Bollobás and Riordan. This polynomial generalizes the classical Tutte polynomial. 
As in the case of Tutte polynomial, this polynomial has different formulations, includig recursive contraction-deletion relations and formulations via subgraph states
(subgraphs having the same set of vertices as the original graph, but probably less edges). In the following definition a formulation of the ribbon graph polynomial is presented in which all possible subgraph states are considered.  

**Definition.** The ribbon *graph polynomial*, $R(G;x,y,z,w)\in Z[x,y,z,w]/ \langle w^{2}-w \rangle$, is defined for a ribbon graph $G$ as

$$ R(G; x,y,z,w)=\sum_{\substack{E(A)\subseteq E(G),V(A)=V(G)}} (x-1)^{r(G)-r(A)}y^{n(A)}z^{k(A)-f(A)+n(A)}w^{t(A)}  $$

where the sum is taken overall subgraph states of $G$. The values $e(A), v(A), f(A), k(A)$ are the numbers of edges, vertices, boundary components and connected components of the graph state $A$, respectively, while $r(A) = v(A)- k(A)$ and  $n(A) = e(G)- r(G)$ are the rank and nullity of $A$. The value $t(A)$ is equal to $0$ or $1$ if the ribbon graph $A$ is orientable or not,
respectively. 

The ribbon graph polynomial is a polynomial in four variables in the general case, however, when the containing surface is orientable, every subgraph state is 
orientable and then, the exponent of $w$ is allways equal to $0$ and the ribbon graph polynomial has only three variables. 
The Tutte polynomial can be computed as a particular evaluation of the ribbon graph polynomial. 

**Lemma.** The Tutte polynomial of the graph $G$ can be obtained as $T(G;x,y)=R(G;x,y-1,1,1)$

As a consequence of this relation, both polynomials share topological information concerning  the graphs, such as connectivity and 
the existence of spanning subgraphs, but the ribbon graph polynomial encodes additional information concerning the surface nature of ribbon graphs such as orientability
and the number of boundary components. 
<br>

For the ribbon graph in Figure 2, use the function *RGpolynomial(aristas)* in module *ribbon_graph_polynomial.py* to compute its ribbon graph polynomial, which is presented below.

```python
P=RGpolynomial(aristas)
print(P)
```
<br>

```python
 4*x + y**3*z**2 + 2*y**2*z**2 + y**2*(x - 1) + 3*y**2 + y*(x - 1)**2 + 5*y*(x - 1) + 9*y + (x - 1)**2 + 1 
``` 

<br>

## Penrose Polynomial

The Penrose polynomial is an invariant of graphs embedded in surfaces which is related to important problems in graph theory, including graph colourings. 
This polynomial was suggested as an invariant for planar graphs in 1971  by Roger Penrose  and extended for graphs embedded in surfaces by 
Ellis-Monaghan and Moffat in 2013. 

<figure>
  <img src="/images/medial.JPG" width="600" >
  <figcaption>Figure 4. Construction of medial graph and a checkerboard colouring of its complementary regions. </figcaption>
</figure>
<br>
<br>

In order to compute the Penrose polynomial of a cellularly embedded graph $G$ in a closed surface $F$, first, construct the medial graph $G_m$ associated to $G$, 
which has a vertex of degree $4$ on each edge of $G$. We connect the vertices of $G_{m}$ with arcs running parallel to the face boundaries of $G$ in $F$ to obtain the 
edges of $G_m$. We show in Figure 4 the construction of the medial graph $G_{m}$ for the graph $G$ embedded in a torus from Figure 2. An important property of a medial graph is that its complementary regions in the surface $F$ can be coloured in a checker-board fashion (see Figure 4). 

As in the case of classical polynomial invariants in knot theory, every vertex of the medial graph $G_m$ admits three possible smoothing states: crossing state, white smoothing (connecting the two black regions concerning  the vertex) and black smoothing. A *Penrose state* of $G_m$ is a ribbon graph obtained after smoothing each vertex of $G_m$ as a white smoothing or a crossing state. Then, if $G$ has $n$ edges, its medial graph $G_m$ admits $2^{n}$ Penrose states. The Penrose polynomial of graph $G$ is defined by:

     
**Definition.**  Given a cellularly embedded graph $G$  in a closed surface $F$, with medial graph $G_{m}$, the *Penrose polynomial* of $G$, $P(G; \lambda)\in Z[\lambda]$, is defined by:

$$P(G; \lambda)=\sum_{s\in \mathcal{P}(G_{m})} (-1)^{cr(s)}\lambda^{c(s)} $$ <br>
where the sum is taken over all Penrose states $\mathcal{P}(G_{m})$ of $G_m$, while $cr(s)$ and $c(s)$ are the numbers of crossing vertex states and boundary components of 
Penrose state $s$.

The function *Penrose_polynomial(aristas)* in module *Penrose_polynomial.py* allows us to compute the Penrose polynomial of a given embedded graph. For the graph in 
Figure 2, we apply the function *Penrose_polynomial(aristas)*  to compute its Penrose polynomial:

```python
P = Penrose_polynomial(aristas)
print(P)
```
<br>

```python
 2*z**3 - 6*z**2 + 4*z 
``` 
<br>


## Implementations

We briefly describe each module in this project and the main functions that each one contains.  
<br>
 
### embedded_graph_canvas.py  <br>
This module contains functions to allow the user to draw an embedded graph in  an interactive manner using the mouse or a pointer. The main function is called 
*embed_graph_canvas(g,v)*, whose parameters are the genus ($g$) of the compact orientable surface in which the graph will be embedded and the number of vertices 
of the graph ($v$). The specific tasks in function *embed_graph_canvas(g,v)* are  distributed in the following set of sequential functions:   

```python
def embed_graph_canvas(g,v):
    polygon, vertex = draw_canvas(g,v)
    segments = draw_graph()
    segments1, seg = codify_arcs(segments, polygon, vertex)
    aristas = simplify_arcs(segments1,seg)
    return aristas
```

- *draw_canvas(g,v):* This function receives as input data the genus of the surface ($g$) and the number of vertices of the embedded graph ($v$) and it returns 
the coordinates of the vertices of the polygon that represents the surface and the coordinates of  the vertices of the embedded graph.
- *draw_graph():* In this function the canvas (polygon and 'fat' vertices of the graph) to plot the embedded graph is drawn. The user can draw the embedded 
graph in the canvas using the mouse or pointer. The output of the function contains the endpoints of the arcs provided by user.  
- *codify_arcs(segments, polygon, vertex):* This auxiliary function uses the array of endpoints of arcs provided by user $segments$ and the arrays $polygon$ 
and $vertex$ containing the vertices of the polygon and the vertices of the graph, respectively, to encode each endpoint of an arc via the edge of the polygon or the vertex of the embedded graph to which it is closer, as well as the order in which it appears in the corresponding set.     
- *simplify_arcs(segments1,seg):* The arcs provided by user in the surface canvas, are merged to form the edges of the embedded graph.  
 
The output of function  *embed_graph_canvas(g,v)* is an array called $aristas$ of size $(2e,2)$, where $e$ is the number of edges of the embedded graph and the values 
$aristas[i,:]$ and $aristas[i+e,:]$ represent the endpoints of the same edge of the embedded graph. 
<br>
 
### ribbon_graph_polynomial.py
This module contains functions to compute the ribbon graph polynomial of and embedded graph in a closed orientable surface. The main function in the module is 
*RGpolynomial(aristas)*. The input of this function, $aristas$, is a $(2e,2)$-size array in the format of the output of the function *embed_graph_canvas(g,v)* from 
the *embedded_graph_canvas.py* module, which encodes the endpoints the edges of the embedded graph. The secondary functions in the module are:

- *components(A,v):* This function computes the number of components of graph state $A$.
- *faces(A,v):* The function returns the number of faces or boundary components of graph state $A$.
- *RGpolynomial_state(A,v):*  This function computes the term contribution to the ribbon graph polynomial of the graph state $A$. The function *RGpolynomial(aristas)* 
computes and adds the contribution of every graph state with the function *RGpolynomial_state(A,v)*. 
<br>
  
### Penrose_polynomial.py

The *Penrose_polynomial.py* module is devoted to compute the Penrose polynomial of an embedded graph in an orientable surface. The main function in this module is *Penrose_polynomial(aristas)*, which computes the Penrose polynomial of
the embedded graph described by array *aristas*. Note that in the definition of function *Penrose_polynomial()*, the algorithm iterates over all edges of the embedded graph to determine all possible Penrose states of the medial graph. For each one of these Penrose states its contribution to the Penrose polynomial is computed using the auxiliary function *faces(A,I)* that counts the number of boundary components of Penrose state defined by set of indexes *I*.  

```python
def Penrose_polynomial(aristas):
    z = symbols('z')
    P = 0
    a = int(len(aristas[:,0])/2)

    for i in it.product([0, 1], repeat = a):
        I = np.asarray(i)
        P = P + (-1)**len(I[ I == 1]) * z **faces(aristas,np.concatenate((I,I), axis=0))
    return P
```

Note that as we are using the graph state formulation to compute ribbon graph and Penrose polynomials, these algorithms have exponential complexity on the number of edges of the given embeded graph. In the case of the ribbon graph polynomial of embedded graph $G$, all subgraph states are considered, in total $2^{|E(G)|}$ states. In the case of the Penrose polynomial of $G$, we find all Penrose states of $G_m$, which are indexed by all possible subsets of $E(G)$, again we have $2^{|E(G)|}$  Penrose states.  
<br>

## Use Case: Definition of polynomial invariants of lens spaces

We briefly describe the use of the implementations of ribbon graph polynomials and Penrose polynomial to propose polynomial invariants of lens spaces as presented in

[1] J. Frías, J. C. Gómez-Larrañaga, J. L. León-Medina, F. Manjarrez-Gutiérrez, *3-manifold polynomials* (submitted).

A *lens space* is a compact connected 3-manifold without boundary, which can be decomposed as two  disjoint solid tori $H_1$ and $H_2$, which are *glued* together via an orientation reversing homemorphism between their  boundaries $\varphi: \partial H_1 \rightarrow \partial H_2$. Let $\alpha\in\partial H_1$ and $\beta\in \partial H_2$ be simple meridian curves en each solid torus ($\alpha$ and $\beta$ are the boundaries of meridian discs in $H_1$ and $H_2$, respectively). Up to homeomorphism, the lens space is determined by the isotopy class of $\varphi(\alpha) \in H_2$. Every  simple closed curve $\gamma$ in a torus $T$ is parametrized, up to isotopy, by a pair of integers $(p,q)$, where $p$ and $q$ are the numbers of turns in the longitudinal and meridional directions of $\gamma$ around $T$. For instance, the blue curve in Figure 5, has parameters $(3,1)$.  

<figure>
  <img src="/images/lens.png" width="300" >
  <figcaption>Figure 5. Heegaard graph associated to lens space $L(3,1)$. </figcaption>
</figure>
<br>
<br>

In conlusion, a lens space is completely described by a pair of integers $(p,q)$ which parameterize the curve $\varphi(\alpha)$ in $H_2$. In this case, denote the lens space by $L(p,q)$. It is a classical theorem in low-dimensional topology the classification of lens spaces:   Two lens spaces $L(p,q)$ and $L(p',q')$ are homeomorphic if and only if $p'=p$ and $q'\equiv \pm q ^{\pm 1}(mod p)$. 

If a lens space $L(p,q)=H_1\cup H_2$ decomposes as before, under the assumption of minimality and transversality of $\varphi(\alpha)\cap \beta$, we obtain an embedded graph in the torus $\partial H_2$, whose vertices are the points in $\varphi(\alpha)\cap \beta$ and the edges are the complementary arcs of the vertices in $\varphi(\alpha)\cup \beta$. We call this embedded graph, the *Heegaard graph* associated to lens space $L(p,q)$, and it is essentialy unique.

After computing the ribbon graph, Penrose and Tutte polynomials for few Heegaard graphs of lens spaces, it was noticed that homeomorphic lens spaces had same polynomials, as shown in the following list of Penrose polynomial for the first lens spaces: 

[3,1] z**4 + z**3 - 6*z**2 + 4*z <br>
[3,2] z**4 + z**3 - 6*z**2 + 4*z <br>
[4,1] z**6 - 8*z**5 + 39*z**4 - 88*z**3 + 88*z**2 - 32*z <br>
[4,3] z**6 - 8*z**5 + 39*z**4 - 88*z**3 + 88*z**2 - 32*z <br>
[5,1] z**6 + 21*z**5 - 110*z**4 + 200*z**3 - 160*z**2 + 48*z <br>
[5,2] z**5 - 10*z**4 + 45*z**3 - 64*z**2 + 28*z <br>
[5,3] z**5 - 10*z**4 + 45*z**3 - 64*z**2 + 28*z <br>
[6,1] z**8 - 12*z**7 + 123*z**6 - 532*z**5 + 1140*z**4 - 1312*z**3 + 784*z**2 - 192*z <br>
[6,5] z**8 - 12*z**7 + 123*z**6 - 532*z**5 + 1140*z**4 - 1312*z**3 + 784*z**2 - 192*z <br>
[7,1] z**8 + 113*z**7 - 798*z**6 + 2324*z**5 - 3640*z**4 + 3248*z**3 - 1568*z**2 + 320*z <br>

Thanks to this observation, we were able to prove:

*Theorem.* The Penrose, Tutte and ribbon graph polynomials computed for Heegaard graphs are invariants of lens spaces. 

In order to prove that these invariants are complete, it is necessary to recover the parameters $p$ and $q$ from the polynomials. Once again, the analysis of the reduced databases of polynomials for lens spaces led us to conclude and provide a proof that parameter $p$ can be recovered from Tutte and Penrose polynomial, while parameter $q$ can be characterized in some cases.

We also explore the use of polynomials for more  complex 3-manifolds, in particular, we compute the Penrose and ribbon graph polynomials for the classical Heegaard graph of the Poincaré homology sphere:

<figure>
  <img src="/images/f49.jpg" width="600" >
  <figcaption>Figure 6. Poincaré homology sphere and its Heegaard graph. </figcaption>
</figure>
<br>
<br>

```python
z*12 − 24*z**11 + 553*z**10 − 6186*z**9 + 42664*z**8 − 193904*z**7+ 595168*z**6 − 1238528*z**5 + 1718528*z**4 − 1518592*z**3 + 770816*z**2 − 170496*z
```

It is an interesting question, how much topological and geometric information of a 3-manifold can be recovered from the graph polynomials applied to Heegaard graphs.
