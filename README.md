## Presentation

In the present project some tools are implemented in Python to draw and compute classical polynomials for graphs embedded in orientable surfaces. The implementations could 
be useful for mathematicians and scientists from other research areas interested in computing polynomials associated to graphs embedded in surfaces. 
<br>
<br>
The project is divided into thre main modules *embedded_graph_canvas.py*, *ribbon_graph_polynomial.py* and *Penrose_polynomial.py*. The module *embedded_graph_canvas.py* 
allows the user to draw an embedded graph in an orientable surface, and the main goal of modules *ribbon_graph_polynomial.py* and *Penrose_polynomial.py* is to compute 
the ribbon graph polynomial and Penrose polynomial, respectively, of a given embedded graph.
<br>
<br>

### Embedded Graphs and Ribbon Graphs

A surface is a topological manifold of dimension 2. In the case of compact connected surfaces there is a complete classification theorem: 
every compact connected surface $F$ is homeonmorphic to the $2$-sphere, to the connected sum of $g\geq 1$ 
torus (if $F$ is orientable) or to the connected sum of $g\geq 1$ projective planes (in the non-orientable case).
<br>
<br>
A compact connected surface without boundary can be represented as a polygon with an even number of edges which are pairwise identifyed using the quotient topology. 
In the following picture a torus is decomposed into an octagon after *cutting* along four curves having a common point. Each one of the  cutting curves is represented 
by two edges in the octagon with the same color and type of arrow; the arrows indicate the orientation of the identification to recover the surface. In general, given an
compact orientable surface of genus $g$ (homeomorphic to the connected sum of $g$ tori), it can be represented by a $4g$-gon with pairs of edges identifyed.   

<br>
<figure>
  <img src="/images/identify.JPG" width="700" >
  <figcaption>Figure 1. A compact surface represented by a polygon with edges pairwise identifyed. </figcaption>
</figure>
<br>
<br>

On the other hand, a graph is fundamental combinatorial object defined by two sets: a set of vertices (V), and a set of edges (E). In recent years the study of graphs 
and the spaces in which they can be embedded is taking relevance, particularly, the understanding of  embedding of graphs in surfaces is an active research topic.


**Definition** A graph $G(V,E)$ is *cellularly embedded* in a compact connected surface $F$ if $G\subset F$ is such that the edges of $G$ are simple arcs
which intersect at their endpoints in vertices, and such that each component in $F\setminus G$ is homeomorphic to a disc. Two graphs $G_1$  and $G_2$ cellularly
embedded in a surface $F$ are *isomorphic* as embedded graphs is there is a homeomorphism $\varphi: F\rightarrow F$ such that $\varphi|_{G_1}: G_1\rightarrow G_2$ 
is an isomorphism of abstract graphs.      
<br>

A ribbon graph is a topological and combinatorial object that has a duality as a graph and a surface, namely:
<br>

**Definition** A *ribbon graph* $G$ is a surface that has a handle decomposition: a set of $0$-handles or *fat vertices* $V$, and a set of disjoint $1$-handles or *bands*, 
$E$, attached to the $0$-handles in disjoint arcs on their boundaries.    

There is a close relationship between cellularly embedded graphs in compact surfaces and ribbon graphs. Given a cellularly embedded graph $G\subset F$ in the surface $F$, 
taking a regular neighborhood of $G$ in $F$, we get a ribbon graph $\tilde{G}$. On the other hand, if we start with a ribbon graph $\tilde{G}$, after capping off all 
boundary components of $\tilde{G}$ with discs, using quotient topology, we get  compact surface $F$ without boundary, in which the core $G$ of the initial ribbon graph 
$\tilde{G}$ is cellularly embedded. In the figure below, a cellularly embedded graph in a torus is presented, as well as the ribbon graph associated to it.   

<br>
<figure>
 <img src="/images/embed.JPG" width="700" >
  <figcaption>Figure 2. A cellularly embedded graph and associated ribbon graph. </figcaption>
</figure>
<br>

In the module *ribbon_graph_polynomial.py*, the function *embed_graph_canvas(g,v)* can be used to draw  the embedded graph of the figure above. 

```python
aristas = embed_graph_canvas(1,3)
print(aristas)
```
<br>
### Ribbon Graph Polynomial

The ribbon graph polynomial is an invariant of ribbon graphs introduced  by Bollobás and Riordan. This polynomial generalizes the classical Tutte polynomial. 
As in the case of Tutte polynomial, this polynomial has different formulations, includig recursive contraction-deletion relations and formulations via subgraph states. 
In the following definition a formulation of the Ribbon Graph Polynomial is presented in which all possible subgraph states are considered.  

**Definition** The ribbon *graph polynomial*, $R(G;x,y,z,w)\in Z[x,y,z,w]/ \langle w^{2}-w \rangle$, is defined for a ribbon graph $G$ as

$$ R(G; x,y,z,w)=\sum_{\substack{E(A)\subseteq E(G),V(A)=V(G)}} (x-1)^{r(G)-r(A)}y^{n(A)}z^{k(A)-f(A)+n(A)}w^{t(A)}  $$

where the sum is taken overall subgraph states of $G$, given by all subsets of edges of $G$ and preserving the total of vertices of $G$.
The values $e(A), v(A), f(A), k(A)$ are the numbers of edges, vertices, boundary components and connected components of the graph state $A$, respectively, 
while $r(A) = v(A)- k(A)$ and  $n(A) = e(G)- r(G)$ are the rank and nullity of $A$. The value $t(A)$ is equal to $0$ or $1$ if the ribbon graph $A$ is orientable or not,
respectively. 


The Tutte polynomial can be computed as a particular evaluation of the ribbon graph polynomial \cite{books/daglib/0037859}. 
\begin{lem}The Tutte polynomial of the graph $G$ can be obtained as $T(G;x,y)=R(x,y-1,1,1)$
\end{lem}
As a consequence of this relation, both polynomials share topological information concerning topological information of the graphs, such as connectivity and 
the existence of spanning subgraphs, but the ribbon graph polynomial encode additional information concerning the surface nature of ribbon graphs such as orientability
and the number of boundary components. 
\[
4*x + y**3*z**2 + 2*y**2*z**2 + y**2*(x - 1) + 3*y**2 + y*(x - 1)**2 + 5*y*(x - 1) + 9*y + (x - 1)**2 + 1
\}


### Penrose Polynomial

The Penrose polynomial is an invariant of graphs embedded in surfaces which is related to important problems in graph theory, including graph colourings. 
This polynomial was suggested as an invariant for planar graphs in 1971  by Roger Penrose \cite{Penrose1971} and extended for graphs embedded in surfaces by 
Ellis-Monaghan and Moffat in 2013 \cite{ELLISMONAGHAN2013424}. 
![Logo](/images/medial.JPG)

In order to compute the Penrose polynomial of a cellularly embedded graph $G$ in a closed surface $F$, first, construct the medial graph $G_m$ associated to $G$, 
which has a vertex of degree $4$ on each edge of $G$. We connect the vertices of $G_{m}$ with arcs running parallel to the face boundaries of $G$ in $F$ to obtain the 
edges of $G_m$. We show in Figure \ref{} the construction of the medial graph $G_{m}$ for the graph $G$ embedded in a torus. An important property of a medial graph is
that its complementary regions in the surface $F$ can be coloured in a checker-board fashion (see Figure \ref{fig:f4}). As in the case of classical polynomials invariants
in knto theory, every vertex of the medial graph $G_m$ admits three possible smoothing states: crossing state, white smoothing (connecting the two white regions concerning 
the vertex) and black smoothing. A \emph{Penrose state} of $G_m$ is a ribbon graph obtained after smoothing each vertex of $G_m$ as a white smoothing or a crossing state. 
Then, if $G$ has $n$ edges, its medial graph $G_m$ admits $2^{n}$ Penrose states. The Penrose polynomial of graph $G$ is defined by:

     
\begin{df}[Penrose polynomial]\label{df:pen} Given a cellularly embedded graph $G$  in a closed surface $F$, with medial graph $G_{m}$, the \emph{Penrose polynomial} of $G$, 
$P(G; \lambda)\in Z[\lambda]$, is defined by:

\[ P(G; \lambda)=\sum_{s\in \mathcal{P}(G_{m})} (-1)^{cr(s)}\lambda^{c(s)} \]
where the sum is taken over all Penrose states $\mathcal{P}(G_{m})$ of $G_m$, while $cr(s)$ and $c(s)$ are the numbers of crossing vertex states and boundary components of 
Penrose state $s$.
\end{df}

\[
2*z**3 - 6*z**2 + 4*z
\]

\section{Implementations}


We propose a Python implementation to compute two polynomials of ribbon and embedded graphs, namely, the ribbon graph polynomial and the Penrose polynomial. 
The implementations are divided into three main modules: $embedded\_graph\_canvas.py$,
 $ribbon\_graph\_polynomial.py$ and $Penrose\_polynomial.py$. We briefly describe each module.  
 


\begin{itemize}

\item $\textbf{embedded\_graph\_canvas.py}$ \\
This module contains functions to allow the user to draw an embedded graph in interactive manner using the mouse or pointer. The main function is called 
$embed\_graph\_canvas(g,v)$, whose parameters are the genus ($g$) of the compact orientable surface in which the graph will be embedded and the number of vertices 
of the graph($v$). The specific tasks in function $embed\_graph\_canvas(g,v)$ are  distributed in the following set of sequential functions:   

\begin{itemize}

\item[] $draw\_canvas(g,v):$ This function receives as input data the genus of the surface ($g$) and the number of vertices of the embedded graph ($v$) and it returns 
the coordinates of the vertices of the polygon that represents the surface and the coordinates of  the vertices of the embedded graph.
\item[] $draw\_graph()$: In this function the canvas (polygon and 'fat' vertices of the graph) to plot the embedded graph is drawn. The user can draw the embedded 
graph in the canvas using the mouse or pointer. The output of the function are the endpoints of the arcs provided by user.  
\item[] $codify\_arcs(segments, polygon, vertex):$ This auxiliary function uses the array of endpoints of arcs provided by user $segments$ and the arrays $polygon$ 
and $vertex$ containing the vertices of the polygon and the vertices of the graph, respectively, to codify each endpoint of an arc indicating the edge or the polygon 
or the vertex of the embedded graph to which it is closer, as well as the order in which it appears in the corresponding set.     
\item[] $simplify\_arcs(segments1,seg):$ The arcs provided by user in the surface canvas, are merged to form the edges of the embedded graph.  
\end{itemize} 
The output of function  $embed\_graph\_canvas(g,v)$ is an array called $aristas$ of size $(2e,2)$, where $e$ is the number of edges of the embedded graph and the values 
$aristas[i,:]$ and $aristas[i+e,:]$ represent the endpoints of the same arc of the embedded graph. The array $aristas$ corresponding to the embedded graph in Figure \ref{} 
are shown bellow. For instance, in this example, the entrances  $[0,0]$ and  $[1,2]$ represent the two endpoint of the same arc, and indicate that the $0$-th point around 
the $0$-th vertex of the graph is connected to the second point around the vertex with label $1$.    
\[ [[1 ,1],[0, 0],[0, 1],[2, 2],[2, 1],[1, 3],[1, 2],[2, 0],[0, 2],[1, 0]] \]

%def embed_graph_canvas(g,v):
%    polygon, vertex = draw_canvas(g,v)
%    segments = draw_graph()
%    segments1, seg = codify_arcs(segments, polygon, vertex)
%    aristas = simplify_arcs(segments1,seg)
%    return aristas
      
\item $\textbf{ribbon\_graph\_polynomial.py}$\\
In the present module contains functions to compute the Penrose polynomial of and embedded graph in a closed orientable surface. The main function in the module is 
$RGpolynomial(aristas)$. The input of this function, $aristas$, is a $(2e,2)$-size array in the format of the output of the function $embed\_graph\_canvas(g,v)$ from 
the $\textbf{embedded\_graph\_canvas.py}$ module, which encodes the endpoints of the embedded graph. The secondary functions in the module are:
\begin{itemize}
\item[] $components(A,v):$ This function computes the number of components of graph state $A$.
\item[] $faces(A,v):$ The function returns the number of faces or boundary components of graph state $A$.
\item[]  $RGpolynomial\_state(A,v):$  This function computes the term contribution to the Penrose polynomial of the graph state $A$. The function $RGpolynomial(aristas)$ 
computes and adds the contribution of every graph state with the function $RGpolynomial\_state(A,v)$. 
\end{itemize} 
\item $\textbf{Penrose\_polynomial.py}$

\end{itemize}


\section{Aplications}


\begin{figure}
  \centering
    \includegraphics[width=12cm]{f49}%[width=8cm]{b2}
  \caption{Poincaré homology sphere}
  \label{fig:poinc}
\end{figure}



