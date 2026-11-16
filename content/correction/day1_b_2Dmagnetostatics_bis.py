# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.6
#   kernelspec:
#     display_name: ngsolve_latest
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Day I - B) 2D Magnetostatics
#
# - **<font color='green'>[RUN & OBSERVE]</font>** $\rightarrow$ the cell should be run directly without modification
# - **<font color='orange'>[RUN & PLAY]</font>** $\rightarrow$ the cell can be run directly, but some parameters should be changed interactively
# - **<font color='red'>[FILL & RUN]</font>**    $\rightarrow$ the cell should be filled before being run
# - **<font color='magenta'>[FILL & PLAY]</font>** $\rightarrow$ the cell should be filled, and then some parameters should be changed interactively.
#

# %% [markdown]
# ---
# ## Outlines:
#
# <a id='outlines'></a>
#
# - [**1) Linear problem**](#1_linear_problem)
# - [**2) Geometry**](#2_geometry)
# - [**3) Linear problem**](#3_linear_problem)
#     - [a) Physical properties](#3a_physical_properties)
#     - [b) Function space](#3b_function_space)
#     - [c) Integral equation](#3c_integral_equation)
#     - [d) Assemble and solve the system](#3d_assemble_and_solve)
#     - [e) Post-processing: flux density](#3e_post_processing_flux_density)
# - [**4) Nonlinear problem**](#4_nonlinear_problem)
#     - [a) Material properties](#4a_material_properties)
#     - [b) Picard iteration with Aitken acceleration](#4b_picard_iteration)
#     - [c) Simple Newton method](#4c_simple_newton)
#         - [i) Directional derivative](#4ci_directional_derivative)
#         - [ii) Newton loop](#4cii_newton_loop)
#     - [f) Post-processing: flux density](#4f_post_processing_flux_density)
# - [**5) More elaborated Newton**](#5_more_elaborated_newton)
#
# |**<u>Objective</u>** : Solve 2D magnetostatics problems (linear and nonlinear)|
# |---|
# |*Indicative duration* : 2 hours|
#
# <a id='1_linear_problem'></a>
#
# _______
# ## 1) Linear problem
#
# We want to solve the following weak formulation, i.e., find $a\in H^1_0 (\Omega)$ such that
#
# $$ \forall v\in H^1_0(\Omega) , \quad \int_\Omega \text{curl}~ v  \cdot \nu \text{curl}~ a  = \int_\Omega v j $$
#
# - in our 2D problem: $\text{curl}(\cdot) =  \begin{pmatrix} 0 & 1 \\-1 & 0 \end{pmatrix} \text{grad}(\cdot)$
#
# - $\nu$ is the magnetic reluctivity :
#     - $\nu = 1/\mu_0 $ in air and conductors
#     - $\nu = 1/(1000 \mu_0) $ in iron; it might even depend on $B = \text{curl}~ a$
#
# - $j$ is the current density :
#     - $ j = 0 $ in air and iron
#     - $ j= +J $ in positive conductor
#     - $ j = -J$ in negative conductor
#   
# with $J = 10\;A/mm^2$. 

# %% [markdown]
# ________
# ## 2) Geometry
#
# <a id='2_geometry'></a>
#
# We consider a closed iron core with a winding. Without airgap, there is a risk of saturation !
#
# |**<font color='green'>[RUN & OBSERVE]</font>**|
# |---|

# %%
from utils.myGeometries import gaped_inductor
from ngsolve.webgui import Draw
mesh = gaped_inductor(airgap = 0, h = 0.01) # no airgap
Draw(mesh.MaterialCF({"condN":1, "iron":2, "condP":3}), mesh)

print(f"Region names : {mesh.GetMaterials()}")  # region names
print(f"Boundaries names : {mesh.GetBoundaries()}")  # boundaries names

# %% [markdown]
# ________
# ## 3) Linear problem
#
# <a id='3_linear_problem'></a>
#
# ### a) Physical properties
#
# <a id='3a_physical_properties'></a>
#
# Assume first $\nu_{iron} = \frac{1}{1000 \mu0}$, and $J = 10 A/mm^2$.
#
# |**<font color='green'>[RUN & OBSERVE]</font>**|
# |---|

# %%
from numpy import pi
mu0 = 4e-7 * pi   # H/m, void magnetic permeability
mur = 1000        # relative permeability of iron
J = 10e6          # A/m²

# %% [markdown]
# ### b) Function space
#
# <a id='3b_function_space'></a>
#
# |**<font color='red'>[FILL & RUN]</font>**|
# |---|

# %%
from ngsolve import H1
fes = H1(mesh, order = 1, dirichlet = "out")

# %% [markdown]
# ### c) Integral equation
#
# <a id='3c_integral_equation'></a>
#
# We want to solve the following weak formulation, i.e., find $a\in H^1_0 (\Omega)$ such that
#
# $$ \forall v\in H^1_0(\Omega) , \quad \int_\Omega \text{curl}~ v  \cdot \nu \text{curl}~ a  = \int_\Omega v j $$
#
#
# |**<font color='red'>[FILL & RUN]</font>**|
# |---|

# %%
from ngsolve import grad, CoefficientFunction
R = CoefficientFunction( ((0,1),(-1,0)), dims = (2,2))   # -pi/2 rotation matrix
curl = lambda u : R * grad(u)


from ngsolve import BilinearForm, LinearForm, dx
a, v = fes.TnT()
A = BilinearForm(fes)
A += 1/mu0 * curl(a) * curl(v) * dx("air|condP|condN")
A += 1/(mur*mu0) * curl(a) * curl(v) * dx("iron")
l = LinearForm(fes)
l += v * J * dx("condP") 
l += - v * J * dx("condN")

# %% [markdown]
# ### d) Assemble and solve the system
#
# <a id='3d_assemble_and_solve'></a>
#
# |**<font color='green'>[RUN & OBSERVE]</font>**|
# |---|

# %%
K = A.Assemble().mat
f = l.Assemble().vec

from ngsolve import GridFunction
sol = GridFunction(fes)
sol.vec.data = K.Inverse(freedofs = fes.FreeDofs()) * f
Draw(sol)

# %% [markdown]
# ### e) Post-processing : flux density
# We can compute the flux density $B = \text{curl}~a$.
# <a id='3e_post_processing_flux_density'></a>
# |**<font color='red'>[FILL & RUN]</font>**|
# |---|

# %%
B = curl(sol)
Draw(B, mesh, 
     vectors={"grid_size" : 50, "offset" : 0.5 },
     settings = { "Objects" : { "Wireframe" : False, "Surface" : False }})

# %% [markdown]
# The flux density is far too high! The iron saturates at lower induction levels. We should therefore consider the nonlinear material model.


# %% [markdown]
# _____
# ## 4) Nonlinear problem
#
# <a id='4_nonlinear_problem'></a>

# %% [markdown]
# ### a) Material properties
#
# <a id='4a_material_properties'></a>
#
# The magnetic reluctivity in iron depends on the square of the flux density and reads
# $$ \nu_{iron}(b^2) = 100 + 10 \exp(1.8 b^2 ) $$
#
# |**<font color='red'>[FILL & RUN]</font>**|
# |---|

# %%
from ngsolve import exp, pi

def nu_iron(b2):
    """ magnetic reluctivity of iron w.r.t b2"""  
    return 100 + 10 * exp(1.8 * b2)

def dnu_iron_db2(b2):
    """ derivative of the magnetic reluctivity of iron w.r.t b2 """  
    return 18 * exp(1.8 * b2)

#-----------------------------------------------------------------------
# Plot the BH curve
import numpy as np
import matplotlib.pyplot as plt
b = np.linspace(0,2,100)
plt.plot(nu_iron(b**2) * b, b); plt.xlabel("H (A/m)"); plt.ylabel("B (T)")
plt.title("BH curve of iron"); plt.grid(); plt.show()


# %% [markdown]
# We want then to cancel the residual:
#
# $$ R(a; v) = \int_\Omega \nu(|\text{curl}~a|^2) \, \text{curl}~a \cdot \text{curl}~v 
#             - \int_\Omega  j \; v$$
#
# Finding $a \in H^1_0(\Omega)$ such that $R(a; v) = 0 \;\;\forall v \in H^1_0(\Omega)$ is exactly our nonlinear weak formulation.

# %%
def residual(a, 
             v  # test function
            ):  # returns symbolic expression of the residual 
    b2 =   curl(a)**2
    res =  curl(v) * 1/mu0 * curl(a) * dx("air|condP|condN")
    res += curl(v) * nu_iron(b2) * curl(a) * dx("iron")
    res += -v * J * dx("condP") + v * J * dx("condN")
    return res


# %% [markdown]
# ### b) Picard iteration with Aitken acceleration
#
# <a id='4b_picard_iteration'></a>
#
# A simple strategy for nonlinear problems is the **Picard iteration**: freeze the reluctivity at the previous iterate $a_{k}$ and solve a linear problem for $a_{k+1}$:
#
# $$ \forall v \in H^1_0(\Omega), \qquad
#    \int_\Omega \text{curl}~v \cdot \nu(|\text{curl}~a_k|^2) \, \text{curl}~a_{k+1}
#    = \int_\Omega v \, j $$
#
#
# Convergence can be accelerated using **Aitken's $\Delta^2$ method**:
# compute the fixed-point direction $d_k = a_{k+1} - a_k$, then
#
# $$\alpha_k = -\alpha_{k-1}\,\frac{d_{k-1}\!\cdot\!(d_k-d_{k-1})}{\|d_k-d_{k-1}\|^2}.$$
#
# |**<font color='red'>[FILL & RUN]</font>**|
# |---|

# %%
from ngsolve import BilinearForm

def K_operator(a,v,aOld):
    b2 = curl(aOld)**2
    bf =  BilinearForm(curl(v) * 1/mu0 * curl(da) * dx("air|condP|condN"))
    bf += curl(v) * nu_iron(b2) * curl(da) * dx("iron")
    return bf.Assemble().mat

from ngsolve import LinearForm

def F(v):
    lf = LinearForm(v * J * dx("condP") - v * J * dx("condN"))
    return lf.Assemble().vec


# %%
#-----------------------------------------------------------------------
# initialization
from ngsolve import H1, GridFunction
fes = H1(mesh, order = 1, dirichlet = "out") # finite element space
da, v = fes.TnT()
sol = GridFunction(fes)
resList = []
#-----------------------------------------------------------------------
# loop
from ngsolve import InnerProduct, Integrate

flag_Aitken = True
step = 0.01

tol = 1e-8
maxit = 1500
for i in range(maxit):
    # compute and assemble residual and derivatives
    res = LinearForm(residual(sol, v)).Assemble().vec

    # compute norm of the residual
    resList.append(np.linalg.norm(res.FV().NumPy()[fes.FreeDofs()]))

    # compute fixed-point proposal
    a_new = GridFunction(fes)                                             # fixed-point proposal
    a_new.vec.data = K_operator(a,v,sol).Inverse(freedofs = fes.FreeDofs()) * F(v)

    if flag_Aitken:
        d = a_new.vec.FV().NumPy() - sol.vec.FV().NumPy()      # fixed-point direction d_i
        if i > 0:
            dd = d - d_old
            step = -step * np.dot(d_old, dd) / np.dot(dd, dd)
            step = min(max(step, 1e-3), 1.0)

    # update
    sol.vec.data = (1 - step) * sol.vec + step * a_new.vec

    print(f"{i = } | res = {resList[-1]:.5e} | {step = :.2e}")

   
    if flag_Aitken:
        # keep d_i for the Aitken step of the next iteration
        d_old = d.copy()
    else:
        # simple step control
        if i>0 and resList[-1] < resList[-2]:
            step = min(1, 1.2*step)
        else :
            step = step/2

    # simple stop criterion based on the norm of residual
    if resList[-1] < tol:
        break

# %%
plt.semilogy(resList)
plt.grid()
plt.xlabel("Iterations"); plt.ylabel("Norm of the residual");
plt.show()


# %% [markdown]
# ### c) Simple Newton method
#
# <a id='4c_simple_newton'></a>
#
# Now it is time to apply the Newton method. We have to choose the function space and set up the initial guess.
#
# Then, Newton method consists of three steps:
# 1. Find $\delta a\in H(\Omega)$ such that $R'(a_{n}, v ; \delta a) = - R(a_{n}, v) \quad \forall v\in H(\Omega)$
# 2. Update $a_{n+1} = a_{n} + \alpha \delta a $
# 3. Stop criteria : $ |R(a_{n}, \delta a)| \leq \epsilon_1$,  $ |R(a_{n}, v_{i})| \leq \epsilon_1$, $n\geq n_{max}$,...
#
#

# %% [markdown]
# #### i) Directional derivative
#
# <a id='4ci_directional_derivative'></a>
#
# The directional derivative of the residual (needed for Newton, optional for fixed-point):
#
# $$ R'(a; v ; \delta a) = \frac{d}{dt} R(a+t\,\delta a; v)\Big|_{t=0} $$
#
# |**<font color='red'>[FILL & RUN]</font>**|
# |---|

# %%
def d_residual(da,  # trial function
               v,   # test function
               aOld, # previous a
              ):    # returns symbolic expression of the residual directional derivative w.r.t aOld + t*da
    b2 = curl(aOld)**2
    dres =  curl(v) * 1/mu0 * curl(da) * dx("air|condP|condN")
    dres += curl(v) * nu_iron(b2) * curl(da) * dx("iron")
    dres += 2 * curl(v) * dnu_iron_db2(b2) * (curl(aOld) * curl(da)) * curl(aOld) * dx("iron")
    return dres


# %% [markdown]
# #### iii) Newton loop
#
# <a id='4cii_newton_loop'></a>
#
# |**<font color='red'>[FILL & RUN]</font>**|
# |---|

# %%
#-----------------------------------------------------------------------
# initialization
from ngsolve import H1, GridFunction
fes = H1(mesh, order = 1, dirichlet = "out") # finite element space
da, v = fes.TnT()
sol = GridFunction(fes)
resList = []
#-----------------------------------------------------------------------
# loop
from ngsolve import LinearForm, BilinearForm
step = 0.01
tol = 1e-8
maxit = 100
for i in range(maxit):
    # compute and assemble residual and derivatives
    res = LinearForm(residual(sol, v)).Assemble().vec
    dres = BilinearForm(d_residual(da, v, sol)).Assemble().mat

    # compute norm of the residual
    resList.append(np.linalg.norm(res.FV().NumPy()[fes.FreeDofs()]))
    
    # compute descent direction
    d = dres.Inverse(freedofs = fes.FreeDofs()) * res

    # update
    sol.vec.data += - step * d
    print(f"{i = } | res = {resList[-1]:.5e} | {step = :.2e}")
    
    # simple step control
    if i>0 and resList[-1] < resList[-2]:
        step = min(1, 1.2*step)
    else : 
        step = step/2

    # simple stop criterion based on the norm of residual
    if resList[-1] < tol:
        break

# %% [markdown]
# Is it converging? 
# - If not how to make it converge?
# - If yes, at what rate?
#   
# |**<font color='green'>[RUN & OBSERVE]</font>**|
# |---|

# %%
plt.semilogy(resList)
plt.grid()
plt.xlabel("Iterations"); plt.ylabel("Norm of the residual");
plt.show()

# %% [markdown]
# ### f) Post-processing : flux density
# We can compute the flux density $B = \text{curl}~a$.
# <a id='4f_post_processing_flux_density'></a>
# |**<font color='red'>[FILL & RUN]</font>**|
# |---|

# %%
B = curl(sol)
Draw(B, mesh, 
     vectors={"grid_size" : 50, "offset" : 0.5 },
     settings = { "Objects" : { "Wireframe" : False, "Surface" : False }})

# %% [markdown]
# ## 5) More elaborated Newton
# <a id='5_more_elaborated_newton'></a>
# We can improve by : 
# - considering more stopping criteria
# - adding proper linesearch
# - using NGSolve automatic differentiation capabilities
#
# You can analyze the algorithm in the file `utils\mySolver`, and make some comments.
#
# |**<font color='green'>[RUN & OBSERVE]</font>**|
# |---|

# %%
from utils.mySolvers import solveNL

result = solveNL(fes, residual, verbosity = 3)

# %%
plt.semilogy(result["residual"])
plt.grid()
plt.xlabel("Iterations"); plt.ylabel("Norm of the residual");
plt.show()

# %%
Draw(curl(result["solution"]), mesh, 
     vectors={"grid_size" : 50, "offset" : 0.5 },
     settings = { "Objects" : { "Wireframe" : False, "Surface" : False }})

# %% [markdown]
# [*Back to top - outlines*](#outlines)
