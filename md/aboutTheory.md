#### Theory

The ordinary differential equations (ODE’s) we use in GeneEx are chemical rate equations that model the expression of a gene (node) in a GRN. We then use numerical methods to solve these ODE’s and output the steady state gene expression of each gene in the GRN. Each equation takes on the following general form:

<img src="images/theory_1.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:400px; height:40px" />



As an elementary example, the dynamics of an isolated gene (*B*) is modeled using the differential equation:

<img src="images/theory_2.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:400px; height:50px" />

where *B<sub>t</sub>* is defined as the expression of gene *B* at time *t*, and *g<sub>B</sub>* and *k<sub>B</sub>*, are parameters that represent the gene's basal production and degradation rates, respectively. In this straightforward case, the steady state is the expression level of gene *B* when *dB<sub>t</sub>/dt=0*:

<img src="images/theory_3.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:400px; height:50px" />

If gene *B* is regulated by another gene *A*, then their interaction is modeled using the shifted Hill function, an equation that approximates the chemical kinetics of ligand bonding. Here we provide an heuristic explanation. Given an activator (*A*), the gene product of gene A, which can potentially bind *n* times to a promotor site on the DNA of gene B (*P<sub>BA</sub>*), we can use the principle of conservation of mass to state:


<img src="images/theory_4.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:400px; height:50px" />


where *P<sub>T</sub>* is the total concentration of bound and unbound promotor. We are neglecting intermediate states where fewer than *n* activators are bound. Subsequently, using the theory of mass-action kinetics, we can state that the rate of concentration change for the activator-promotor complex is:

<img src="images/theory_5.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:400px; height:50px" />

where *k<sub>on</sub>* and *k<sub>off</sub>* are the association and dissociation rate constants, respectively. Using the steady-state solution at which  *dnAP<sub>BA</sub>/dt = 0*, as well as the conservation of mass equation mentiond earlier, we can obtain **the fraction of unbound promotor**, also known as the Hill-Langmuir function:

<img src="images/theory_6.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:400px; height:50px" />

where *n* is the Hill coefficient and *T<sub>BA</sub>* is the threshold constant and is defined as:

<img src="images/theory_7.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:400px; height:50px" />


Although we derived the equations by assuming that *n* is the number of binding sites, *n* is commonly defined as the cooperativity of the regulation in the Hill function, and will be unique for every interaction (*n<sub>BA</sub>*, in the example above). 


In the simple case where gene *B* is activated by gene *A*, we can use a variation of the Hill function, known as the shifted Hill function(*H<sub>s</sub>*), to model *the production rate* of gene B:


<img src="images/theory_8.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:700px; height:90px" />


where *λ<sup>+<sup><sub>BA</sub>* is the fold change of the regulation, or in other words, the degree to which the basal production rate is magnified, in the case of an activator, or diminished, in the case of an inhibitor (*λ<sup>-<sup><sub>BA</sub>*). 

Using the shifted Hill function, we can now construct the ODE that models the rate of expression change of gene B at time t (*B<sub>t</sub>*) as a function of activator A at time t (*A<sub>t</sub>*):

<img src="images/theory_9.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:400px; height:50px" />


The various tools in GeneEx implement a slightly modified version of the above rate equation that put the basal production *g<sub>B</sub>* in terms of the **maximal** production, *G<sub>B</sub>*, where:

<img src="images/theory_10.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:400px; height:50px" />

For the more typical case where gene *B* has multiple regulators, the production parameter is multiplied by the shifted Hill function of each regulator. Here we can see the ODE for the rate of expression change of gene B at time t, when regulated by activator A and inhibitor C:

<img src="images/theory_11.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:400px; height:50px" />


To account for stochastic effects due to cell to cell variation and low copy numbers in individual cells, GeneEx provides an option to include a noise term based on a Wiener process (*W<sub>t</sub>*) with variance *σ<sup>2</sup>*:

<img src="images/theory_12.png"
     style="display: block;
  margin-left: auto;
  margin-right: auto; width:400px; height:50px" />



The ODE’s built in GeneEx are solved using fourth order Runge-Kutta method and SDE’s are solved using Euler-Maruyama method.

