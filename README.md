#    <div style="text-align: center"> The Network Epidemiologist's DIY Handbook </div>
######  <div style="text-align: right"> - Ashwin Karichannavar, hal9k000, Vishwesha Guttal </div>
---

## A Soft Prelude:

The agent-model used here is based on *Steven Strogatz's and Duncan Watt's* ubiquitous Small World Network architecture <sup> [1] </sup>. Having created such a network with an *edge rewiring probability (p) = 1* and a *K value (number of adjacent nodes on average to any given node) (K) = 18-20*. The nodes are updated as per the classic SEIRD protocol, additionally accounting for whether the node is presently contagious or not as well as the severity of the infection.

To learn more about the updatation rules, as well as the various parameters used please refer to the PDF document linked to in the references below <sup> [2] </sup>.

The Indian Goverment has taken various steps since the start of the outbreak (lockdowns, a limited contact tracing and testing policy, ban on international arrivals etc), which usually have been initiated in different phases (corresponding to different time points). Pre-lockdown measures have been noted in the following document <sup> [3] </sup>. These measures in turn have been modelled in *COVID19_SmWorld_Intervention.py*.

## The Starter Kit:

To run the code above, you'll need *Python 3.X* installed on your workstation, together with the following third-party Python packages:
* NetworkX <sup>[4] </sup>
* NumPy
* Pandas
* Matplotlib

Additionally, applications such as [Cytoscape](https://cytoscape.org/ "Cytoscape") or [Gephi](https://gephi.org/ "Gephi") might be useful for graph visualisation.


## A Prepper Stash For COVID-19:

---

**Note the directory hierarchies and dependencies in the various scripts before altering them. If in doubt, revert to the hierarchy present above.**

----

The tutorial here will pre-eminently focus on the *COVID19_SmWorld_Basic1.py* script present in the above repository. We will be going over the major methods and their functionality one by one. Make sure to pay attention to the comments peppered throughout the scripts. Again, to learn More about the updatation rules and what the various rate parameters stand for, please refer to the PDF document linked to in the references below <sup> [2] </sup>.

* Firstly, we should note that *COVID19_SmWorld_Basic1.py* script is used to simply simulate the natural disease progression in the *absence of any preventive measures (lockdowns, contact tracing, ban on international arrivals etc)*, which we shall henceforth refer to as the *base case*.

* The first method to take stock of in the aforementioned script is `rates(self)`. This is used to initialise all the parameters in abeyance with experimental and clinical evidence <sup> [2] </sup>.

* There are two contropanel method calls. `controlpanel(self)` as well as `time_evolution(self, t)` are atavistic. `controlpanel2(self)` and `gillespie_time_evolution(self)` are the newer method calls wih similar functionality used in their place, having been formulated to make the entire simulation in accordance the Gillespie stochastic protocol (Section 1.4 in [2] ).

* The `node_edge_annotation(self)` method (called from within `controlpanel2(self)`), is to initialize the various states (*"infected", "susceptible", "moderate"* etc) at the start of the simulation, ie *t = 0*. **Note that all these states are stored throughout the simulation as node attributes or labels**. So a "node B", that at t=T is in infected ( I ), contagious (transmitting) ( I<sub>t</sub> ) and quarantined ( I<sub>q</sub> ) states  will have *node labels "I", "T" and "QNS" associated with it* at t=T in the code. These node labels are stored independently of each other and are saved as attributes called *"state", "transtate" and "type"* respectively. The other node attributes *"tstate", "t_trans" and "t_type"* store the time and random numbers used to determine transition probabilites for *"state", "transtate" and "type"* attributes respectively. There is one further node attribute *"source"* which only plays a role in *COVID19_SmWorld_Intervention.py*. This is used to mark whether a particular node is a local citizen ("L") or an international arrival ("F").

* The variable `self.time_pool` (initialised in `__init__(self)`) acts a nodal clock (randoom numbers from an exponential distribution are assigned to each and every node). In each iteration inside `gillespie_time_evolution(self)`, the node corresponding to the minimal nodal clock value ( `self.event_node = self.time_pool.argmin()` ) ( *let us refer to it as "Node A"* ) is updated according to previously established rules. At the end of the updation we have:

      self.clock = self.clock + self.time_pool[self.event_node]
      self.time_pool = self.time_pool - self.time_pool[self.event_node]

 where the current time of the simulation is advanced by the minimum nodal clock value, which is also subtrcted from the other nodal clock values. Another random number drawn from an exponential distribution is then assigned to *"Node A"*. The rate of decay of this exponetial distribution is arrived at via the function `getRate(self)`.

* The `susceptible_updater(self)` and `exposed_updater(self)` methods are used to update *"Node A"* if it happens to be in the susceptible or the exposed pools respectively, as per updation rules mentioned in the reference <sup> [2] </sup>.

* Similarly the `infected_updater(self)` method is used to update *"Node A"* if it happens to be in the infected pool. There are multiple checks (and potentially multiple transitions) that need to be checked for here, which are listed out in Sections 1.2.4 - 1.2.7 in [2].

* The functions `stat_gen(self,t)` and `statistics(self)` are involved in generating various agglomerate data (such as, for example the total number of contagious individuals at the current time step of the simulation), and plotting out this agglomerate data respectively.

* The `output_graph(self, t)` method can optionally be called from within `stat_gen(self, t)` and is used to save the network structure at each time step of the simulation as a *.graphml* file, along with all the node attributes (labels) and can be plotted in **Cytoscape**.

* Moving onwards, the *COVID19_SmWorld_Intervention.py* script borrows various methods from *COVID19_SmWorld_Basic1.py* and is used to model the effects of various interventions imposed in India. Note that the script **does not** currently incorporate the impact of the lockdown on disease spread. The interventions that have been imposed, follow figures laid out in [3].

* In *COVID19_SmWorld_Intervention.py* we have a new function `rates2(self)` that stores various parameters associated with the interventions imposed <sup> [3] </sup>.

* There are three more new methods present in *COVID19_SmWorld_Intervention.py*, which are namely:
     - `int_arrivals(self, t)`, which randomly turns a certain number of susceptible nodes in the population into infected international arrivals (who may be in I<sub>e</sub>, I<sub>a</sub> or I<sub>m</sub> states) with every passing day, until the ban on all international flights goes into effect.

     - `int_contact_tracing(self, t)`, which starting from Phase II mimics a contact tracing and quarantining/hospitilisation policy for symptomatic international arrivals.

     - `loc_contact_tracing(self,t)`, which starting from Phase III mimics a contact tracing and quarantining/hospitilisation policy for symptomatic local cases.

## Currently Known Issues With The Code:

* The `loc_contact_tracing(self,t)` and `int_contact_tracing(self, t)` methods haven't been suitably modified for the application of Gillespie protocol.

* The `beta_net` variable computed in `susceptible_updater(self)` could change for a fixed node with each time step. Hence this rate is dynamic and cannot be used to


## References:

[1]: http://worrydream.com/refs/Watts-CollectiveDynamicsOfSmallWorldNetworks.pdf

[2]: https://drive.google.com/file/d/1d3RbUkIKKK2sEzABsyWcK8gXobyums-w/view?usp=sharing "SEIRD Model Documentation"

[3]: https://docs.google.com/document/d/1Jdwl1-IlXidUt7Zm3lVYgugzzodeVirVLcowETVzYDM/edit?usp=sharing "Documents On Indian interventions"

[4]: https://networkx.github.io/documentation/networkx-2.3/reference/introduction.html "NetworkX Introduction"
