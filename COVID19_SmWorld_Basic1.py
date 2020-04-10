# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 12:04:55 2020

@author: Koustav
"""

import networkx as nex
import numpy as np
import matplotlib.pyplot as plt
import os
import random as ran
import math
import time
import copy

class COVID_19_Basic():
    
    def __init__(self):
        
        p=0.1 #Probability of rewiring each edge
        self.k=18 #Number of closest neighbours per node
        self.n=10000 #Number of nodes in the small world graph to begin with
        self.SmWorldGr= nex.watts_strogatz_graph(self.n, self.k, p)
        #Created Small World Graph.
        #self.attributes()
        self.rates() #Assigns various transmission and clinical rates to patients.
        
        self.time= 150 #Stores the number of days for which the simulation is carried out.

        self.time_pool = np.random.exponential(10, self.n)
        self.time_steps = []
        self._time_steps = []
        self.clock = 0
        
        self.sus_size=[]
        self.exp_size=[]
        self.inf_size=[]
        self.trans_size=[]
        self.rec_size=[]
        self.dead_size=[]
        
        self.sev_size=[]
        self.hosp_size=[]
        self.qsev_size=[]
        self.qnonsev_size=[]

        self._sus_size=[]
        self._exp_size=[]
        self._inf_size=[]
        self._trans_size=[]
        self._rec_size=[]
        self._dead_size=[]
        self._sev_size=[]
        self._hosp_size=[]
        self._qsev_size=[]
        self._qnonsev_size=[]


        self.verbose_time = 0.1
        self.str="dedug"
        self.log_file = 'debug_log.txt'
        with open(self.log_file, 'w') as l:
            pass
        
        #print("Hello")
        
        #for n in range(1,10):
            #print(ran.random())
        
        
    def controlpanel(self):
        
        t_SEIR=nex.get_node_attributes(self.SmWorldGr, 't_SEIR') #Time at wich node last updated it's SEIR value.
        t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
        state=nex.get_node_attributes(self.SmWorldGr, 'state')
        tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')
        
        self.node_edge_annotation() #Annotates nodes & edges initially (at t=0)
        
        for t in range(0, self.time):
            #continue
            self.time_evolution(t)
            state=nex.get_node_attributes(self.SmWorldGr, 'state')
            print(set(state.values()))
            #print("Hashim")
        
        self.statistics()
        #Time to crunch the final numbers

    def controlpanel2(self):
        t_SEIR=nex.get_node_attributes(self.SmWorldGr, 't_SEIR') #Time at wich node last updated it's SEIR value.
        t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
        state=nex.get_node_attributes(self.SmWorldGr, 'state')
        tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')

        self.node_edge_annotation()

        self.gillespie_time_evolution()
        self.statistics()
            
    
    def time_evolution(self, t):   #At each time point, responsible for updating the SEIR and related attributes of the nodes.
        
        '''t_SEIR=nex.get_node_attributes(self.SmWorldGr, 't_SEIR') #Time at wich node last updated it's SEIR value.
        t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
        state=nex.get_node_attributes(self.SmWorldGr, 'state')
        tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')'''
        
        
        '''Dealing with updating exposed individuals to infected status'''
        self.exposed_updater(t)
        '''Dealing with moving infected individuals upto recovery & associated clinical and transmission dynamics'''
        self.infected_updater(t)
        '''Updation from susceptible to exposed'''
        self.susceptible_updater(t)
        
        self.stat_gen(t)         #Compiling data for statistics'''

    def gillespie_time_evolution(self):
        counter = 0
        init = 0
        while self.clock <= self.time:
            self.state=nex.get_node_attributes(self.SmWorldGr, 'state')
            event_node = self.time_pool.argmin()
            event = self.state[event_node]
            #print("Event node: {}".format(event_node))
            if event == 'S':
                #print("Event {}".format(event))
                #self.trans = [event_node]
                self.susceptible=[event_node]
                self.susceptible_updater(self.clock)
                #self.stat_gen(self.clock)
            elif event == 'E':
                #print("Event {}".format(event))
                self.exposed = [event_node]
                self.exposed_updater(self.clock)
                #self.stat_gen(self.clock)
            elif event == 'I':
                #print("Event {}".format(event))
                self.trans = [event_node]
                print(self.trans)
                print("Goober")
                print(event_node)
                self.infected_updater(self.clock)
                #self.stat_gen(self.clock)
            else:
                #print("Event {}".format(event))
                pass

            self.stat_gen(self.clock)
            #print(self.state[event_node])
            self.clock = self.clock + self.time_pool[event_node]
            self.time_pool = self.time_pool - self.time_pool[event_node]
            self.time_pool[event_node] = np.random.exponential(10, 1)[0]
            self.time_steps.append(self.clock)
            #print(self.clock)
            #verbose_time = 1
            if self.clock > self.verbose_time*counter:
                #print(self.clock)
                #self.__reset_memory(0)
                print("Est. Simulation Time : {} hours".format((((time.time()-init)/(60*60))*(self.time-self.clock))/self.verbose_time))
                init = time.time()
                log = """\
                    Time : {} days\n\
                    Susceptible : {}\n\
                    Infected : {}\n\
                    Exposed : {}\n\
                    Recovered : {}\n\
                    Dead : {}\n\n\n""".format(
                        round(self.clock,2), 
                        self.sus_size[-1], 
                        self.inf_size[-1],
                        self.exp_size[-1],
                        self.rec_size[-1],
                        self.dead_size[-1]
                    )
                print(log)
                self.__reset_memory(0)
                with open(self.log_file, 'a') as f:
                    f.write(log)
                counter = counter + 1
        self.__reset_memory(1)


    def __reset_memory(self, s):
        """
        gillespie time evolution methods results in lot of iterations, so for large simulations,
        as the size of the arrays storing the variables increase, the simulation time goes on 
        increasing, so at intervals, either store the data to file, or clear it.

        this method stores the state if the sustem when called into different array
        when passed 0
        """
        if s == 0:
            if len(self.sus_size): self._sus_size.append(self.sus_size[-1])
            if len(self.exp_size): self._exp_size.append(self.exp_size[-1])
            if len(self.inf_size): self._inf_size.append(self.inf_size[-1])
            if len(self.trans_size): self._trans_size.append(self.trans_size[-1])
            if len(self.rec_size): self._rec_size.append(self.rec_size[-1])
            if len(self.dead_size): self._dead_size.append(self.dead_size[-1])
            if len(self.sev_size): self._sev_size.append(self.sev_size[-1])
            if len(self.hosp_size): self._hosp_size.append(self.hosp_size[-1])
            if len(self.qsev_size): self._qsev_size.append(self.qsev_size[-1])
            if len(self.qnonsev_size): self._qnonsev_size.append(self.qnonsev_size[-1])
            if len(self.time_steps): self._time_steps.append(self.time_steps[-1])

            self.sus_size=[]
            self.exp_size=[]
            self.inf_size=[]
            self.trans_size=[]
            self.rec_size=[]
            self.dead_size=[]
            self.sev_size=[]
            self.hosp_size=[]
            self.qsev_size=[]
            self.qnonsev_size=[]
            self.time_steps = []
        else:
            self.sus_size=self._sus_size
            self.exp_size=self._exp_size
            self.inf_size=self._inf_size
            self.trans_size=self._trans_size
            self.rec_size=self._rec_size
            self.dead_size=self._dead_size
            self.sev_size=self._sev_size
            self.hosp_size=self._hosp_size
            self.qsev_size=self._qsev_size
            self.qnonsev_size=self._qnonsev_size
            self.time_steps = self._time_steps


        
        
    def exposed_updater(self,t):
        #print("Genda")
        
        t_SEIR=nex.get_node_attributes(self.SmWorldGr, 't_SEIR') #Time at wich node last updated it's SEIR value.
        t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
        #state=nex.get_node_attributes(self.SmWorldGr, 'state')   #Dict of all nodes with their SEIR status
        tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')
        
        '''ASSUMTION: Individual only becomes transmitting (contagious) after onset of symptoms (I)'''
        
        for n in self.exposed:
            ch=t_SEIR[n][1]
            
            if(typo[n]=='A'):
                pass
                    #print(" Asymptomatic Stoad\t:%d" %(t_typo[n][0]))
            
            #print("Toad:\t%d" %(t_SEIR[n][0]))
            #print("Noad:\t %f" %(t_SEIR[n][1]))
            if( ch> math.exp(-self.gamma_e*(t- t_SEIR[n][0]))):
                #Person becomes infected from exposed
                
                self.SmWorldGr.nodes[n]['state']= 'I'
                self.SmWorldGr.nodes[n]['transtate']= 'T'       
                #Person becomes transmitting ( or only after onset of symptoms)
                self.SmWorldGr.nodes[n]['t_SEIR']= (t, ran.random())
                self.SmWorldGr.nodes[n]['t_trans']= (t, ran.random())
                
                self.infektion.append(n)
                self.trans.append(n)
                if n in self.exposed: self.exposed.remove(n)
                
                '''tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
                t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')'''
                
                ch1=ran.random()
                ch2=ran.random()
                self.SmWorldGr.nodes[n]['t_type']=(t,ch2)
                        
                '''Choosing whether the exposed will go on to develop severe, non-severe, or asymptomatic cases'''
                        
                if( ch1<=self.p_asymp):        #Individual is classified as aymptomatic/mild
                    self.SmWorldGr.nodes[n]['type']= 'A'
                    self.asymp.append(n)
                
                #Chronicles absolute time (day) at which node attained it's current clinical type
                elif (ch1>self.p_asymp and ch1<= (self.p_asymp +self.p_sev)): #Individual will go on to develop a severe case
                    self.SmWorldGr.nodes[n]['type']= 'S'
                    self.sev.append(n)
                else:
                    self.SmWorldGr.nodes[n]['type']= 'NS'
                            
                            
                self.nonsev.append(n)
                            
                ch2=ran.random()
                if( ch2 <= self.P_crit):
                    #NS case will go on to become severe later.
                    self.nonsevcrt.append(n)
                    self.SmWorldGr.nodes[n]['t_type']=(t,ch2, ran.random())
                            
                    '''Third value in tuple is assigned only to NS cases that will turn to S and used to determine when these individuals 
                    will end up in a 'S' case'''
                
                '''if(typo[n]=='A'):
                    print(" Asymptomatic Stoad\t:%d" %(t_typo[n][0]))
                    print(" Asymptomatic Time Transition\t:%d\t%f" %(t_trans[n]))
                    print(" Asymptomatic Transition\t:%s %s" %(tstate[n], n))'''
                
    def susceptible_updater(self, t):
        
        t_SEIR=nex.get_node_attributes(self.SmWorldGr, 't_SEIR') #Time at wich node last updated it's SEIR value.
        #t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
        #state=nex.get_node_attributes(self.SmWorldGr, 'state')   #Dict of all nodes with their SEIR status
        tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        #t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')
        
        '''Checking for transmission of disease from contagious infected persons to susceptible persons'''
        
        new_exp=[]
        #To keep  track of newly minted "E" nodes in each call of this function
        
        
        for n in self.susceptible:
            #print('Hobo')
            beta_net=0;
            count_q= 0 #Stores number of quarantined cases (QS or QNS)
            count_h=0   #Stores number of hospitalsied cases.
            count_a=0       #Stores number of asymptomatic cases.
            count_ms=0      #Stores number of non-quarantined/non-hospitalised severe or moderate cases.
            
            ch=t_SEIR[n][1]
            
            for r in self.SmWorldGr.neighbors(n):
              if(tstate[r]=='T'):
                #Neighbour needs to be transmitting first and foremost.  
                if (typo[r]== 'S' or typo[r]=='NS'):
                    count_ms+=1
                elif(typo[r]== 'QS' or typo[r]=='QNS'):
                    count_q+=1
                elif(typo[r]== 'H'):
                    count_h+=1
                elif(typo[r]== 'A'):
                    count_a+=1
            
            beta_net= count_a*(self.beta_a) + count_ms*(2*self.beta_a) + count_h*(0.5*self.beta_a) + count_q*(0.9*self.beta_a)
            
            if(beta_net==0):
                continue
            
            if( ch> math.exp(beta_net*(t- t_SEIR[n][0]))):
                self.SmWorldGr.nodes[n]['state']= 'E'
                self.SmWorldGr.nodes[n]['t_SEIR']= (t, ran.random())
                self.exposed.append(n)
                self.susceptible.remove(n)
                
            
            
            
        
                        
            
        
                
            
    def infected_updater(self, t):
        #Updating clinical dynamics, transmission cases and recovery/death toll at each step.
        
        t_SEIR=nex.get_node_attributes(self.SmWorldGr, 't_SEIR') #Time at wich node last updated it's SEIR value.
        t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
        state=nex.get_node_attributes(self.SmWorldGr, 'state')   #Dict of all nodes with their SEIR status
        tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')
        
        '''Updating potential T cases to NT based on prob check'''
        #print("Staid")
        print(self.trans)
        _n = self.trans[-1]
        print("Grey")
        for n in self.trans:
            if( state[n] != 'I'):
                #If node is not infectious then it can no longer be contagious.
                self.trans.remove(n)
            else:
                ch=t_trans[n][1]
                #print("Andim %1.4f" %(ch))
                if( ch> math.exp(-self.gamma_t*(t- t_trans[n][0]))):
                    #Person becomes non-transmitting and lists are updated accordingly.
                    
                    self.nontrans.append(n)
                    self.SmWorldGr.nodes[n]['transtate']= 'NT'
                    self.SmWorldGr.nodes[n]['t_trans']= (t, 0)
                    #Chronicles absolute time (day) at which node attained it's current transmission status
                    #print("Lensmart:\t Inf %d\t Cont %d" %(len(self.infektion), len(self.trans)))
                    self.trans.remove(n)
                    #print("Walmart:\t Inf %d\t Cont %d" %(len(self.infektion), len(self.trans)))
                    tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
                    t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
                
                    if(typo[n]=='A'):
                        pass
                        #print(" Asymptomatic Stoad\t:%d" %(t_typo[n][0]))
                        #print(" Asymptomatic Time Transition\t:%d\t%f" %(t_trans[n]))
                        #print(" Asymptomatic Transition\t:%s %s" %(tstate[n], n))
                    '''elif(typo[n]=='NS'):
                        print(" Non Sev Stoad\t:%d" %(t_typo[n][0]))
                        print(" Non Sev Time Transition\t:%d\t%f" %(t_trans[n]))
                        print(" Non Sev Transition\t:%s %s" %(tstate[n], n))
                    elif(typo[n]=='S'):
                        print(" Severe Stoad\t:%d" %(t_typo[n][0]))
                        print(" Severe Time Transition\t:%d\t%f" %(t_trans[n]))
                        print(" Severe Transition\t:%s %s" %(tstate[n], n))'''
                        
                        
                        
        '''Checking if NS persons slated to turn critical get moved up to QS or S, dependinng on circumstances'''
        #for n in self.nonsevcrt:
        _k = [_n] if _n in self.nonsevcrt else []
        for n in _k:
            ch=t_typo[n][2]
            if (ch> math.exp(-self.gamma_h*(t-t_SEIR[n][0]))):
                # Takes a median of 5.0 days for a non-severe case to turn severe.
                
                if( typo[n]== 'NS'):
                    # Person hasn't been diagnosed yet.
                    self.SmWorldGr.nodes[n]['type']='S'
                    #print("Rashid")
                    
                elif (typo[n]==  'QNS'):
                    #Person was diagnosed and was in self-isolation when it turned  severe.
                    self.SmWorldGr.nodes[n]['type']='QS'
                    self.qnonsev.remove(n)      # Removed from list of quarantined non-severe persons.
                    self.qsev.append(n)         #Added to list of quarantined severe persons.
                    #print("Ibn Batuta")
                
                self.SmWorldGr.nodes[n]['t_type']= (t, ran.random())
                self.sev.append(n)      #Person added to list of severe cases.
                self.nonsev.remove(n)
                self.nonsevcrt.remove(n)
                
        
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')             
            
        '''Checking if a particular infectious has gotten back positive diagnostic for COVID-19
        & is hospitalised/quarantined as the case may be'''
        #for n in self.infektion:
        _k = [_n] if _n in self.infektion else []
        for n in _k:
            if (typo[n]== "S" or typo[n]== "NS"):
                ch=t_SEIR[n][1]
                if( ch> math.exp(-self.gamma_q*(t- t_SEIR[n][0]))):
                    #Takes a median of five days from the onset of the symptoms for a person to be diagnosed.
                    
                    if(typo[n]== "S"):
                        #Person has a severe case
                        if (len(self.hosp)< self.hospbeds):
                            '''NEWLY DIAGNOSED SEVERE CASES ARE OFFERED HOSPITAL BEDS, IF ANY ARE AVALIABLE'''
                            self.SmWorldGr.nodes[n]['type']= 'H'
                            self.hosp.append(n)
                        else:
                            #Due to lack of hospital beds, person enters self-quarantine
                            self.SmWorldGr.nodes[n]['type']= 'QS'
                            self.qsev.append(n)
                            
                        self.SmWorldGr.nodes[n]['t_type']= (t, ran.random())
                        #Time at which patient entered hospital (is diagnosed). Second co-ord associated with date of Recovery/Death
                        
                    elif (typo[n]== "NS"):
                        # Less severe cases to self-quarantine.
                        self.SmWorldGr.nodes[n]['type']= 'QNS'
                        self.qnonsev.append(n)
                        self.SmWorldGr.nodes[n]['t_type']= (t, ran.random())
                        
                        if (len(t_typo[n])==3):
                            # If these NS individuals are slated to go rogue.
                            ch2=t_typo[n][2]
                            self.SmWorldGr.nodes[n]['t_type']= (t, ran.random(), ch2)
                        
                        #ch1=ran.random() ; ch2=ran.random()
                        
                            
                        
        
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type') 

        
        t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
        tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')  
        
        
        '''Checking if non-severe quarantined (QNS) persons have recovered.
        '''
        
        
        #for n in self.qnonsev:
        _k = [_n] if _n in self.qnonsev else []
        for n in _k:
            if( state[n] != 'I' ):
                #If person is not infectious yet, they cannot recover.
                #print("Roger Wabbit")
                continue
            #print("Baloo")
            #print(typo[n])
            ch=t_typo[n][1]
            if( ch> math.exp(-self.gamma_a*(t- t_SEIR[n][0]))):
                #Takes a mean of fourteen days from the onset of the symptoms for a person to recover.
                if(len(t_typo[n])==3):
                    # If these NS individuals are slated to go rogue. but have somehow reached the end of their convalescence.
                    #print("Deitrus")
                    continue

                #Recovered patient can no longer be transmitting/non-transmitting
                if(tstate[n]== 'T'):
                    self.SmWorldGr.nodes[n]['t_trans']= (t,0)
                    if n in self.trans: self.trans.remove(n)
                    #Removing from transmitter list
                
                self.SmWorldGr.nodes[n]['transtate']= ''
                
                self.SmWorldGr.nodes[n]['type']= ''
                self.SmWorldGr.nodes[n]['state']='R'
                self.SmWorldGr.nodes[n]['t_SEIR']=(t,0)
                self.SmWorldGr.nodes[n]['t_type']= (t, 0)
                #Time at which patient recovered.
                
                self.recovered.append(n)
                self.qnonsev.remove(n)
                self.nonsev.remove(n)
                self.infektion.remove(n)
        
                
        
        
        '''Now accounting for severe cases'''
        
        t_SEIR=nex.get_node_attributes(self.SmWorldGr, 't_SEIR') #Time at wich node last updated it's SEIR value.
        t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
        state=nex.get_node_attributes(self.SmWorldGr, 'state')   #Dict of all nodes with their SEIR status
        tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')
           
        #for n in self.sev:
        _k = [_n] if _n in self.sev else []
        for n in _k:
            if( state[n] != 'I' ):
                '''If person is not infectious yet, they cannot recover. This can specifically occur in the context of "S" patients who have 
                been detected by contact testing while in their Exposed state and packed off to "QS"'''
                #print("Woger Dabbit")
                continue
            
            if typo[n]=='H':
                #Hospitalised patient.
                ch=t_typo[n][1]
                if( ch> math.exp(-self.gamma_s*(t- t_typo[n][0]))):
                    #Takes a median of 28.6 days after hospitalisation to recover.
                    #Patient is at the end of hospistalisation period.
                    '''15% of those who go on to develop a severe case (and are hospitalised) die'''
                    ch1=ran.random()
                    if( ch1 <= self.p_critical*self.death_hos_crt):
                        self.dead.append(n)
                        self.SmWorldGr.nodes[n]['state']='D'
                        
                    else:
                        # Patient recovers
                        self.recovered.append(n)
                        self.SmWorldGr.nodes[n]['state']='R'
                        
                    #Recovered patient can no longer be transmitting/non-transmitting
                    if(tstate[n]== 'T'):
                        self.SmWorldGr.nodes[n]['t_trans']= (t,0)
                        if n in self.trans: self.trans.remove(n)
                        #Removing from transmitter list
                
                    self.SmWorldGr.nodes[n]['transtate']= ''
                
                    self.SmWorldGr.nodes[n]['type']= ''    
                        
                    self.SmWorldGr.nodes[n]['t_SEIR']=(t,0)
                    self.SmWorldGr.nodes[n]['t_type']= (t, 0)
                    self.hosp.remove(n)         #No longer hospitalised
                    self.sev.remove(n)          # No longer a severe case
                    self.infektion.remove(n)    # No longer infected.
                else:
                    continue
            if typo[n]=='QS':
                #Quarantined severe patient.
                if (len(self.hosp)< self.hospbeds):
                    '''There are empty hospital beds'''
                    self.SmWorldGr.nodes[n]['type']='H'
                    #Patient is promptly hospitalised.
                    self.qsev.remove(n)
                    self.hosp.append(n)
                    continue
                
                ch=t_typo[n][1]
                if( ch> math.exp(-self.gamma_s*(t- t_typo[n][0]))):
                    #Takes a median of 28.6 days after hospitalisation to recover.
                    #Patient is at the end of hospistalisation period.
                    '''30% of those who go on to develop a severe case (and are NOT hospitalised) die'''
                    ch1=ran.random()
                    if( ch1 <= self.p_critical*self.death_qs_crt):
                        self.dead.append(n)
                        self.SmWorldGr.nodes[n]['state']='D'
                        
                    else:
                        # Patient recovers
                        self.recovered.append(n)
                        self.SmWorldGr.nodes[n]['state']='R'
                    
                    #Recovered patient can no longer be transmitting/non-transmitting
                    if(tstate[n]== 'T'):
                        self.SmWorldGr.nodes[n]['t_trans']= (t,0)
                        if n in self.trans: self.trans.remove(n)
                        #Removing from transmitter list
                
                    self.SmWorldGr.nodes[n]['transtate']= ''
                
                    self.SmWorldGr.nodes[n]['type']= ''    
                        
                    self.SmWorldGr.nodes[n]['t_SEIR']=(t,0)
                    self.SmWorldGr.nodes[n]['t_type']= (t, 0)
                    self.qsev.remove(n)         #No longer in quarantine
                    self.sev.remove(n)          # No longer a severe case
                    self.infektion.remove(n)    # No longer infected.
                    
                    
        
        
                    
        ''' Asymptomatic cases to be dealt with finally. N.B. Again we don't consider to be capable of transmission.'''
        
        t_SEIR=nex.get_node_attributes(self.SmWorldGr, 't_SEIR') #Time at wich node last updated it's SEIR value.
        t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
        state=nex.get_node_attributes(self.SmWorldGr, 'state')   #Dict of all nodes with their SEIR status
        tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')

        #for n in self.asymp:
        _k = [_n] if _n in self.asymp else []
        for n in _k:
            if( state[n] != 'I' ):
                #If person is not infectious yet, they cannot recover.
                continue
            
            #print("Asympt t_trans: %d\t%f" %(t_trans[n]))
            ch=t_typo[n][1]
            if( ch> math.exp(-self.gamma_a*(t- t_SEIR[n][0]))):
                #Takes a mean of fourteen days from the onset of the symptoms for an asymptomatic to recover. (ASSUMPTION)
                
                #print("Asymptototic: %s" %(n))
                #print("Asympt t_type: %d\t%f" %(t_typo[n]))
                #print("Asympt typo: %s" %(typo[n]))
                #print("Asympt t_SEIR: %d\t%f" %(t_SEIR[n]))
                #print("Asympt SEIR State: %s" %(state[n]))
                #print("Asympt t_trans: %d\t%f" %(t_trans[n]))
                #print("Asympt transtate: %d\t%f" %(t_trans[n]))
                
                
                #Recovered patient can no longer be transmitting/non-transmitting
                if(tstate[n]== 'T'):
                    self.SmWorldGr.nodes[n]['t_trans']= (t,0)
                    if n in self.trans: self.trans.remove(n)
                    #Removing from transmitter list
                
                self.SmWorldGr.nodes[n]['transtate']= ''
                
                self.SmWorldGr.nodes[n]['type']= ''
                self.SmWorldGr.nodes[n]['state']='R'
                self.SmWorldGr.nodes[n]['t_SEIR']=(t,0)
                self.SmWorldGr.nodes[n]['t_type']= (t, 0)
                #Time at which patient recovered.
                
                self.recovered.append(n)
                self.asymp.remove(n)
                self.infektion.remove(n)
                
                
    def stat_gen(self, t):
        
        self.sus_size.append(len(self.susceptible))
        self.exp_size.append(len(self.exposed))
        self.inf_size.append(len(self.infektion))
        self.trans_size.append(len(self.trans))
        self.rec_size.append(len(self.recovered))
        self.dead_size.append(len(self.dead))
        
        self.sev_size.append(len(self.sev))
        self.hosp_size.append(len(self.hosp))
        self.qsev_size.append(len(self.qsev))
        self.qnonsev_size.append(len(self.qnonsev))

        """
        print("Number of susceptible persons:\t%d" %(len(self.susceptible)))
        print("Number of exposed persons:\t%d" %(len(self.exposed)))
        print("Number of infected persons:\t%d" %(len(self.infektion)))
        print("Number of infected (transmitting) persons:\t%d" %(len(self.trans)))
        print("Number of dead persons:\t%d" %(len(self.dead)))
        
        print("Number of severe cases:\t%d" %(len(self.sev)))
        print("Number of severe quarantined cases:\t%d" %(len(self.qsev)))
        """
        
        #self.ouputgraph(t)
        
        
    def ouputgraph(self, t):
        
        self.SmWorldGrCopy= copy.deepcopy(self.SmWorldGr)
        
        t_SEIR=nex.get_node_attributes(self.SmWorldGrCopy, 't_SEIR')
        t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')
        
        print(t_SEIR)
        
        for (key, val) in t_SEIR.items():
            self.SmWorldGrCopy.nodes[key]['t_SEIR']= val[0]
        for (key, val) in t_trans.items():
            self.SmWorldGrCopy.nodes[key]['t_trans']= val[0]
        for (key, val) in t_typo.items():
            self.SmWorldGrCopy.nodes[key]['t_type']= val[0]
            
        
        
        if(os.path.isdir("Results")==False):
            os.mkdir("Results")
        os.chdir("Results")
        #Changing directory

        if(os.path.isdir("Inter")==False):
            os.mkdir("Inter")
        os.chdir("Inter")
        #Changing directory

        if(os.path.isdir("%s" %(self.str))==False):
            os.mkdir("%s" %(self.str))
        os.chdir("%s" %(self.str))
        
        #Log making
        
        if(os.path.isdir("Graph Data")==False):
            os.mkdir("Graph Data")
        os.chdir("Graph Data")
        
        t_SEIR=nex.get_node_attributes(self.SmWorldGrCopy, 't_SEIR')
        print("t_SEIR:\n\n")
        print(t_SEIR)
        
        t_typo=nex.get_node_attributes(self.SmWorldGrCopy, 't_type')
        print("\nt_typo:\n")
        print(t_typo)
        
        t_trans=nex.get_node_attributes(self.SmWorldGrCopy, 't_trans')
        print("\nt_trans:\n")
        print(t_trans)
        
        state=nex.get_node_attributes(self.SmWorldGrCopy, 'state')
        print("\nState:\n")
        print(state)
        
        typo=nex.get_node_attributes(self.SmWorldGrCopy, 'type')
        print("\nType:\n")
        print(typo)
        
        tstate=nex.get_node_attributes(self.SmWorldGrCopy, 'transtate')
        print("\nTransmission State:\n")
        print(tstate)
        
        print('\nEverything:\n')
        print(self.SmWorldGrCopy.nodes(data=True))
        print('\nEverything:\n')
        #print(self.SmWorldGr.nodes(data="True"))
        
        print("Time:\t%d" %(t))
        
        nex.write_graphml(self.SmWorldGrCopy, '%d_Annotated.graphml' %(t))
        
        os.chdir("../../../../")
        
        
        
    def statistics(self):
        
        sus_size=np.array(self.sus_size)
        exp_size=np.array(self.exp_size)
        inf_size=np.array(self.inf_size)
        trans_size=np.array(self.trans_size)
        rec_size=np.array(self.rec_size)
        dead_size=np.array(self.dead_size)
        
        sev_size=np.array(self.sev_size)
        hosp_size=np.array(self.hosp_size)
        qsev_size=np.array(self.qsev_size)
        qnonsev_size=np.array(self.qnonsev_size)
        
        if(os.path.isdir("Results")==False):
            os.mkdir("Results")
        os.chdir("Results")
        #Changing directory
        
        if(os.path.isdir("Gillespie")==False):
            os.mkdir("Gillespie")
        os.chdir("Gillespie")
        #Changing directory'''

        if(os.path.isdir("%s" %(self.str))==False):
            os.mkdir("%s" %(self.str))
        os.chdir("%s" %(self.str))
        
        #Log making
        
        if(os.path.isdir("Logs")==False):
            os.mkdir("Logs")
        os.chdir("Logs")
        
        #tseries= np.array(list(range(0,self.time)))
        self.data= np.zeros([len(self.time_steps), 11])
        
        self.data[:,0] = self.time_steps
        self.data[:,1] = sus_size
        self.data[:,2] = exp_size
        self.data[:,3] = inf_size
        self.data[:,4] = trans_size
        self.data[:,5] = rec_size
        self.data[:,6] = dead_size
        
        self.data[:,7] = sev_size
        self.data[:,8] = hosp_size
        self.data[:,9] = qsev_size
        self.data[:,10] = qnonsev_size
        
        '''Exporting data as CSV'''
        
        hd_txt= "Time, 1: Suscep Indv, 2: No. Exp Indv, 3: No. Inf Indv, 4: No. Cont Indv, 5: No. Rec Indv, 6: No. Dead Indv, 7: Sev Cases, 8: Hosp Cases, 9:Q Sev Cases, 10: QNS Cases " 
        np.savetxt('Pandemic.csv', self.data, delimiter=",", header=hd_txt,comments="#")
                   
        f=open('log_unabridged.txt','w')
        f.write(" K:\t%d, N:\t%d, Initial Caseload:\t %d, R0:\t%3.2f, p:\t %5.3f" %(self.k, self.n, self.int_caseload, self.R0, self.p))
        
                   
        os.chdir("../")
        
        #Plotting and saving data here.
        
        labels=["Time", "# Susceptible Individuals", "# Exposed Individuals", "# Infected Individuals (Total)", 
                "# Infected Individuals (Contagious)", "# Recovered Individuals",  "# Dead Individuals",
               "# Severe Cases", "# Hosp Individuals", "# Quarantined Severe Cases",
               "# Quarantined Non-Severe Cases"]
        
        for x in range(1,7):
            #Plotting basic SEIRD data.
            
            plt.plot(self.data[:,0], self.data[0:,x],  marker='.', markerfacecolor='none', 
                     label="%s" %(labels[x]))
        plt.xlabel("Time")
        plt.ylabel("Current Individual Count")
        
        plt.legend()
        plt.savefig("SEIRD.png", dpi=300)
        plt.show()
        plt.clf()
        plt.close()
    
        
        for x in range(6,10):
            
            plt.plot(self.data[:,0], self.data[0:,x],  marker='.', markerfacecolor='none', 
                     label="%s" %(labels[x]))
        plt.xlabel("Time")
        plt.ylabel("Current Individual Count")
        
        plt.legend()
        plt.savefig("Clinical Dynamics.png", dpi=300)
        plt.show()
        plt.clf()
        plt.close()
            
        
            
        
        
    def node_edge_annotation(self): #Initial node & edge annototation.
        
        '''Node Annotation:
        Nodes are classified as either "S" (Susceptible), "E" (Exposed), "I" (Infected), "R" (Recovered) or "D" (Dead)
        
        Under another classification system for transmission dynamics, we have
        Nodes during the "I" (infected) stage classified as: 
        "T"- Transmitting (Infectious), "NT" - Non-Transmitting, but still infected, 
        
        
        Based on severity of the disease we have a third classification system (applies for "E" & "I" nodes):
        "A"- Asymptomatic/Mild, "NS"- Non-severe (Moderate), "S"- Severe, "H"- Hospitalised (Severe Cases), "QS"- Self-Quarantined (Severe Cases), 
        "QNS"- Self-Quarantined (Non-Severe Cases)
        '''
        
        nex.set_node_attributes(self.SmWorldGr, 'UNASSIGNED', 'type')
        nex.set_node_attributes(self.SmWorldGr, 'UNASSIGNED', 'transtate')

        L= list(self.SmWorldGr) #List of nodes
        for n in L:
            ch=ran.random()
            self.SmWorldGr.nodes[n]['state']= 'S' #The whole population is initially susceptible.
            self.SmWorldGr.nodes[n]['t_SEIR']= (0,ch) #Chronicles absolute time (day) at which node attained it's current SEIRD status
            self.SmWorldGr.nodes[n]['source']= 'L' #The whole population is initially comprised of local nodes
        
        self.starterpack= ran.choices(L, k = self.int_caseload) #Returns the indices of the intiall cases at time t=0
        
        for n in self.starterpack:
            
            ch1=ran.random()
            ch2=ran.random()
            
            self.SmWorldGr.nodes[n]['state']= 'I' #Starter pack kids are initially infected.
            self.SmWorldGr.nodes[n]['transtate']= 'T' #Starter pack kids are initially transmitting.
            self.SmWorldGr.nodes[n]['t_trans']= (0,ch1) 
            '''Chronicles absolute time (day) at which node attained it's current transmission status
            Second co-ordinate stores a random number associated with the switching from T to NT'''
            self.SmWorldGr.nodes[n]['t_type']= (0,ch2) 
            ch= ran.random()
            if ch<=self.p_asymp:        #Individual is classified as aymptomatic
                self.SmWorldGr.nodes[n]['type']= 'A'
                
                #Chronicles absolute time (day) at which node attained it's current clinical type
            elif (ch>self.p_asymp and ch<= (self.p_asymp +self.p_sev)): #Individual will go on to develop a severe case
                self.SmWorldGr.nodes[n]['type']= 'S'
            else:
                self.SmWorldGr.nodes[n]['type']= 'NS'
            
        '''Edge Annotation
           Edges b/w "S" & "T" are tagged suspectible ("s")
        '''
        
        '''t_SEIR=nex.get_node_attributes(self.SmWorldGr, 't_SEIR') #Tuple of time at wich node last updated it's SEIR value & random number
        t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')'''
        state=nex.get_node_attributes(self.SmWorldGr, 'state')
        tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')
        
        for n in self.starterpack:
            for r in self.SmWorldGr.neighbors(n):
                if(state[r]=="S"): #If neighbors are susceptible.
                    self.SmWorldGr[n][r]== "s"
                    
        #print("Bokaro")
        
        '''Setting up various trackers'''

        self.infektion=copy.copy(self.starterpack)       #Stores indices of infected nodes at any given time step.
        self.susceptible=[]     #Stores indices of susceptible nodes at any given time step.
        self.exposed=[]         #Stores indices of exposed nodes at any given time step.
        self.recovered=[]       #Stores indices of recovered nodes at any given time step.
        self.dead=[]            #Stores indices of putrefying nodes at any given time step.
        
        self.foreigner=[]       # Register of foreign citizens, which  doesn't discriminate.
        
        '''Counts of stages of transmission'''
        
        self.trans=copy.copy(self.starterpack)       #Stores indices of infected, transmitting nodes at any given time step.
        self.nontrans=[]        #Stores indices of infected, non-transmitting nodes at any given time step.
        
        
        '''Counts clinical types'''
        
        self.asymp=[] #Stores indices of infected/exposed, asymptotic nodes at any given time step.
        self.sev=[]     #Stores indices of infected/exposed, severe nodes at any given time step.
        self.nonsev=[]  #Stores indices of infected/exposed, non-severe nodes at any given time step.
        self.qsev=[]  #Stores indices of infected/exposed, quarantined severe nodes at any given time step.
        self.hosp=[]  #Stores indices of infected/exposed, hospitalised nodes at any given time step.
        self.qnonsev=[]  #Stores indices of infected/exposed, quarantined non-severe nodes at any given time step.
        
        self.nonsevcrt=[] #Stores indices of infected non-severe nodes that are slated to become severe somewhere down the line.
        
        self.susceptible= list(self.SmWorldGr)
        #print("Roger")
        #print(self.susceptible[10:36])
        #Updates various empty lists initialised above.
        
        for n in self.infektion:
            if n in self.susceptible: self.susceptible.remove(n)
            # At t=0, all nodes are either infectious or susceptible
            if typo[n]=='A':
                self.asymp.append(n)
            elif typo[n]=='S':
                self.sev.append(n)
            elif typo[n]=='NS':
                self.nonsev.append(n)          

    
    def rates(self): #Assigns various transmission and clinical rates to patients. 
        
        self.gamma_e= math.log(2)/5.2        # Median incubation period of 5.2 days to symptom onset. (Kucharski et al)
        self.gamma_t= math.log(2)/6.5       
        # Median time of infectiousness (2.9 days) post symptom onset (Kucharski et al)
        # Imperial College data from the 16th of March assumes 
        self.gamma_q= 1.0/5.0 ; self.gamma_h = 1.0/5.0     
        # 5 days after symptom onset. (WHO Report). Also serves as time taken to hospitalise a patient post symptom onset.
        
        self.p_sev=0.19         #Percentage of severe cases (WHO Report).
        self.p_asymp=0.5       # Percentage that remain asymptomatic/mild throughout the course of their disease. (ECDC Report)
        
        
        self.gamma_a= 1.0/9.0  # Recovery in 9 days from symptom onset to recovery in mild/asymptomati cases. (WHO Report)
        self.gamma_s= 1.0/10.6       #Post-hospitalisation recovery time of 10.6. (WHO Report)
        self.hospbeds= (0.55*self.n)/1000
        #India has 0.55 hosital beds per thousand of the population.
        
        self.hospbeds=4000
        
        self.R0= 2.4 #Probability of transmission per edge between infected and suspectible persons 
        ''' Reproduction Rate =Beta/Gamma =3.0 (Kucharski et al)
        Gamma (Mean Recovery Rate)= 1/(0.19*28.6+ 0.81*14)'''
        
        self.beta_a= self.R0/(self.k*(math.log(2)/self.gamma_t)*(self.p_asymp + 2*(1-self.p_asymp)))
        
        #Probability of transmission per infected person (asymptomatic) per contact per unit time.
        #Probablitiy of transmission for 'NS' (Moderate) & 'S' (Severe) cases assumed to be twice that of asymptomatic.
        
        #print("Boka:\t %f" )
        
        self.int_caseload=2 #Initial number of infections (at t=0)
        
        self.p_critical=0.3         # 30% of severe cases turn critical
        self.death_hos_crt=0.5      # 50% chance critical patients will die will hopitalised.
        self.death_qs_crt=1         # 100% chance non-hospitalised critical care patients die.
        
        self.P_crit= 0.4
        
        '''Percentage of people who turn severe from moderate condition (from total moderate pool). Based on the WHO-China Joint Mission data,
           which claims 10% of all moderate cases end up in death. Assuming that 30% of all severe cases turn critical
           and that 80% of such critical cases overall end in death, we end at a figure of 40% severe cases turning critical'''
           
        self.gamma_h= 1.0/5.0
        # ASSUMPTION: MEDIAN OF 5.O DAYS FOR NS CASE TO CONVERT TO S.
        
        
obj=COVID_19_Basic()

if __name__ == '__main__':
    obj.controlpanel2()
        