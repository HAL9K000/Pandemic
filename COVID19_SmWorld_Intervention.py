# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 16:34:45 2020

@author: Koustav
"""

import networkx as nex
import numpy as np
import matplotlib.pyplot as plt
import os
import random as ran
import math
import copy
from COVID19_SmWorld_Basic import COVID_19_Basic


class COVID_19_Int(COVID_19_Basic):
    
    def __init__(self):
        
        p=0.1 #Probability of rewiring each edge
        self.k=20 #Number of closest neighbours per node
        self.n=50000 #Number of nodes in the small world graph to begin with
        self.SmWorldGr= nex.watts_strogatz_graph(self.n, self.k, p)
        #Created Small World Graph.
        #self.attributes()
        self.rates() #Assigns various transmission and clinical rates to patients.
        
        self.rates2()       #Assigns additional ratez associated with interventions in the country.
        
        self.time= 200 #Stores the number of days for which the simulation is carried out.
        
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
        
        
        self.official_count=[]      #Official count of cases
        self.int_arrival_quar=[]        #List of numbers of people quarantined at each step.
        self.inf_int_size=[]
        
        self.cumulative_inf=[]
        
        self.str="Trial XIII"
        
    def controlpanel(self):
        
        
        self.node_edge_annotation() #Annotates nodes & edges initially (at t=0)
        
        for t in range(0, self.time):
            #continue
            self.time_evolution(t)
            print("Hashim")
        
        self.statistics()
        #Time to crunch the final numbers
        
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
        '''Contact tracing for international arrivals'''
        if( t >= self.t_ph2):
            self.int_contact_tracing(t)
        '''Contact tracing for local severe cases'''
        if( t >= self.t_cont):
            self.loc_contact_tracing(t)
        '''Updation from susceptible to exposed'''
        self.susceptible_updater(t)
        '''International arrivals occur at the end of the turn'''
        self.int_arrivals(t)
        
        self.stat_gen(t)         #Compiling data for statistics'''
        
        
    def int_arrivals(self, t):
        # Updates new infected international arrivals that have passed the turnstiles unchecked, at each time step.
        
        z1=0 ; z2 = 0
        if (t<= self.t_ph1 and  t>0):
            # We are currently in Phase I.
            
            z1 = round(self.ph1_u[0]*t) - round(self.ph1_u[0]*(t-1))
            # Number of undetected new arrivals added to the network at this time step.
            
            z2=round(self.ph1_d[0]*t) - round(self.ph1_d[0]*(t-1))
            
        elif (t<= self.t_ph2 and t>self.t_ph1):
            # We are in Phase II
            
            z1 = round(self.ph2_u[0]*t) - round(self.ph2_u[0]*(t-1))
            # Number of undetected new arrivals added to the network at this time step.
            
            z2=round(self.ph2_d[0]*t) - round(self.ph2_d[0]*(t-1))
            
        elif (t<= self.t_ph3 and t>self.t_ph2):
            # We are in Phase III
            
            z1 = round(self.ph3_u[0]*t) - round(self.ph3_u[0]*(t-1))
            # Number of undetected new arrivals added to the network at this time step.
            
            z2=round(self.ph3_d[0]*t) - round(self.ph3_d[0]*(t-1))
            
        new_inf= ran.choices(self.susceptible, k=z1)
        #Choosing indices of susceptible nodes that will randomly turn into infected nodes.
        
        self.int_arrival_quar.append(z2)
        
        for n in new_inf:
            
            self.susceptible.remove(n)      #Removing node from list of susceptible nodes.
            self.SmWorldGr.nodes[n]['source']='F'       #Designate source of infection as foreign.
            self.SmWorldGr.nodes[n]['t_source']=(t, ran.random())
            self.SmWorldGr.nodes[n]['t_SEIR']=(t, ran.random())
            
            self.foreigner.append(n)        #Appending n to list of foreigners in the system.
            
            ch=ran.random()
            a=math.log(2)/self.incubation
            b=math.log(2)/self.infectious
            
            if(ch <= (a/(a+b))):
                # Foreign node turns exposed.
                self.SmWorldGr.nodes[n]['state']='E'
                self.exposed.append(n)
            else:
                # Foreign node turns infected
                self.SmWorldGr.nodes[n]['state']='I'
                self.SmWorldGr.nodes[n]['transtate']='T'        #Infected node is transmitting.
                self.SmWorldGr.nodes[n]['t_trans']=(t, ran.random())
                
                self.trans.append(n)
                self.infektion.append(n)
            
            '''Finally to designate type of infection'''
            
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
                if( ch2 <= self.p_ns_s):
                    #NS case will go on to become severe later.
                    self.nonsevcrt.append(n)
                    ch3=ran.random()
                    self.SmWorldGr.nodes[n]['t_type']=(t,ch3, ran.random())
                            
            '''Third value in tuple is assigned only to NS cases that will turn to S and used to determine when these individuals 
            will end up in a 'S' case'''
        
        print("Undetected Arrivals at Time %d are %d" %(t, z1))
        print("Detected Arrivals at Time %d are %d" %(t, z2))
            
    def int_contact_tracing(self, t):
        #Begins starting at Phase II. This contact tracing module is only applicable to international arrivals.
        
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type') 
        state=nex.get_node_attributes(self.SmWorldGr, 'state') 
        
        if( t == self.t_ph2):
            #Day on which phase II kicks off.
            '''On this day, all foreigners currently in 'QS', 'QNS' or 'H' states undergo thorough contact tracing.
            Note that Asymptomatic foreigners are never reported and their corresponding secondary infections aren't either.'''
            new_quar=[]
            #Stores list of neigbouring nodes that have already undergone +ve contact testing in this round and are in quarantine/hospitalised
            for n in self.foreigner:
                
                if( typo[n]== "QNS" or typo[n]== "QS" or typo[n]== "H"):
                  print("Roman")
                    
                  if((t - t_typo[n][0])>=1):
                    print("California Soul")  
                    for r in self.SmWorldGr.neighbors(n):
                        #Neighbors of a foerign infected node undergo contact testing with a certain frequency.
                        
                        #if (state[r]== 'I' or state[r]== 'E' and r not in new_quar):
                        if (state[r]== 'I' and r not in new_quar):
                            #Neighbours must be infected.
                            ch=ran.random()
                            if(typo[r]== "S"):
                                #Person has a severe case
                                if (ch > self.p_contact):
                                        # Person wasn't detected by contact testing.
                                        continue
                                
                                if ((len(self.hosp)< self.hospbeds) and state[r]=="I"):
                                    '''NEWLY DIAGNOSED SEVERE CASES ARE OFFERED HOSPITAL BEDS, IF ANY ARE AVALIABLE, IF THEY SHOW SYMPTOMS.'''
                                    self.SmWorldGr.nodes[r]['type']= 'H'
                                    self.hosp.append(r)
                                else:
                                    #Due to lack of hospital beds, person enters self-quarantine
                                    self.SmWorldGr.nodes[r]['type']= 'QS'
                                    self.qsev.append(r)
                            
                                new_quar.append(r)
                                self.SmWorldGr.nodes[r]['t_type']= (t, ran.random())
                                #Time at which patient entered hospital (is diagnosed). Second co-ord associated with date of Recovery/Death
                            
                            elif(typo[r]== "NS"):
                                # Less severe cases to self-quarantine.
                                if (ch > self.p_contact):
                                        # Person wasn't detected by contact testing.
                                        continue
                                print("Irfan")
                                self.SmWorldGr.nodes[r]['type']= 'QNS'
                                self.qnonsev.append(r)
                                self.SmWorldGr.nodes[r]['t_type']= (t, ran.random())
                                new_quar.append(r)
                        
                                if (len(t_typo[r])==3):
                                    # If these NS individuals are slated to go rogue.
                                    ch2=t_typo[r][2]
                                    self.SmWorldGr.nodes[r]['t_type']= (t, ran.random(), ch2)
                                    
        elif( t>self.t_ph2):
            '''Phase II contact tracing policy for foreigners in full force. All newly (in this turn) hospitalised/quarantined foreigners undergo contact tracing.'''
            new_quar=[]
            #Stores list of neigbouring nodes that have already undergone +ve contact testing and are in quarantine/hospitalised
            
            for n in self.foreigner:
                
                if(typo[n]== "QNS" or typo[n]== "QS" or typo[n]== "H"):
                    
                    
                    if((t - t_typo[n][0])==1):
                      #Contact tracing and round up of suspected cases occurs a day after foreigner is first reported as positive.
                      '''This condition needs to be made in the event there happens to exist two foreign nodes adjacent to each other.
                      If t_typo[n][0] == 0, then say if one of the nodes gets quarantined, then a search will immediately occur for other 
                      infected adjacent node, in which case, the other foreign node might also get put in quarantine. In this event,
                      no contact tracing will be done for the other node that gets put into quarantine'''
                      print("Beastie Boys")
                      for r in self.SmWorldGr.neighbors(n):
                        #Neighbors of a foerign infected node undergo contact testing with a certain frequency.
                        
                        #if (state[r]== 'I' or state[r]== 'E' and r not in new_quar):
                        if (state[r]== 'I' and r not in new_quar):
                            #Neighbours must be infected.
                            print("Yabadaba")
                            ch=ran.random()
                            if(typo[r]== "S"):
                                #Person has a severe case
                                if (ch > self.p_contact):
                                        # Person wasn't detected by contact testing.
                                        continue
                                print("Doolittle")
                                if ((len(self.hosp)< self.hospbeds) and state[r]=="I"):
                                    '''NEWLY DIAGNOSED SEVERE CASES ARE OFFERED HOSPITAL BEDS, IF ANY ARE AVALIABLE, IF THEY SHOW SYMPTOMS.'''
                                    self.SmWorldGr.nodes[r]['type']= 'H'
                                    self.hosp.append(r)
                                else:
                                    #Due to lack of hospital beds, person enters self-quarantine
                                    self.SmWorldGr.nodes[r]['type']= 'QS'
                                    self.qsev.append(r)
                            
                                new_quar.append(r)
                                self.SmWorldGr.nodes[r]['t_type']= (t, ran.random())
                                #Time at which patient entered hospital (is diagnosed). Second co-ord associated with date of Recovery/Death
                            
                            elif(typo[r]== "NS"):
                                # Less severe cases to self-quarantine.
                                if (ch > self.p_contact):
                                        # Person wasn't detected by contact testing.
                                        print("Eggman")
                                        continue
                                print("Bloodclot")
                                self.SmWorldGr.nodes[r]['type']= 'QNS'
                                self.qnonsev.append(r)
                                self.SmWorldGr.nodes[r]['t_type']= (t, ran.random())
                                new_quar.append(r)
                        
                                if (len(t_typo[r])==3):
                                    # If these NS individuals are slated to go rogue.
                                    ch2=t_typo[r][2]
                                    self.SmWorldGr.nodes[r]['t_type']= (t, ran.random(), ch2)
                        
                        
            
            
            
            
            
    def loc_contact_tracing(self,t):
         '''Begins on day 33 (18th March when counted from 14th Feb). This contact tracing module is only applicable 
         to local severe cases.'''
         
         typo=nex.get_node_attributes(self.SmWorldGr, 'type')
         t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type') 
         state=nex.get_node_attributes(self.SmWorldGr, 'state')
         source=nex.get_node_attributes(self.SmWorldGr, 'source')
         
         if (t == self.t_cont):
             new_quar=[]
             for n in self.sev:
                '''Only local severe cases are considered'''
             
                if(source[n]== "F"):
                    #Only local severe cases are considered in this module.
                    continue
                if(typo[n]== "H"):
                    # As per current protocol only hospitalised local cases are tessted and may undergo contact tracing.
                    
                    if((t - t_typo[n][0])>=1):
                      '''Must have benn hospitalised for  a day or more.'''
                      print("Firhous")
                      for r in self.SmWorldGr.neighbors(n):
                      #Neighbors of a foerign infected node undergo contact testing with a certain frequency.
                        
                        if (state[r]== 'I' and r not in new_quar):
                            #Neighbours must be infected.
                            ch=ran.random()
                            print("Bahaus")
                            if(typo[r]== "S"):
                                #Person has a severe case
                                if (ch > self.p_contact):
                                        # Person wasn't detected by contact testing.
                                        continue
                                print("Marley")
                                if ((len(self.hosp)< self.hospbeds) and state[r]=="I"):
                                    '''NEWLY DIAGNOSED SEVERE CASES ARE OFFERED HOSPITAL BEDS, IF ANY ARE AVALIABLE, IF THEY SHOW SYMPTOMS.'''
                                    self.SmWorldGr.nodes[r]['type']= 'H'
                                    self.hosp.append(r)
                                else:
                                    #Due to lack of hospital beds, person enters self-quarantine
                                    self.SmWorldGr.nodes[r]['type']= 'QS'
                                    self.qsev.append(r)
                            
                                new_quar.append(r)
                                self.SmWorldGr.nodes[r]['t_type']= (t, ran.random())
                                #Time at which patient entered hospital (is diagnosed). Second co-ord associated with date of Recovery/Death
                            
                            elif(typo[r]== "NS"):
                                # Less severe cases to self-quarantine.
                                if (ch > self.p_contact):
                                        # Person wasn't detected by contact testing.
                                        continue
                                
                                self.SmWorldGr.nodes[r]['type']= 'QNS'
                                self.qnonsev.append(r)
                                self.SmWorldGr.nodes[r]['t_type']= (t, ran.random())
                                new_quar.append(r)
                                print("Hamzah")
                        
                                if (len(t_typo[r])==3):
                                    # If these NS individuals are slated to go rogue.
                                    ch2=t_typo[r][2]
                                    self.SmWorldGr.nodes[r]['t_type']= (t, ran.random(), ch2)
             
         elif( t > self.t_cont):
             new_quar=[]
             for n in self.sev:
                '''Only local severe cases are considered'''
             
                if(source[n]== "F"):
                    #Only local severe cases are considered in this module.
                    continue
                if(typo[n]== "H"):
                    # As per current protocol only hospitalised local cases are tessted and may undergo contact tracing.
                    
                    if((t - t_typo[n][0])==1):
                      '''Contact tracing as in previous cases, only takes place a day after the hospitalisation.'''
                      print("Yolanda")
                      for r in self.SmWorldGr.neighbors(n):
                      #Neighbors of a foerign infected node undergo contact testing with a certain frequency.
                        
                       if (state[r]== 'I' and r not in new_quar):
                            print("Gerbil")
                            #Neighbours must be infected.
                            ch=ran.random()
                            if(typo[r]== "S"):
                                #Person has a severe case
                                if (ch > self.p_contact):
                                        # Person wasn't detected by contact testing.
                                        continue
                                print("Von Lent")
                                if ((len(self.hosp)< self.hospbeds) and state[r]=="I"):
                                    '''NEWLY DIAGNOSED SEVERE CASES ARE OFFERED HOSPITAL BEDS, IF ANY ARE AVALIABLE, IF THEY SHOW SYMPTOMS.'''
                                    self.SmWorldGr.nodes[r]['type']= 'H'
                                    self.hosp.append(r)
                                    #print("Von Lent")
                                else:
                                    #Due to lack of hospital beds, person enters self-quarantine
                                    self.SmWorldGr.nodes[r]['type']= 'QS'
                                    self.qsev.append(r)
                                    #print("Von Lent")
                            
                                new_quar.append(r)
                                self.SmWorldGr.nodes[r]['t_type']= (t, ran.random())
                                #Time at which patient entered hospital (is diagnosed). Second co-ord associated with date of Recovery/Death
                            
                            elif(typo[r]== "NS"):
                                # Less severe cases to self-quarantine.
                                if (ch > self.p_contact):
                                        # Person wasn't detected by contact testing.
                                        print("War Chest")
                                        continue
                                print("Rastafari")
                                self.SmWorldGr.nodes[r]['type']= 'QNS'
                                self.qnonsev.append(r)
                                self.SmWorldGr.nodes[r]['t_type']= (t, ran.random())
                                new_quar.append(r)
                        
                                if (len(t_typo[r])==3):
                                    # If these NS individuals are slated to go rogue.
                                    ch2=t_typo[r][2]
                                    self.SmWorldGr.nodes[r]['t_type']= (t, ran.random(), ch2)
                        
                    
                    
                    
    def stat_gen(self, t):
        
        self.sus_size.append(len(self.susceptible))
        self.exp_size.append(len(self.exposed))
        
        self.trans_size.append(len(self.trans))
        self.rec_size.append(len(self.recovered))
        self.dead_size.append(len(self.dead))
        self.inf_int_size.append(len(self.foreigner))
        
        self.sev_size.append(len(self.sev))
        self.hosp_size.append(len(self.hosp))
        self.qsev_size.append(len(self.qsev))
        self.qnonsev_size.append(len(self.qnonsev))
        
        if( len(self.int_arrival_quar)< 14 and len(self.int_arrival_quar)>0):
            self.inf_size.append(len(self.infektion) + sum(self.int_arrival_quar))
        elif( len(self.int_arrival_quar)>= 14):
            self.inf_size.append(len(self.infektion) + sum(self.int_arrival_quar[-14:-1]))
            
            
        self.cumulative_inf.append(self.n - len(self.susceptible) + sum(self.int_arrival_quar))
        
        print("Number of susceptible persons:\t%d" %(len(self.susceptible)))
        print("Number of exposed persons:\t%d" %(len(self.exposed)))
        print("Number of infected persons:\t%d" %(len(self.infektion)))
        print("Cumulative number of infected international arrivals:\t%d" %(len(self.foreigner)))
        print("Number of infected (transmitting) persons:\t%d" %(len(self.trans)))
        
        #print("Total number of quarantined (on screening) international arrivals:\t%d" %(len(self.inf_size[-1])))
        print("Total Number of currently infected persons:\t%d" %(self.inf_size[-1]))
        
        
        print("Number of dead persons:\t%d" %(len(self.dead)))
        
        print("Number of severe cases:\t%d" %(len(self.sev)))
        print("Number of severe quarantined cases:\t%d" %(len(self.qsev)))
        
        #self.ouputgraph(t)
        
        
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
        
        cum_inf_size=np.array(self.cumulative_inf)
        
        
        
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
        
        if(os.path.isdir("Logs")==False):
            os.mkdir("Logs")
        os.chdir("Logs")
        
        tseries= np.array(list(range(0,self.time)))
        self.data= np.zeros([self.time, 12])
        
        self.data[:,0] = tseries
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
        
        self.data[:,11] = cum_inf_size
        
        '''Exporting data as CSV'''
        
        hd_txt= "Time, 1: Suscep Indv, 2: No. Exp Indv, 3: No. Inf Indv, 4: No. Cont Indv, 5: No. Rec Indv, 6: No. Dead Indv, 7: Sev Cases, 8: Hosp Cases, 9:Q Sev Cases, 10: QNS Cases, 11: Cumulative Inf Count " 
        np.savetxt('Pandemic.csv', self.data, delimiter=",", header=hd_txt,comments="#")
                   
        f=open('log_unabridged.txt','w')
        f.write(" K:\t%d, N:\t%d, Initial Caseload:\t %d, R0:\t%3.2f, p:\t %5.3f" %(self.k, self.n, self.int_caseload, self.R0, self.p))
        
                   
        os.chdir("../")
        
        #Plotting and saving data here.
        
        labels=["Time", "# Susceptible Individuals", "# Active Exposed Individuals", "# Active Infected Individuals", 
                "# Active Infected Individuals (Contagious)", "# Recovered Individuals",  "# Dead Individuals",
               "# Active Severe Cases", "# Active Hosp Individuals", "# Active Quarantined Severe Cases",
               "# Active Quarantined Non-Severe Cases", "Cumulative Infection Count"]
        
        for x in range(1,7):
            #Plotting basic SEIRD data.
            
            plt.plot(self.data[:,0], self.data[0:,x],  marker='o', markerfacecolor='none', 
                     label="%s" %(labels[x]))
        plt.xlabel("Time")
        plt.ylabel("Current Individual Count")
        
        plt.legend()
        plt.savefig("SEIRD.png", dpi=300)
        plt.show()
        plt.clf()
        plt.close()
    
        
        for x in range(6,10):
            
            plt.plot(self.data[:,0], self.data[0:,x],  marker='o', markerfacecolor='none', 
                     label="%s" %(labels[x]))
        plt.xlabel("Time")
        plt.ylabel("Current Individual Count")
        
        plt.legend()
        plt.savefig("Clinical Dynamics.png", dpi=300)
        plt.show()
        plt.clf()
        plt.close()
        
        
        plt.plot(self.data[:,0], self.data[0:,11],  marker='o', markerfacecolor='none', 
                     label="%s" %(labels[11]))
        plt.plot(self.data[:,0], self.data[0:,3],  marker='o', markerfacecolor='none', 
                     label="%s" %(labels[3]))
        plt.plot(self.data[:,0], self.data[0:,2],  marker='o', markerfacecolor='none', 
                     label="%s" %(labels[2]))
        plt.xlabel("Time")
        plt.ylabel("Infection Stats")
        
        plt.legend()
        plt.savefig("Cumulative Infection Count.png", dpi=300)
        plt.show()
        plt.clf()
        plt.close()

        
    def rates(self): #Assigns various transmission and clinical rates to patients. 
        
        self.incubation= math.log(2)/5.2        # Median incubation period of 5.2 days to symptom onset. (Kucharski et al)
        self.infectious= math.log(2)/6.5       
        # Median time of infectiousness (2.9 days) post symptom onset (Kucharski et al)
        # Imperial College data from the 16th of March assumes 
        self.lab_detection= 1.0/5.0     
        # 5 days after symptom onset. (WHO Report). Also serves as time taken to hospitalise a patient post symptom onset.
        
        self.p_sev=0.19         #Percentage of severe cases (WHO Report).
        self.p_asymp=0.33       # Percentage that remain asymptomatic/mild throughout the course of their disease. (ECDC Report)
        
        
        self.recoverymild= 1.0/14.0  # 14 days from symptom onset to recovery in mild/asymptomati cases. (WHO Report)
        self.recoverysev= 1.0/28.6       #Post-hospitalisation recovery time. (WHO Report)
        self.hospbeds= (0.55*self.n)/1000
        #India has 0.55 hosital beds per thousand of the population.
        self.hospbeds= 2000
        
        self.R0= 2.4 #Probability of transmission per edge between infected and suspectible persons 
        ''' Reproduction Rate =Beta/Gamma =3.0 (Kucharski et al)
        Gamma (Mean Recovery Rate)= 1/(0.19*28.6+ 0.81*14)'''
        
        self.p= self.R0/(self.k*(math.log(2)/self.infectious)*(self.p_asymp + 2*(1-self.p_asymp)))
        
        #Probability of transmission per infected person (asymptomatic) per contact per unit time.
        #Probablitiy of transmission for 'NS' (Moderate) & 'S' (Severe) cases assumed to be twice that of asymptomatic.
        
        #print("Boka:\t %f" )
        
        self.int_caseload=0 #Initial number of infections (at t=0)
        
        self.p_critical=0.3         # 30% of severe cases turn critical
        self.death_hos_crt=0.5      # 50% chance critical patients will die will hopitalised.
        self.death_qs_crt=1         # 100% chance non-hospitalised critical care patients die.
        
        self.p_ns_s= 0.4
        
        '''Percentage of people who turn severe from moderate condition (from total moderate pool). Based on the WHO-China Joint Mission data,
           which claims 10% of all moderate cases end up in death. Assuming that 30% of all severe cases turn critical
           and that 80% of such critical cases overall end in death, we end at a figure of 40% severe cases turning critical'''
           
        self.conversion= math.log(2)/5.0
        # ASSUMPTION: MEDIAN OF 5.O DAYS FOR NS CASE TO CONVERT TO S.
        
    def rates2(self):
        
        '''self.t_ph1=18       #Day number at which phase I ends (3rd March when starting from 14th Feb)
        self.t_ph2=28       #Day number at which phase II ends (13nd March when starting from 14th Feb) [Mandatory quarantine for some EU countries]'''
        
        self.t_ph1=18
        self.t_ph2=28
        
        
        '''self.t_ph3=37       #Day number at which phase III ends (22nd March starting from Valentine's Day) [Moratium on all international flights]
        self.t_ph4= 39      # Start day for Phase IV (24th March)
        
        self.t_cont=33      #Day number (18th March) from which all severe hospitalised cases are tested additionally.'''
        
        self.t_ph3=37
        self.t_ph4= 39
        self.t_cont=33
        
        '''Phase I: (14th Feb-3rd March). Thermal screening of passengers from incoming flights at Delhi & Mumbai airports from select countries
        No contact tracing program is initiated.'''
        
        dp=87000     #Daily number of incoming passagners.
        infp=0.00005 #Probability that any given passanger is infected/got exposed before or during the flight.
        pthermchecks=0.85       #Percentage of disease carriers coming from thermally screened category at the two airports
        pmumbdel=0.6            #Percentage of people landing at Delhi and Mumbai
        
        undet=dp*infp*self.t_ph1*(((((math.log(2)/self.incubation)/11.7)+ ((math.log(2)/self.infectious)/11.7)*(0.33/0.81))*pmumbdel +(1-pmumbdel))*0.85 +0.15)
        print("Undetected Entries in Phase I:\t%f" %(undet))
        det=dp*infp*self.t_ph1*((((math.log(2)/self.infectious)/11.7)*(0.48/0.81))*pmumbdel)*0.85
        print("Nullified COVID-19 cases on entry in Phase I:\t%f" %(det))
        
        
        self.caseload_ph1= int(round(det)) 
        
        self.undet_ph1= int(round(undet)) #Caseload (undetected interational arrival cases popping up) in phase I (14th Feb- 3rd March)
        
        m=float(self.caseload_ph1)/self.t_ph1
        #print(m)
        
        self.ph1_d=(m, 0)       #Stores the coeffiecients of the linear spline for new undetected cases.
        m=float(self.undet_ph1)/self.t_ph1
        #print(m)
        self.ph1_u=(m, 0)
        
        
        '''Phase II: (3rd March-13th March). Thermal screening of passengers from all incoming flights at all airports.
        Dedicated contact tracing program is initiated. Assume efficieny of contact tracing to be 70%.'''
        
        dp=70000
        infp=0.00015        #Possibility of being infected has shot up.
        undet=dp*infp*10*(((math.log(2)/self.incubation)/11.7)+ ((math.log(2)/self.infectious)/11.7)*(0.33/0.81))
        print("Undetected Entries in Phase II:\t%f" %(undet))
        det=dp*infp*10*(((math.log(2)/self.infectious)/11.7)*(0.48/0.81))
        print("Nullified COVID-19 cases on entry in Phase II:\t%f" %(det))
        
        self.caseload_ph2= int(round(det)) 
        
        self.undet_ph2= int(round(undet)) #Caseload (undetected cases popping up) in Phase II (3rd-13th March)
        
        m=float(self.caseload_ph2)/(float(self.t_ph2 - self.t_ph1))
        #print(m)
        self.ph2_d=(m, 0)       
        
        m=float(self.undet_ph2)/(self.t_ph2 - self.t_ph1)
        #print(m)
        self.ph2_u=(m, 0)       #Stores the coeffiecients of the linear spline for new undetected cases.
        
        '''Phase III: (13th March-22nd March). Thermal screening of passengers from all incoming flights at all airports.
        Increasing restrictions on overseas flights. Mandatory quarantine for all passengers disembarking from 
        Italy, Germany, France, Spain, Korea, China, Iran.
        Dedicated contact tracing program continues. Assume efficieny of contact tracing to be 70%.'''
        
        dp=30000
        infp=0.00015        
        qfrac= 0.5  #Percentage of incoming passengers that have been subjected to madatory quarantine.
        
        undet=dp*infp*9*(((math.log(2)/self.incubation)/11.7)+ ((math.log(2)/self.infectious)/11.7)*(0.33/0.81))*qfrac
        print("Undetected Entries in Phase III:\t%f" %(undet))
        det=dp*infp*9*(qfrac+((math.log(2)/self.infectious)/11.7)*(0.48/0.81)*(1-qfrac))
        print("Nullified COVID-19 cases on entry in Phase III:\t%f" %(det))
        
        self.caseload_ph3= int(round(det)) 
        
        self.undet_ph3= int(round(undet)) #Caseload (undetected cases popping up from international arrivals) in Phase III (13rd-22nd March)
        
        m=float(self.caseload_ph3)/(self.t_ph3 - self.t_ph2)
        #print(m)
        self.ph3_d=(m, 0)       
        
        m=float(self.undet_ph3)/(self.t_ph3 - self.t_ph2)
        self.ph3_u=(m, 0)  #Stores the coeffiecients of the linear spline for new undetected cases.
        #print(m)
        
        
        
        '''Setting parameters for contact testing'''
        
        self.p_contact=0.8
        
        
obj=COVID_19_Int()

if __name__ == '__main__':
    obj.controlpanel()
        
        