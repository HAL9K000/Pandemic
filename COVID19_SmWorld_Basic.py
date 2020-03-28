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
import copy

class COVID_19_Basic():
    
    def __init__(self):
        
        p=0.006 #Probability of rewiring each edge
        self.k=16 #Number of closest neighbours per node
        self.n=10000 #Number of nodes in the small world graph to begin with
        self.SmWorldGr= nex.watts_strogatz_graph(self.n, self.k, p)
        #Created Small World Graph.
        #self.attributes()
        self.rates() #Assigns various transmission and clinical rates to patients.
        
        self.time= 300 #Stores the number of days for which the simulation is carried out.
        
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
        
        print("Hello")
        
        for n in range(1,10):
            print(ran.random())
        
        
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
        '''Updation from susceptible to exposed'''
        self.susceptible_updater(t)
        
        self.stat_gen()         #Compiling data for statistics'''
        
        
    def exposed_updater(self,t):
        print("Genda")
        
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
                    print(" Asymptomatic Stoad\t:%d" %(t_typo[n][0]))
            
            #print("Toad:\t%d" %(t_SEIR[n][0]))
            #print("Noad:\t %f" %(t_SEIR[n][1]))
            if( ch> math.exp(-self.incubation*(t- t_SEIR[n][0]))):
                #Person becomes infected from exposed
                
                self.SmWorldGr.nodes[n]['state']= 'I'
                self.SmWorldGr.nodes[n]['transtate']= 'T'       
                #Person becomes transmitting ( or only after onset of symptoms)
                self.SmWorldGr.nodes[n]['t_SEIR']= (t, ran.random())
                self.SmWorldGr.nodes[n]['t_trans']= (t, ran.random())
                
                self.infektion.append(n)
                self.trans.append(n)
                self.exposed.remove(n)
                
                tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
                t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
                
                if(typo[n]=='A'):
                    print(" Asymptomatic Stoad\t:%d" %(t_typo[n][0]))
                    print(" Asymptomatic Time Transition\t:%d\t%f" %(t_trans[n]))
                    print(" Asymptomatic Transition\t:%s %s" %(tstate[n], n))
                
    def susceptible_updater(self, t):
        
        #t_SEIR=nex.get_node_attributes(self.SmWorldGr, 't_SEIR') #Time at wich node last updated it's SEIR value.
        #t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
        state=nex.get_node_attributes(self.SmWorldGr, 'state')   #Dict of all nodes with their SEIR status
        tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        #t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')
        
        '''Checking for transmission of disease from contagious infected persons to susceptible persons'''
        
        new_exp=[]
        #To keep  track of newly minted "E" nodes in each call of this function
        
        for n in self.trans:
            
            #state=nex.get_node_attributes(self.SmWorldGr, 'state')
            #Needed to ensure that state variable is updated and that r belongs only to 'S' & not some 'E'
            
            
            
            p=self.p
            if( typo[n]=='H'):
                p=0.33*self.p
                '''ASSUMPTION: IF PATIENT IS STILL CONTAGIOUS WHEN HOSPITALISED, TRANSMISSION RATE DROPS TO 33% 
                of it's value'''
                
            elif( typo[n]=='QS' or typo[n]=='QNS'):
                p=0.5*self.p
                '''ASSUMPTION: IF PATIENT IS STILL CONTAGIOUS WHEN HOSPITALISED, TRANSMISSION RATE DROPS TO 50% 
                of it's value'''
            
                
            for r in self.SmWorldGr.neighbors(n):
                if state[r]=='S' and r not in new_exp:
                    #If a neighbour is susceptible and that neighbor hasn't already tranformed to an 'E' node
                    ch=ran.random()
                    if( ch <= p):
                        #Susceptible individual has become exposed.
                        self.SmWorldGr.nodes[r]['state']='E'
                        self.SmWorldGr.nodes[r]['t_SEIR']=(t, ran.random())
                        print("Susceptible Transformation to exposed %d\t%f" %(self.SmWorldGr.nodes[r]['t_SEIR']))
                        self.exposed.append(r)
                        
                        new_exp.append(r)
                        
                        print(len(self.susceptible))
                        if r not in self.susceptible:
                            print("Eureka")
                            #print(state[r])
                            #print(r)
                        self.susceptible.remove(r)
                        
                        ch1=ran.random()
                        ch2=ran.random()
                        self.SmWorldGr.nodes[r]['t_type']=(t,ch2)
                        
                        '''Choosing whether the exposed will go on to develop severe, non-severe, or asymptomatic cases'''
                        
                        if( ch1<=self.p_asymp):        #Individual is classified as aymptomatic
                            self.SmWorldGr.nodes[r]['type']= 'A'
                            self.asymp.append(r)
                
                            #Chronicles absolute time (day) at which node attained it's current clinical type
                        elif (ch1>self.p_asymp and ch1<= (self.p_asymp +self.p_sev)): #Individual will go on to develop a severe case
                            self.SmWorldGr.nodes[r]['type']= 'S'
                            self.sev.append(r)
                        else:
                            self.SmWorldGr.nodes[r]['type']= 'NS'
                            self.nonsev.append(r)
                            print("SOja")
                        
            
        
                
            
    def infected_updater(self, t):
        #Updating clinical dynamics, transmission cases and recovery/death toll at each step.
        
        t_SEIR=nex.get_node_attributes(self.SmWorldGr, 't_SEIR') #Time at wich node last updated it's SEIR value.
        t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
        state=nex.get_node_attributes(self.SmWorldGr, 'state')   #Dict of all nodes with their SEIR status
        tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')
        
        '''Updating potential T cases to NT based on prob check'''
        
        for n in self.trans:
            if( state[n] != 'I'):
                #If node is not infectious then it can no longer be contagious.
                self.trans.remove(n)
            else:
                ch=t_trans[n][1]
                #print("Andim %1.4f" %(ch))
                if( ch> math.exp(-self.infectious*(t- t_trans[n][0]))):
                    #Person becomes non-transmitting and lists are updated accordingly.
                    
                    self.nontrans.append(n)
                    self.SmWorldGr.nodes[n]['transtate']= 'NT'
                    self.SmWorldGr.nodes[n]['t_trans']= (t, 0)
                    #Chronicles absolute time (day) at which node attained it's current transmission status
                    print("Lensmart:\t Inf %d\t Cont %d" %(len(self.infektion), len(self.trans)))
                    self.trans.remove(n)
                    print("Walmart:\t Inf %d\t Cont %d" %(len(self.infektion), len(self.trans)))
                    tstate=nex.get_node_attributes(self.SmWorldGr, 'transtate')
                    t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')
                
                    if(typo[n]=='A'):
                        print(" Asymptomatic Stoad\t:%d" %(t_typo[n][0]))
                        print(" Asymptomatic Time Transition\t:%d\t%f" %(t_trans[n]))
                        print(" Asymptomatic Transition\t:%s %s" %(tstate[n], n))
        
        t_trans=nex.get_node_attributes(self.SmWorldGr, 't_trans')            
        
                    
        '''Checking if a particular infectious has gotten back positive diagnostic for COVID-19
        & is hospitalised/quarantined as the case may be'''
        
        for n in self.infektion:
            
            if (typo[n]== "S" or typo[n]== "NS"):
                ch=t_SEIR[n][1]
                if( ch> math.exp(-self.lab_detection*(t- t_SEIR[n][0]))):
                    #Takes a median of five days from the onset of the symptoms for a person to be diagnosed.
                    
                    if(typo[n]== "S"):
                        #Person has a severe case
                        if (len(self.hosp)< self.hospbeds):
                            '''ONLY NEWLY DIAGNOSED SEVERE CASES ARE OFFERED HOSPITAL BEDS, IF ANY'''
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
                        #Time at which patient is diagnosed.
        
        typo=nex.get_node_attributes(self.SmWorldGr, 'type')
        t_typo=nex.get_node_attributes(self.SmWorldGr, 't_type')                
        
                        
        '''Checking if non-severe quarantined (QNS) persons have recovered.
        N.B: We don't accont for death in non-severe cases, not yet.'''
        
        
        
        for n in self.qnonsev:
            print("Baloo")
            print(typo[n])
            ch=t_typo[n][1]
            if( ch> math.exp(-self.recoverymild*(t- t_SEIR[n][0]))):
                #Takes a mean of fourteen days from the onset of the symptoms for a person to recover.
                
                #Recovered patient can no longer be transmitting/non-transmitting
                if(tstate[n]== 'T'):
                    self.SmWorldGr.nodes[n]['t_trans']= (t,0)
                    self.trans.remove(n)
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
                
        for n in self.sev:
            
            if typo[n]=='H':
                #Hospitalised patient.
                ch=t_typo[n][1]
                if( ch> math.exp(-self.recoverysev*(t- t_typo[n][0]))):
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
                        self.trans.remove(n)
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
                if( ch> math.exp(-self.recoverysev*(t- t_typo[n][0]))):
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
                        self.trans.remove(n)
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

        for n in self.asymp:
            
            if( state[n] != 'I' ):
                #If person is not infectious yet, they cannot recover.
                continue
            
            #print("Asympt t_trans: %d\t%f" %(t_trans[n]))
            ch=t_typo[n][1]
            if( ch> math.exp(-self.recoverymild*(t- t_SEIR[n][0]))):
                #Takes a mean of fourteen days from the onset of the symptoms for an asymptomatic to recover. (ASSUMPTION)
                
                print("Asymptototic: %s" %(n))
                print("Asympt t_type: %d\t%f" %(t_typo[n]))
                print("Asympt typo: %s" %(typo[n]))
                print("Asympt t_SEIR: %d\t%f" %(t_SEIR[n]))
                print("Asympt SEIR State: %s" %(state[n]))
                print("Asympt t_trans: %d\t%f" %(t_trans[n]))
                print("Asympt transtate: %d\t%f" %(t_trans[n]))
                
                
                #Recovered patient can no longer be transmitting/non-transmitting
                if(tstate[n]== 'T'):
                    self.SmWorldGr.nodes[n]['t_trans']= (t,0)
                    self.trans.remove(n)
                    #Removing from transmitter list
                
                self.SmWorldGr.nodes[n]['transtate']= ''
                
                self.SmWorldGr.nodes[n]['type']= ''
                self.SmWorldGr.nodes[n]['state']='R'
                self.SmWorldGr.nodes[n]['tstate']=(t,0)
                self.SmWorldGr.nodes[n]['t_type']= (t, 0)
                #Time at which patient recovered.
                
                self.recovered.append(n)
                self.asymp.remove(n)
                self.infektion.remove(n)
                
                
    def stat_gen(self):
        
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
        
        print("Number of susceptible persons:\t%d" %(len(self.susceptible)))
        print("Number of exposed persons:\t%d" %(len(self.exposed)))
        print("Number of infected persons:\t%d" %(len(self.infektion)))
        print("Number of infected (transmitting) persons:\t%d" %(len(self.trans)))
        print("Number of dead persons:\t%d" %(len(self.dead)))
        
        print("Number of severe cases:\t%d" %(len(self.sev)))
        print("Number of severe quarantined cases:\t%d" %(len(self.qsev)))
        
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

        if(os.path.isdir("Trial IV")==False):
            os.mkdir("Trial IV")
        os.chdir("Trial IV")
        
        #Log making
        
        if(os.path.isdir("Logs")==False):
            os.mkdir("Logs")
        os.chdir("Logs")
        
        tseries= np.array(list(range(0,self.time)))
        self.data= np.zeros([self.time, 11])
        
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
        
        '''Exporting data as CSV'''
        
        hd_txt= "Time, 1: Suscep Indv, 2: No. Exp Indv, 3: No. Inf Indv, 4: No. Cont Indv, 5: No. Rec Indv, 6: No. Dead Indv, 7: Sev Cases, 8: Hosp Cases, 9:Q Sev Cases, 10: QNS Cases " 
        np.savetxt('Pandemic.csv', self.data, delimiter=",", header=hd_txt,comments="#")
                   
        os.chdir("../")
        
        #Plotting and saving data here.
        
        labels=["Time", "# Susceptible Individuals", "# Exposed Individuals", "# Infected Individuals (Total)", 
                "# Infected Individuals (Contagious)", "# Recovered Individuals",  "# Dead Individuals",
               "# Severe Cases", "# Hosp Individuals", "# Quarantined Severe Cases",
               "# Quarantined Non-Severe Cases"]
        
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
            
        
            
        
        
    def node_edge_annotation(self): #Initial node & edge annototation.
        
        '''Node Annotation:
        Nodes are classified as either "S" (Susceptible), "E" (Exposed), "I" (Infected), "R" (Recovered) or "D" (Dead)
        
        Under another classification system for transmission dynamics, we have
        Nodes during the "I" (infected) stage classified as: 
        "T"- Transmitting (Infectious), "NT" - Non-Transmitting, but still infected, 
        
        
        Based on severity of the disease we have a third classification system (applies for "E" & "I" nodes):
        "A"- Asymptomatic, "NS"- Non-severe, "S"- Severe, "H"- Hospitalised (Severe Cases), "QS"- Self-Quarantined (Severe Cases), 
        "QNS"- Self-Quarantined (Non-Severe Cases)
        '''
        
        L= list(self.SmWorldGr) #List of nodes
        for n in L:
            ch=ran.random()
            self.SmWorldGr.nodes[n]['state']= 'S' #The whole population is initially susceptible.
            self.SmWorldGr.nodes[n]['t_SEIR']= (0,ch) #Chronicles absolute time (day) at which node attained it's current SEIRD status
        
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
                    
        print("Bokaro")
        
        '''Setting up various trackers'''

        self.infektion=copy.copy(self.starterpack)       #Stores indices of infected nodes at any given time step.
        self.susceptible=[]     #Stores indices of susceptible nodes at any given time step.
        self.exposed=[]         #Stores indices of exposed nodes at any given time step.
        self.recovered=[]       #Stores indices of recovered nodes at any given time step.
        self.dead=[]            #Stores indices of putrefying nodes at any given time step.
        
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
        
        self.susceptible= list(self.SmWorldGr)
        print("Roger")
        print(self.susceptible[2211:2231])
        #Updates various empty lists initialised above.
        
        for n in self.infektion:
            self.susceptible.remove(n)
            # At t=0, all nodes are either infectious or susceptible
            if typo[n]=='A':
                self.asymp.append(n)
            elif typo[n]=='S':
                self.sev.append(n)
            elif typo[n]=='NS':
                self.nonsev.append(n)          

    
    def rates(self): #Assigns various transmission and clinical rates to patients. 
        
        self.incubation= math.log(2)/5.2        # Median incubation period of 5.2 days to symptom onset. (Kucharski et al)
        self.infectious= math.log(2)/2.9        # Median time of infectiousness (2.9 days) post symptom onset (Kucharski et al)
        self.lab_detection= 1.0/5.0     # 5 days after symptom onset. (WHO Report). Also serves as time taken to hospitalise a patient post symptom onset.
        
        self.p_sev=0.19         #Percentage of severe cases (WHO Report).
        self.p_asymp=0.03       # Percentage that remain asymptomatic throughout the course of their disease. (ECDC Report)
        
        
        self.recoverymild= 1.0/14.0  # 14 days from symptom onset to recovery in mild cases. (WHO Report)
        self.recoverysev= 1.0/28.6       #Post-hospitalisation recovery time. (WHO Report)
        self.hospbeds= (0.55*self.n)/1000
        #India has 0.55 hosital beds per thousand of the population.
        
        self.R0= 3.0 #Probability of transmission per edge between infected and suspectible persons 
        ''' Reproduction Rate =Beta/Gamma =3.0 (Kucharski et al)
        Gamma (Mean Recovery Rate)= 1/(0.19*28.6+ 0.81*14)'''
        
        self.p= self.R0/(self.k*(math.log(2)/self.infectious))
        #Probability of transmission per infected person per contact per unit time
        
        self.int_caseload=1 #Initial number of infections (at t=0)
        
        self.p_critical=0.3         # 30% of severe cases turn critical
        self.death_hos_crt=0.5      # 50% chance critical patients will die will hopitalised.
        self.death_qs_crt=1         # 100% chance non-hospitalised critical care patients die.
        
        
obj=COVID_19_Basic()

if __name__ == '__main__':
    obj.controlpanel()
        
        
        
        
        