#!/usr/bin/env python
# Created: Sat Jul 14 14:20:55 2001
# Last changed: Time-stamp: <01/07/16 13:45:19 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas
# File: xkMeans.py

import string, re
import os, sys, commands
from Tkinter import *
from Numeric import *
from Bio.Tools.Clustering import kMeans 
from tkFileDialog import askopenfilename, asksaveasfilename
import Pmw

"""
    This is a small xwrapper for Jeff's kMeans module. It displays the
    ongoing change of cluster sizes and allows easy retrieval of the observations
    per cluster (Button-1 on cluster box).
    Input is a cvs file

    """

class xkMeans:
    def __init__(self, parent = None):
        self.Init_vars()
        self.Init_graphics(parent)

    def Init_vars(self):
        self.x0, self.y0, self.xdiff, self.xwidth = 50, 300, 50, 20
        self.k = 3
        self.values = {
            'format':'CSV',
            'k':2
            }
        
    def Init_graphics(self, parent = None):
        self.parent = parent
        self.main = Frame(parent)
        self.buttons = Frame(self.main, width = 30)
        self.Init_Buttons()
        self.scanvas = Pmw.ScrolledCanvas(self.main)
        self.canvas = self.scanvas.interior()
        self.canvas.configure(bg = 'white')

        self.main.pack(fill = BOTH, expand = 1)
        self.buttons.pack(side = LEFT, fill = Y)
        self.scanvas.pack(side = RIGHT, fill = BOTH, expand = 1)

    def Init_Buttons(self):
        f = self.buttons
        b = []
        width = 20
        self.info = Label(f, width = width, fg = 'blue', relief = 'ridge')
        b.append(self.info)
        b.append(Button(f, text = 'Open', command = self.ReadDialog, width = width))
        self.k_entry = Pmw.EntryField(f,
                                      labelpos = 'w',
                                      value = str(self.values['k']),
                                      label_text = 'k:',
                                      validate = {'validator' : 'integer'},
                                      )
        b.append(self.k_entry)
                            
        b.append(Button(f, text = 'Run', command = self.Run, width = width))
        b.append(Button(f, text = 'Exit', command = self.Quit, width = width))
        
        for widget in b: widget.pack(side=TOP)

    def Quit(self):
        sys.exit(0)

    def Info(self, txt):
        self.info.configure(text = txt)
        self.info.update()
        
    def ReadDialog(self):
        file = askopenfilename()
        if not file: return
        self.Read(file)
        
    def Read(self, file):
        self.Info('Reading %s' % os.path.split(file)[-1])
        read_func = getattr(self, 'Read%s' % self.values['format'])
        read_func(file)
        self.Info('%d observations' % self.n_observations)
        
    def ReadCSV(self, file):
        fid = open(file)
        line = fid.readline().strip()
        self.attributes = line.split(',')
        data = []
        raw = []
        start_field = (self.attributes[0].lower() == 'label')

        n=0
        lines = fid.readlines()
        print len(lines)
        for line in lines:
            self.Info('%d observation' % n)
            data.append([float(x) for x in line.split(',')[start_field:]])
            raw.append(line)
            n+=1

            
        self.data = data
        self.raw_data = raw
        self.M = array(data)
        self.n_observations = len(data)
        self.n_atrributes = len(self.attributes) - start_field

    def Gather(self, clusters):
        dict = {}
        for i in range(0, self.k): dict[i] = 0
        for i in clusters: dict[i]+= 1
        return dict
    
    def GatherMembers(self, clusters):
        dict = {}
        for i in range(self.k): dict[i] = []
        for i in range(len(clusters)):
            dict[clusters[i]].append(i)

        return dict
    
    def Update_fn(self, i, centroids, clusters):
        " visualize the change of the cluster sizes, pass to kMeans.cluster"

        x0, y0, xdiff, xwidth = self.x0, self.y0, self.xdiff, self.xwidth
        data_points = len(clusters)
        k = self.k
        c = self.canvas
        
        # first round has empty clusters
        if clusters[0] == None: return

        # clear the canvas
        c.delete('deletable')
        self.Info('Round %2.2d' % i)

        dict = self.Gather(clusters)

        # draw the new cluster sizes (easier and ?faster? than resizing
        # old rectangles
        x = x0
        for i in range(k):
            x += xdiff
            y = y0 - dict[i]*100.0/data_points
            c.create_rectangle(x,y0-1,x+xwidth, y, fill = 'magenta', tag=('deletable', 'boxes'))
            c.create_text(x+xwidth/2, y0 + 20, text = str(dict[i]), tag = 'deletable')
            
        c.update()


    def Run(self):
        self.ChooseData()
        
    def ChooseData(self):
        #self.top = Toplevel(self.parent)
        self.top = Frame(self.buttons, bg = 'dimgrey')
        self.top.pack(side = TOP)
        choose = Pmw.RadioSelect(self.top,
                                 buttontype = 'checkbutton',
                                 orient = 'vertical',
                                 hull_borderwidth = 2,
                                 hull_relief = 'ridge',
	)
	choose.pack(side = TOP, expand = 1, padx = 10, pady = 10)

	for text in (self.attributes):
	    choose.add(text)
            choose.invoke(text)

        Button(self.top, text = 'Run', command = self.RunCluster).pack(side = TOP)
        self.choose = choose
        
    def RunCluster(self):
        self.top.destroy()
        selected_data = [self.choose.index(x) for x in self.choose.getcurselection()]
        M = take(self.M, selected_data,1)
        self.k = int(self.k_entry.get())

        c = self.canvas
        x0, y0, xdiff, xwidth = self.x0, self.y0, self.xdiff, self.xwidth
        c = self.canvas
        c.create_line(x0+xdiff/2, y0, x0 + xdiff*(self.k+1), y0)
        c.update()
        self.scanvas.resizescrollregion()
        
        if self.Cluster(M):
            self.FinishClusters()
        
    def FinishClusters(self):
        c = self.canvas
        k = self.k
        
        c.itemconfigure('boxes', fill = 'green3')
        c.update()

        cluster_ids = list(c.find_withtag('boxes'))
        cluster_ids.sort()
        
        dict = self.GatherMembers(self.clusters)


        for i in range(k):
            tk_id = cluster_ids[i]
            members = dict[i]
            c.tag_bind(tk_id, '<1>',  lambda x, i=i, m=members,
                       f = self.ShowCluster:f(i, m))
        c.focus()
        
    def ShowCluster(self, i, members):
        top = Toplevel(self.parent)
        tid = Pmw.ScrolledText(top)
        tid.pack(fill = BOTH, expand = 1)
        res = 'Cluster %d (%d members)\n' % (i, len(members))
        for member in members:
            res +=self.raw_data[member]
        tid.insert(END, res)
            
    def Cluster(self, M):
        k = self.k
        try:
            self.centroids, self.clusters = kMeans.cluster(M, k, update_fn = self.Update_fn)

        except:
            self.canvas.delete(ALL)
            txt = '\t\t\tempty cluster -> bad choice of initial centroids\n\t\t\ttry again'
            self.canvas.create_text(self.x0,self.y0, text = txt, fill = 'red')
            return 0
        
        return 1

if __name__ == '__main__':
    win = Tk()
    win.geometry('1000x600')
    main = xkMeans(parent = win)
    if len(sys.argv) == 2:
        file = sys.argv[1]
        main.Read(file)


    win.mainloop()
