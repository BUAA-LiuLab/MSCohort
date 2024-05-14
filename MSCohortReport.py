

import numpy as np
import math
from math import pi
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.backends.backend_svg
from collections import Counter
from matplotlib.pyplot import MultipleLocator
import os



class CFunctionPlotChromatography:
    def __init__(self, inputDP):
        self.dp = inputDP

    def draw(self, figfolder):
        save_path = figfolder + IO_FILENAME_PIC_EXPORT[0]
        print("Generating image " + str(save_path) + ".png ......")
        self.init_fig()
        self.draw_Peakwidth(self.fig, self.gs, self.dp.myID.PSM29_PeakWidth)
        self.draw_RT_Peakwidth(self.fig, self.gs, self.dp.myID.PSM3_RT, self.dp.myID.PSM29_PeakWidth,self.dp.myID.R2_GRADIENT)

        plt.savefig(save_path + '.png', format='png')
        plt.close()

    def init_fig(self):
        fig = plt.figure()
        fig.set_size_inches(20, 5, forward=True)
        self.gs = GridSpec(6, 20)  
        self.fig = fig

    def draw_Peakwidth(self, fig, gs, peakwidth):
        ax1 = fig.add_subplot(gs[0:6, 0:8])
        duration = [round(x, 2) for x in peakwidth]
        ax1.hist(duration, edgecolor='k', color="steelblue", bins=np.arange(0, max(duration) + 1, 0.1))
        max_duration = max(duration) if max(duration)<3.0 else 3.0
        ax1.set_xlim(0, max_duration)
        x_locator = MultipleLocator(0.2)
        ax1.xaxis.set_major_locator(x_locator)
        ax1.set_xlabel("Full Width (min.)", fontsize=15, fontname="arial", fontweight='bold')
        ax1.set_ylabel("Frequency", fontsize=15, fontname="arial", fontweight='bold')
        plt.tick_params(labelsize=13)
        plt.title("Histogram of Full Width of Chromatography", fontsize=15, fontname="arial", fontweight='bold')
        chrom_median = np.median(duration)
        plt.axvline(chrom_median, 0, 1, color="r", linestyle="dashed")  
        ymaxlim = ax1.get_ylim()  
        ax1.text(chrom_median + 0.5, 0.85 * ymaxlim[1], ha="left", fontsize=12,
                 s="Median of full width \n= " + str(chrom_median))  

        notail = 0
        for i in self.dp.myID.PSM29_PeakWidth:
            if i <= 1:
                notail += 1
        tail_score = round((notail / len(self.dp.myID.PSM29_PeakWidth)) * 100,
                           2)  

        
        plt.axvline(1, 0, 0.4, color="r", linestyle="dashed")  
        ax1.annotate("", xytext=(1, 0.3 * ymaxlim[1]), xy=(1.5, 0.3 * ymaxlim[1]), 
                     arrowprops=dict(arrowstyle="->, head_length = 1.5, head_width = .5", facecolor="red"))
        ax1.text(1.4, 0.35 * ymaxlim[1], ha="center", fontsize=12,
                 s="Proportion of peptides\nwith full width > 1 min. = " +
                   str(round(100 - tail_score, 2)) + "%")

     
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

    def draw_RT_Peakwidth(self,fig, gs,RT,Peakwidth,Gradient):
        ax2 = fig.add_subplot(gs[0:6, 12:20])
        peakwidth = [round(x, 2) for x in Peakwidth]
        timesum = Gradient[0]
        if len(peakwidth) > 100000:
            step = int(len(RT)/5)
        elif len(peakwidth) > 50000:
            step = int(len(RT) / 3)
        else:
            step = 1

        ax2.scatter(RT[::step], peakwidth[::step], s=2)
        ax2.set_xlim(0, timesum + 3)
        ax2.set_xticks(np.arange(0, timesum + 3, int(0.1 * round(timesum))))
        ax2.set_xticklabels(np.arange(0, timesum + 3, int(0.1 * round(timesum))), fontsize=13,
                               fontname="Times New Roman")
        ax2.set_xlabel("Retention Time (min.)", fontsize=15, fontname="Times New Roman", fontweight='bold')

        ax2.set_ylim(0, 1.1)
        major_ticks = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
       
        ax2.set_yticks(major_ticks)
        ax2.set_yticklabels(major_ticks, fontsize=13, fontname="Times New Roman")
        
        ax2.set_ylabel("Peakwidth(min)", fontsize=15, fontname="Times New Roman", fontweight='bold')

        
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)


class CPlotMS2Scans:
    def __init__(self, inputDP, inputDataMS1, inputDataMS2):
        self.dp = inputDP
        self.dataMS1 = inputDataMS1
        self.dataMS2 = inputDataMS2

    def draw(self, figfolder):
        save_path = figfolder + IO_FILENAME_PIC_EXPORT[5]
        print("Generating image " + str(save_path) + ".png ......")
        self.init_fig()
        ax_marg_x, ax_marg_y, ax_joint = self.draw_duration(self.fig, self.gs, self.dp.myCYCLE.LIST_CYCLES_TIME,
                                                            self.dp.myCYCLE.LIST_CYCLES_N_MS2)
        self.set_ax(ax_marg_x, ax_marg_y, ax_joint, self.dp.myCYCLE.LIST_CYCLES_TIME, self.dp.myCYCLE.LIST_CYCLES_N_MS2)
        self.draw_acc(self.fig, self.gs, self.dp.myMS1.INDEX_SCAN_TIME_MS1, self.dp.myMS2.INDEX_SCAN_TIME_MS2,
                      self.dataMS1.INDEX_RT, self.dp.myCYCLE.LIST_CYCLES_N_MS2)
        plt.savefig(save_path + '.png', format='png')
        plt.close()

    def init_fig(self):
        self.fig = plt.figure()
        self.fig.set_size_inches(20, 6, forward=True)
        self.gs = GridSpec(4, 24)

    def draw_duration(self, fig, gs, ms3_duration, ms3_pic):
        ax_joint = fig.add_subplot(gs[1:4, 0:8]) 
        ax_marg_x = fig.add_subplot(gs[0, 0:8]) 
        ax_marg_y = fig.add_subplot(gs[1:4, 8:10])

        ax_joint.scatter(ms3_duration, ms3_pic, s=10)
        x_ticks1 = np.arange(0, max(max(ms3_duration) + 0.5, 3.6), 0.25)
        y_ticks = np.arange(0, max(max(ms3_pic) + 3, 51), 5)
        ax_marg_x.hist(ms3_duration, rwidth=1, alpha=0.8, edgecolor='k', bins=x_ticks1)  
        ax_marg_y.hist(ms3_pic, orientation="horizontal", rwidth=1, alpha=0.8, edgecolor='k', bins=y_ticks)

        return ax_marg_x, ax_marg_y, ax_joint

    def set_ax(self, ax_marg_x, ax_marg_y, ax_joint, ms3_duration, ms3_pic):
        x_ticks = np.arange(0, max(max(ms3_duration) + 0.5, 3.6), 0.5)
        y_ticks = np.arange(0, max(max(ms3_pic) + 3, 51), 5)
        ax_joint.set_xticks(x_ticks)
        ax_joint.set_yticks(y_ticks)
        ax_joint.grid(which='major', alpha=0.5)

        y1 = ax_marg_x.get_yticks()
        x1 = ax_marg_y.get_xticks()
        ax_marg_x.set_xticks(x_ticks)
        ax_marg_x.set_yticks(np.arange(0, max(y1 + 20), max(y1) / 2))
        ax_marg_y.set_xticks(np.arange(0, max(x1 + 20), max(x1) / 2))
        ax_marg_y.set_yticks(y_ticks)

        plt.setp(ax_marg_x.get_xticklabels(), visible=False)  
        plt.setp(ax_marg_y.get_yticklabels(), visible=False)

        ax_joint.tick_params("x", labelsize=13, labelrotation=315)
        ax_joint.tick_params("y", labelsize=13)
        ax_marg_x.tick_params(labelsize=13)
        ax_marg_y.tick_params(labelsize=13, labelrotation=315)

        ax_joint.set_xlabel("Time of One Cycle (Sec.)", fontname="arial", fontsize=15, fontweight='bold')
        ax_joint.set_ylabel("MS2 Scans in One Cycle", fontname="arial", fontsize=15, fontweight='bold')

        ax_marg_y.spines['top'].set_visible(False)
        ax_marg_y.spines['right'].set_visible(False)
        ax_marg_x.spines['top'].set_visible(False)
        ax_marg_x.spines['right'].set_visible(False)

    def draw_acc(self, fig, gs, ms1_scantime, ms2_scantime, ms2_RT, N_MS2):

        ax_acc2 = fig.add_subplot(gs[0:4, 14:24])
        ms2_RT = [x / 60 for x in ms2_RT]
        ms3_RT = ms2_RT[0:-1]
        ax_acc2.scatter(ms3_RT, N_MS2, s=3) 

        
        timesum1 = np.sum(ms1_scantime)
        timesum2 = np.sum(ms2_scantime)
        timesum = int((timesum2 + timesum1) / 60)

        ax_acc2.set_xlim(0, timesum + 3)
        if timesum > 10:
            ax_acc2.set_xticks(np.arange(0, timesum + 3, int(0.1 * round(timesum))))  #
            ax_acc2.set_xticklabels(np.arange(0, timesum + 3, int(0.1 * round(timesum))), fontsize=13)
        else:
            ax_acc2.set_xticks(np.arange(0, timesum + 6, 1))  #
            ax_acc2.set_xticklabels(np.arange(0, timesum + 6, 1), fontsize=13)
        ax_acc2.set_xlabel("Retention Time (min.)", fontsize=15, fontname="arial", fontweight='bold')

        
        ax_acc2.set_ylim(0, max(N_MS2) + 2)
        major_ticks2 = np.arange(0, max(N_MS2) + 2, 10, dtype=int)
        minor_ticks2 = np.arange(0, max(N_MS2) + 2, 2)
        ax_acc2.set_yticks(major_ticks2)
        ax_acc2.set_yticks(minor_ticks2, minor=True)
        ax_acc2.set_yticklabels(major_ticks2, fontsize=13)
        ax_acc2.set_ylabel("Acquired MS2 Scans in One Cycle", fontsize=15, fontname="arial", fontweight='bold')
     
        ax_acc2.grid(which='both')
        ax_acc2.grid(which='minor', alpha=0.2)
        ax_acc2.grid(which='major', alpha=0.5)

        ax_acc2.spines['top'].set_visible(False)
        ax_acc2.spines['right'].set_visible(False)

