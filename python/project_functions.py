import pickle
import numpy as np
import copy
from sklearn.metrics import roc_auc_score
import pandas as pd
import subprocess
import seaborn as sns
import os
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.metrics import roc_curve



def plot_concordance(df,xaxis,measure,ylabel,ylimits,title):

       
    my_palStrip = {"Coverage": "lightskyblue", "No Coverage": "orangered"}
    my_palBox = {"Coverage": "lightskyblue", "No Coverage": "orangered"}
 

    ax=sns.boxplot(data=df,x=xaxis,y=measure,order=['Coverage','No Coverage'],linewidth=0,boxprops=dict(alpha=0.2,linewidth=0,),width=0.3,\
    fliersize=0,palette=my_palBox)
    plt.setp(ax.collections, alpha=.2)

    ax=sns.swarmplot(data=df,x=xaxis,y=measure,order=['Coverage','No Coverage'],linewidth=0,s=7.5,alpha=0.6,palette=my_palStrip)



    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.spines['left'].set_edgecolor('grey')
    ax.spines['bottom'].set_edgecolor('grey')
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

    plt.xticks(size=14)
    plt.xlabel('',size=14)
    plt.yticks(size=14)
    plt.ylabel(ylabel,size=14)
    plt.ylim(ylimits)

    plt.title(title,size=14)


def plot_AUC(df,xaxis,measure,ylabel,ylimits,title):
    
    my_palStrip = {"ILAE1": "green", "ILAE2+": "red"}
    my_palBox = {"ILAE1": "green", "ILAE2+": "red"}
 

    ax=sns.boxplot(data=df,x=xaxis,y=measure,order=['ILAE1','ILAE2+'],linewidth=0,boxprops=dict(alpha=0.2,linewidth=0,),width=0.3,\
    fliersize=0,palette=my_palBox)
    plt.setp(ax.collections, alpha=.2)

    ax=sns.swarmplot(data=df,x=xaxis,y=measure,order=['ILAE1','ILAE2+'],linewidth=0,s=7.5,alpha=0.6,palette=my_palStrip)



    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.spines['left'].set_edgecolor('grey')
    ax.spines['bottom'].set_edgecolor('grey')
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)

    plt.xticks(size=14)
    plt.xlabel('',size=14)
    plt.yticks(size=14)
    plt.ylabel(ylabel,size=14)
    plt.ylim(ylimits)

    plt.title(title,size=14)
    