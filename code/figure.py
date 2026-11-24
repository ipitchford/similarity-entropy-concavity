"""Draw the exact five-state example. Dependencies: Matplotlib, NumPy."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def draw():
    a=12/17
    A,B=(1+a*a)/2,(1+2*a*a)/3
    x=np.linspace(0,1,401)
    u=x*A+(1-x)*a
    v=x*a+(1-x)*B
    gap=x*np.log(u/A)+(1-x)*np.log(v/B)
    berry,pine="#8b3f68","#267165"
    fig,(graph,ax)=plt.subplots(1,2,figsize=(7.15,2.6),gridspec_kw={'width_ratios':[1,2.05]},layout="constrained")
    left=[(0,.77),(0,.23)]
    right=[(1,.93),(1,.5),(1,.07)]
    for l in left:
        for r in right:
            graph.plot([l[0],r[0]],[l[1],r[1]],color='#b4b7b1',lw=1.6,zorder=1)
    for group,color,labels in [(left,berry,'12'),(right,pine,'345')]:
        for (xx,yy),label in zip(group,labels):
            graph.scatter(xx,yy,s=250,c=color,zorder=2)
            graph.text(xx,yy,label,ha='center',va='center',color='white',fontsize=9,zorder=3)
    graph.text(0,-.12,'Part P',color=berry,ha='center');graph.text(1,-.12,'Part Q',color=pine,ha='center')
    graph.set(xlim=(-.4,1.4),ylim=(-.28,1.13),aspect='equal');graph.axis('off')
    ax.fill_between(x,0,gap*1e4,color=berry,alpha=.10);ax.plot(x,gap*1e4,color=berry,lw=2)
    ax.scatter([.5],[gap[200]*1e4],c=pine,s=30,zorder=3)
    ax.set(xlabel='Weight assigned to P',ylabel='Mixing gap (10⁻⁴ nats)',xlim=(0,1),ylim=(0,3.5),xticks=[0,.25,.5,.75,1],yticks=[0,1,2,3])
    ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.25)
    fig.savefig(Path(__file__).resolve().parent.parent/'figures'/'mixing-example.pdf',bbox_inches='tight',pad_inches=.08)
    plt.close(fig)
if __name__=='__main__':
    draw()
