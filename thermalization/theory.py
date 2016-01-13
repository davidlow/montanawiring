#  Theory for montana wiring, thermalization between stages
#  Copyright (C) 2015 David Low, Nowack Lab
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>
import numpy as np;
import matplotlib.pyplot as plt;

def powerbobbin(deltaT,    # > 0
                 wirelen,
                 k_ins,     #thermal conductivity of polyesterimide W/(m K)
                 k_cu       = 385, #thermal conductivity of cu, W/(m K)
                 k_sty      = 1,#thermal conductivity of stycast W/(m K)
                 wrapnum    = 3,#num wrap loom around bobbin in 1 layer
                 d_bob      = 14.478e-3, #diameter of bobbin in meters
                 d_wire     = 79.9e-6, # diameter of wire in meters 
                 insthick   = 17.78e-6, # thickness of insulation in meters 
                 layerspace = 80e-6,
                 ):
    ns    = range(1, numturns(wirelen, wrapnum, d_bob) + 1);
    UtotA = sum( [ invUnA(n, k_ins, k_sty, wrapnum, wirelen, d_bob, d_wire, 
                          insthick, layerspace)**(-1) for n in ns ]);
    return UtotA * deltaT;

def maketable(mat, deltaTs, lengths, diameters):
    table = [];
    for d in diameters:
        row = [];
        for l in lengths:
            row.append( 
            ( powerlength(deltaTs[0], l, larger(mat[0], mat[1]), d), 
              powerlength(deltaTs[1], l, larger(mat[1], mat[2]), d) 
            )
            )
        table.append(row);
    return table;


def table():
    cu = [400, 600, 300]; # all of these are in W/ (m K)
    pb = [50, 20, 1];
    cn = [11, 10, 0.4];

    allmaterials = [cu, pb, cn];
    matnames     = ['copper', 'phosphor bronze', 'constantan'];
    deltaTs = [300-60, 60-4];

    lengths = [1, .5, .2, .1, .05, .01];    # in m

    gauges = [30, 32, 34, 36, 38, 40];      # in AWG

    diameters = [.127e-3 * 92**((36-n)/39) for n in gauges]; #in m

    for mat in range(len(matnames)):
        print('');
        print(matnames[mat] + '  ' + str(lengths));
        table = maketable(allmaterials[mat], deltaTs, lengths, diameters);
        for r in range(len(table)):
            row = '| {:2.2f} AWG '.format(gauges[r]);
            for l in table[r]:
                row = row + '| {:2.4f} mW '.format(l[1]*1e3);
            row = row + '|'
            print(row);


def larger (a,  b):
    return a if a>b else b;

def powerlength(deltaT,    # > 0
                 wirelen,
                 k_wire = 21, #thermal conductivity of wire in W/(m K)
                 d_wire = 79.9e-6, # in meters
                 ):
    P = k_wire * (np.pi * (d_wire)*(d_wire) / 4.0 ) * deltaT / wirelen;
    return P;


def numturns(wirelen, wrapnum, d_bob):
    # number of layers you need to wrap a wire of wirelen around a bobbin
    # of diameter d_bob if you can have wrapnum wraps per layer.
    n = np.ceil(wirelen / (np.pi * d_bob * wrapnum));
    return int(n);

def wrapwirearea(n, wirelen, wrapnum, d_bobbin, d_wire):
    # how much surface area the wrapped wire has per layer
    return ( ((wirelen / (np.pi * d_bobbin)) - (n-1)) * 
             (np.pi * d_bobbin) * (np.pi * d_wire)
           );

def invUnA(n,
           k_ins, 
           k_sty,
           wrapnum,
           wirelen,
           d_bob,
           d_wire,
           insthick,    # 17.78e-6, thickness of insulation in meters
           layerspace,  # 80e-3 distance between layers in meters
           ):
    invUins= (insthick/k_ins) / (
                wrapwirearea(n, wirelen, wrapnum, d_bob, d_wire));
    invUsty= (n * layerspace / k_sty)  / (
                wrapwirearea(n, wirelen, wrapnum, d_bob, d_wire));
    return invUins + invUsty

def run():
    lengths = 10**(np.linspace(-3,0,100));
    P1 = np.array([powerlength(270, l)       for l in lengths]);
    P2 = np.array([powerbobbin(270, l, 0.1, 
                           layerspace=2e-3, 
                           k_sty=0.064)  
                                         for l in lengths]); #k_ins for 300K
    P3 = np.array([powerlength(27,  l)       for l in lengths]);
    P4 = np.array([powerbobbin(27,  l, 0.05, 
                           layerspace=2e-3, 
                           k_sty=0.064)  
                                         for l in lengths]); #k_ins for 30K
    P5 = np.array([powerbobbin(27,  l, 0.05, 
                           layerspace=5e-3, 
                           k_sty=0.064)  
                                         for l in lengths]); #k_ins for 30K
    fig = plt.figure();
    ax  = fig.add_subplot(1,1,1);
    ax.loglog(lengths, P1*1e3*50, marker='o', 
            label='P1 = 300K -> 30K bobbin');
    ax.loglog(lengths, P2*1e3*50, marker='^', 
            label='P2 = 30K bobbin wrap');
    ax.loglog(lengths, P3*1e3*50, marker='s', 
            label='P3 = 30K -> 4K bobbin');
    ax.loglog(lengths, P4*1e3*50, marker='x', 
            label='P4 = 4K bobbin wrap');
    ax.legend(loc='upper left');
    ax.annotate('P1 < P2, P2 > P3, P3 < P4', xy=(-20,10), 
            xycoords='axes points', horizontalalignment='right',
            verticalalignment='bottom', fontsize=15);
    plt.xlabel('Length of Wire for Given Component (m)');
    plt.ylabel('Max Possible Power Transfer(mW)');
    fig.savefig('powertransver_vs_wirelen_forbobbins.pdf',
                bbox_inches='tight');

if __name__ == '__main__':
    table();

