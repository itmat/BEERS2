from collections import namedtuple
from beers.cluster import BaseCounts
import itertools
import math
import numpy as np

class BridgeAmplificationStep:

    def __init__(self, parameters):
        self.parameters = parameters
        self.cycles = parameters["bridge_amplification_cycles"]
        self.snp_rate = parameters["bridge_amplification_snp_percentage"]/100


    def execute(self, cluster_packet):
        print(f"Cluster Pkt size: {len(cluster_packet.clusters)}")
        for cycle in range(1,self.cycles + 1):
            for cluster in cluster_packet.clusters:
                if cycle < 4:
                    snps = self.determine_snps(cluster, cycle)
                    print(snps)
                #g = cluster.base_counts.G
                #a = cluster.base_counts.A
                #t = cluster.base_counts.T
                #c = cluster.base_counts.C
                #updated_base_counts = {"G":g, "A":a, "T":t, "C":c}
                #for ng,na,nt,nc in itertools.zip_longest(g,a,t,c, fillvalue=0):

    def determine_snps(self, cluster, cycle):
        snp_count = math.floor(len(cluster.molecule.sequence) * self.snp_rate * 2 ** cycle)
        print(f" Cycle: {cycle}, Seq Length: {len(cluster.molecule.sequence)} Snp count: {snp_count}")
        snps = []
        if snp_count > 0:
            snps = np.random.choice(range(0,len(cluster.molecule.sequence)), size=snp_count)
        return sorted(snps)






