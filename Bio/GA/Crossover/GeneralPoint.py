"""
Generalized N-Point Crossover.

For even values of N, perform N point crossover 
  (select N/2 points each in the two genomes, and alternate)
For odd values of N, perform symmetric N+1 point crossover
  (select N/2 points for both genomes)
  
N-Point introduction (my notation):
    genome 1:    A-----B-----C-----D-----E-----F-----G
    genome 2:    a=====b=====c=====d=====e=====f=====g
    
    2-point, symmetric (points=1):
                 A-----B-----C- 1 -D-----E-----F-----G
                 a=====b=====c= 2 =d=====e=====f=====g
    returns: (ABCdefg, abcDEFG)

    2-point, asymmetric (points=2):
                 A-----B- 1 -C-----D-----E-----F-----G
                 a=====b=====c=====d=====e= 2 =f=====g
    returns: (ABfg, abcdeCDEFG)

and for the drastic (n can be arbitrary to the length of the genome!):
    12-point, symmetric (points=11):
                 A- 1 -B- 2 -C- 3 -D- 4 -E- 5 -F- 6 -G
                 a= 7 =b= 8 =c= 9 =d= 10 e= 11 f= 12 g
    returns: (AbCdEfG, aBcDeFg)
    (note that points=12 will yield the same result, but 11
     may be somewhat faster)
        
"""
# standard modules
import random

class GeneralPointCrossover(object):
    """Perform n-point crossover between genomes at some defined rates.

       Ideas on how to use this class:
           - Call it directly ( construct, do_crossover )
           - Use one of the provided subclasses
           - Inherit from it:
               * replace _generate_locs with a more domain 
                 specific technique
               * replace _crossover with a more efficient 
                 technique for your point-count
    """
    def __init__(self, points, crossover_prob = .1):
        """Initialize to do crossovers at the specified probability.
        """
        self._crossover_prob = crossover_prob

        self._sym     = points % 2 # odd n, gets a symmetry flag
        self._npoints = (points + self._sym)//2 # (N or N+1)//2
    
    def do_crossover(self, org_1, org_2):
        """Potentially do a crossover between the two organisms.
        """
        new_org = ( org_1.copy(), org_2.copy() )
        
        # determine if we have a crossover
        crossover_chance = random.random()
        if crossover_chance <= self._crossover_prob:
            
            # pre-compute bounds (len(genome))
            bound  = (len(new_org[0].genome), len(new_org[1].genome))
            
            mbound = min(bound)
            # can't have more than 0,x_0...x_n,bound locations
            if (self._npoints == 0 or self._npoints > mbound-2):
                self._npoints = mbound-2
                
            y_locs = []
            # generate list for the shortest of the genomes
            x_locs = self._generate_locs( mbound )

            if (self._sym != 0):  
                y_locs = x_locs
            else:
                # figure out which list we've generated 
                # for, and generate the other
                if (mbound == bound[0]):
                    y_locs = self._generate_locs( bound[1] )
                else:
                    y_locs = x_locs
                    xlocs  = self._generate_locs( bound[0] )
              
            # copy new genome strings over
            tmp = self._crossover(0, new_org, (x_locs,y_locs))
            new_org[1].genome = self._crossover(1, new_org, (x_locs,y_locs))
            new_org[0].genome = tmp

        return new_org

    def _generate_locs(self, bound):
        """Generalized Location Generator:
            
           arguments:
               bound (int)   - upper bound 
            
           returns: [0]+x_0...x_n+[bound]
             where n=self._npoints-1
               and 0 < x_0 < x_1 ... < bound
        """
        results = []
        for increment in range(self._npoints):
            x = random.randint(1,bound-1)
            while (x in results):  # uniqueness
                x = random.randint(1,bound-1)
            results.append( x )
        results.sort()             # sorted
        return [0]+results+[bound] # [0, +n points+, bound]
            
    def _crossover( self, x, no, locs ):
        """Generalized Crossover Function:
            
           arguments: 
               x (int)        - genome number [0|1]
               no (organism,organism)
                              - new organisms
               locs (int list, int list)
                              - lists of locations, 
                                [0, +n points+, bound]
                                for each genome (sync'd with x)

            return type: sequence (to replace no[x])
        """
        s = no[ x ].genome[ :locs[ x ][1] ]
        for n in range(1,self._npoints):
            # flipflop between genome_0 and genome_1
            mode = (x+n)%2
            # _generate_locs gives us [0, +n points+, bound]
            #  so we can iterate: { 0:loc(1) ... loc(n):bound }
            t = no[ mode ].genome[ locs[mode][n]:locs[mode][n+1] ]
            if (s): s = s + t
            else:   s = t
        return s

    
class TwoCrossover(GeneralPointCrossover):
    """Helper class for Two Point crossovers.
    
    Offers more efficient replacements to the GeneralPoint framework
    for single pivot crossovers
    """
    def _generate_locs(self, bound):
        """Replacement generation.
         
        See GeneralPoint._generate_locs documentation for details
        """
        return [0, random.randint(1,bound-1), bound]
    
    def _crossover( self, x, no, locs ):
        """Replacement crossover
         
           see GeneralPoint._crossover documentation for details
        """
        y = (x+1)%2
        return no[x].genome[ : locs[x][1] ] + no[y].genome[ locs[y][1] : ]

class InterleaveCrossover(GeneralPointCrossover):
    """Demonstration class for Interleaving crossover.
    
    Interleaving:  AbCdEfG, aBcDeFg
    """
    def __init__(self,crossover_prob=0.1):
        GeneralPointCrossover.__init__(self,0,crossover_prob)
    
    def _generate_locs(self,bound):
        return range(-1,bound+1)
    
    def _crossover( self, x, no, locs ):
        s = no[ x ].genome[ 0:1 ]
        for n in range(1,self._npoints+2):
            mode = ( x+n )%2
            s += no[ mode ].genome[ n:n+1 ]
        return s+no[mode].genome[self._npoints+3:]
