# Compiled by Randal S. Olson
# Description: Useful functions for statistics and data parsing that I've written or run across during my research


# Handy function for parsing lists from a csv reader into a string-referenced dictionary
# Input: parsed list from csv.reader()
# Returns: dictionary of data lists, list of header strings in same order as csv file
def parse_csv_data(data_list):
    data_dict = {}
    headers = []
    first_row = True
    
    for row in data_list:
        
        # first row has headers
        if first_row:
            first_row = False
            
            for header in row:
                colname = str(header)
                headers.append(colname)
                data_dict[colname] = []
                
        else:
             colnum = 0
             for datapoint in row:
                data_dict[headers[colnum]].append(int(datapoint))
                colnum += 1
                        
    return data_dict, headers


# Confidence interval bootstrapping function
# Written by: cevans
# URL: https://bitbucket.org/cevans/bootstrap/

from numpy.random import randint
from scipy.stats import norm
from numpy import *

def ci(data, statfun, alpha=0.05, n_samples=10000, method='bca'):
       
        # Ensure that our data is, in fact, an array.
        data = array(data)

        # First create array of bootstrap sample indexes:
        indexes = randint(data.shape[0],size=(n_samples,data.shape[0]))

        # Then apply this to get the bootstrap samples and statistics over them.
        samples = data[indexes]

        stat = array([statfun(x) for x in samples])
        
        # Normal-theory Interval --- doesn't use sorted statistics.
        if method == 'nti':
                bstd = std(stat)
                pass
        
        stat_sorted = sort(stat)
        
        # Percentile Interval
        if method == 'pi':
                return ( stat_sorted[round(n_samples*alpha/2)], stat_sorted[round(n_samples*(1-alpha/2))] )

        # Bias-Corrected Accelerated Interval
        elif method == 'bca':
                ostat = statfun(data)

                z = norm.ppf( ( 1.0*sum(stat < ostat) + 0.5*sum(stat == ostat) ) / (n_samples + 1) )
                
                # Calculate the jackknife distribution and corresponding statistics quickly.
                j_indexes = (lambda n: delete(tile(array(range(0,n)),n),range(0,n*n,n+1)).reshape((n,n-1)))(len(data))
                jstat = [statfun(x) for x in data[j_indexes]]
                jmean = mean(jstat)

                a = sum( (jstat - jmean)**3 ) / ( 6.0 * sum( (jstat - jmean)**2 )**1.5 )

                zp = z + norm.ppf(1-alpha/2)
                zm = z - norm.ppf(1-alpha/2)

                a1 = norm.cdf(z + zm/(1-a*zm))
                a2 = norm.cdf(z + zp/(1-a*zp))

                return (stat_sorted[round(n_samples*a1)],stat_sorted[round(n_samples*a2)])

        else:
                raise "Method %s not supported" % method


# Mangles the CIs from the ci() function into a format that matplotlib's errorbar() function likes
# Input:
#   data        = data to get bootstrapped CIs for
#   statfun     = function to compute CIs over (usually, mean)
#   alpha       = size of CIs (0.05 --> 95% CIs)
#   n_samples   = # of bootstrap populations to construct
# Returns: CI values formatted for matplotlib errorbar() function

def ci_errorbar(data, statfun, alpha=0.05, n_samples=10000, method='bca'):

    data = list(data)
    stat = statfun(data)
    cie = list(ci(data, statfun, alpha, n_samples, method))

    low = stat - float(cie[0])
    high = float(cie[1]) - stat

    return [[low], [high]]
