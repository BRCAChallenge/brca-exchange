

def knapsack(W, wt, val, n):
  """ A Dynamic Programming based Python Program for 0-1 Knapsack problem.
      Returns the maximum value that can be put in a knapsack of capacity W.
  """
    K = [[0 for x in range(W+1)] for x in range(n+1)]
    for i in range(n+1):
        for w in range(W+1):
            if i==0 or w==0:
                K[i][w] = 0
            elif wt[i-1] <= w:
                K[i][w] = max(val[i-1] + K[i-1][w-wt[i-1]],  K[i-1][w])
            else:
                K[i][w] = K[i-1][w]
 
    return K[n][W]

def select_ld(prx_dict, curr_selection):
  """ Greedily select SNPs out of current selection,
      based on proxy dictionary from SNAP output:
      An unselected SNP is selected only if it is not a proxy of any selected SNP.
      Returns the set of selected SNPs.
  """
  selection = set()
  correlated = set()
  for snp in curr_selection:
    if snp in correlated:
      continue
    if snp in prx_dict.keys():
      c = [prx['proxy'] for prx in prx_dict[snp]]
    else:
      c = []
    correlated.update(c)
    selection.add(snp)
  return selection

