import numpy as np
import sys

filename = sys.argv[1]
columns = [3,4]

column_sum = [0 for _ in columns]

with open(filename, 'r') as file:
  for line in file:
    if line[0] == '#':
      # ignore comment lines
      pass
    else:
      values = line.split(' ')
      # get values from columns of interest
      for col in range(len(columns)):
        column_sum[col] += eval(values[columns[col]])

print("column sums:")
for col in range(len(columns)):
    print("column ", columns[col], ": ", column_sum[col])
