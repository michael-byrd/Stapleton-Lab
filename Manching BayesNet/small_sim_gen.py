

from numpy import genfromtxt
my_data = genfromtxt('SimulatedResponse.csv', delimiter=',')


print(my_data[3])
