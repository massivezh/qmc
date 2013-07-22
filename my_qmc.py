from qmc import *


q = QuineMcCluskey([
	QmcFunction([3,4,6,7], [1,5],3),
	QmcFunction([3,4,6], [0,1,7],3),
	QmcFunction([0,2], [1,3,5,7],3)
	], LITERALS_COST_FUNCTION)

q.findPrimeImplicants()
q.printImplicantList()
q.printTable()
print q.simplify()
