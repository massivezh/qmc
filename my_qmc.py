from qmc import *

def QF(ON,OFF,cnt):
	DC=range(2**cnt)
	for i in ON:
		DC.remove(i)
	for i in OFF:
		DC.remove(i)
	QmcFunction(ON,DC)

q = QuineMcCluskey([
	QF([3,4,6,7], [1,5],3),
	QF([3,4,6], [0,1,7],3),
	QF([0,2], [1,3,5,7],3)
	], LITERALS_COST_FUNCTION)

q.findPrimeImplicants()
q.printImplicantList()
q.printTable()
print q.simplify()
