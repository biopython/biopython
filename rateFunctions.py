#use the appropriate function of your choice
#k is always positive
#t is always positive as well
#A0 is the initial amount


#the functions all return the expected concentration/amount of A at time t

def firstOrder(k, A0, t): 
	return (A0*(2.71828182845904523536028747135266249775724709369^(-1*t*k)))

def secondOrder(k,A0,t):
	return (1/((1/A0)+k*t))

def zeroOrder(k,A0,t):
	return ((-1*k*t)+A0)




#this function gives reaction rate k for an Arrenheius reaction
#give T in KELVINS, and Ea in JOULES
def Arrenheius(A,Ea,T):
	return (A*2.71828182845904523536028747135266249775724709369^(-1*Ea/(8.314*T)))

def 