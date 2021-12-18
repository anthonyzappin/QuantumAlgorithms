
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit import execute, IBMQ
from qiskit import BasicAer
import sys
import math
import array
import fractions
import numpy as np


def check_if_power(N):
    b=2
    while (2**b) <= N:
        a = 1
        c = N
        while (c-a) >= 2:
            m = int( (a+c)/2 )

            if (m**b) < (N+1):
                p = int( (m**b) )
            else:
                p = int(N+1)

            if int(p) == int(N):
                print('N is {0}^{1}'.format(int(m),int(b)) )
                return True

            if p<N:
                a = int(m)
            else:
                c = int(m)
        b=b+1

    return False


def get_value_a(N):

    ok='0'
    a=2

    while math.gcd(a,N)!=1:
        a=a+1

    smallest_a = a

    ok = input('Is the number {0} ok for a? Press 0 if not, other number if yes: '.format(a))
    if ok=='0':
        if(N==3):
            print('Number {0} is the only one you can use. Using {1} as value for a\n'.format(a,a))
            return a
        a=a+1

    while ok=='0':
    
        while math.gcd(a,N)!=1:
            a=a+1
    
        ok = input('Is the number {0} ok for a? Press 0 if not, other number if yes: '.format(a))

        if ok!='0':
            break
        
        a=a+1

        if a>(N-1):
            print('You rejected all options for value a, selecting the smallest one\n')
            a=smallest_a
            break

    print('Using {0} as value for a\n'.format(a))

    return a

#Function to apply the continued fractions to find r and the gcd
def get_factors(x_value,t_upper,N,a):

    if x_value<=0:
        print('x_value is <= 0, there are no continued fractions\n')
        return False

    T = pow(2,t_upper)

    x_over_T = x_value/T
    
    i=0
    b = array.array('i')
    t = array.array('f')

    b.append(math.floor(x_over_T))
    t.append(x_over_T - b[i])

    while i>=0:
    
        if i>0:
            b.append( math.floor( 1 / (t[i-1]) ) ) 
            t.append( ( 1 / (t[i-1]) ) - b[i] )

        aux = 0
        j=i
        while j>0:    
            aux = 1 / ( b[j] + aux )      
            j = j-1
        
        aux = aux + b[0]

        frac = fractions.Fraction(aux).limit_denominator()
        den=frac.denominator

        print('Approximation number {0} of continued fractions:'.format(i+1))
        print("Numerator:{0} \t\t Denominator: {1}\n".format(frac.numerator,frac.denominator))

        i=i+1

        if (den%2) == 1:
            if i>=15:
                print('Returning because have already done too much tries')
                return False
            print('Odd denominator, will try next iteration of continued fractions\n')
            continue
    
        exponential = 0

        if den<1000:
            exponential=pow(a , (den/2))
        
        if math.isinf(exponential)==1 or exponential>1000000000:
            print('Denominator of continued fraction is too big!\n')
            aux_out = input('Input number 1 if you want to continue searching, other if you do not: ')
            if aux_out != '1':
                return False
            else:
                continue

        putting_plus = int(exponential + 1)
        putting_minus = int(exponential - 1)
        one_factor = math.gcd(putting_plus,N)
        other_factor = math.gcd(putting_minus,N)
    
        if one_factor==1 or one_factor==N or other_factor==1 or other_factor==N:
            print('Found just trivial factors, not good enough\n')

            if t[i-1]==0:
                print('The continued fractions found exactly x_final/(2^(2n)) , leaving funtion\n')
                return False
            if i<15:
                aux_out = input('Input number 1 if you want to continue searching, other if you do not: ')
                if aux_out != '1':
                    return False       
            else:
                print('Returning because have already done too many tries\n')
                return False         
        else:
            print('The factors of {0} are {1} and {2}\n'.format(N,one_factor,other_factor))
            print('Found the desired factors!\n')
            return True

#Functions that calculate the modular inverse using Euclid's algorithm
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)
def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m

#Quantum Fourier Transformation Method
def create_QFT(circuit,up_reg,n,with_swaps):
    i=n-1

    while i>=0:
        circuit.h(up_reg[i])        
        j=i-1  
        while j>=0:
            if (np.pi)/(pow(2,(i-j))) > 0:
                circuit.cu1( (np.pi)/(pow(2,(i-j))) , up_reg[i] , up_reg[j] )
                j=j-1   
        i=i-1  

    if with_swaps==1:
        i=0
        while i < ((n-1)/2):
            circuit.swap(up_reg[i], up_reg[n-1-i])
            i=i+1

#Inverse QFT
def create_inverse_QFT(circuit,up_reg,n,with_swaps):
    if with_swaps==1:
        i=0
        while i < ((n-1)/2):
            circuit.swap(up_reg[i], up_reg[n-1-i])
            i=i+1

    i=0
    while i<n:
        circuit.h(up_reg[i])
        if i != n-1:
            j=i+1
            y=i
            while y>=0:
                 if (np.pi)/(pow(2,(j-y))) > 0:
                    circuit.cu1( - (np.pi)/(pow(2,(j-y))) , up_reg[j] , up_reg[y] )
                    y=y-1   
        i=i+1

def getAngle(a, N):
    s=bin(int(a))[2:].zfill(N) 
    angle = 0
    for i in range(0, N):
        if s[N-1-i] == '1': 
            angle += math.pow(2, -(N-i))
    angle *= np.pi
    return angle

def getAngles(a,N):
    s=bin(int(a))[2:].zfill(N) 
    angles=np.zeros([N])
    for i in range(0, N):
        for j in range(i,N):
            if s[j]=='1':
                angles[N-i-1]+=math.pow(2, -(j-i))
        angles[N-i-1]*=np.pi
    return angles


def ccphase(circuit, angle, ctl1, ctl2, tgt):
    circuit.cu1(angle/2,ctl1,tgt)
    circuit.cx(ctl2,ctl1)
    circuit.cu1(-angle/2,ctl1,tgt)
    circuit.cx(ctl2,ctl1)
    circuit.cu1(angle/2,ctl2,tgt)

def phiADD(circuit, q, a, N, inv):
    angle=getAngles(a,N)
    for i in range(0,N):
        if inv==0:
            circuit.u1(angle[i],q[i])
        else:
            circuit.u1(-angle[i],q[i])

def cphiADD(circuit, q, ctl, a, n, inv):
    angle=getAngles(a,n)
    for i in range(0,n):
        if inv==0:
            circuit.cu1(angle[i],ctl,q[i])
        else:
            circuit.cu1(-angle[i],ctl,q[i])
        

def ccphiADD(circuit,q,ctl1,ctl2,a,n,inv):
    angle=getAngles(a,n)
    for i in range(0,n):
        if inv==0:
            ccphase(circuit,angle[i],ctl1,ctl2,q[i])
        else:
            ccphase(circuit,-angle[i],ctl1,ctl2,q[i])
        
def ccphiADDmodN(circuit, q, ctl1, ctl2, aux, a, N, n):
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 0)
    phiADD(circuit, q, N, n, 1)
    create_inverse_QFT(circuit, q, n, 0)
    circuit.cx(q[n-1],aux)
    create_QFT(circuit,q,n,0)
    cphiADD(circuit, q, aux, N, n, 0)
    
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 1)
    create_inverse_QFT(circuit, q, n, 0)
    circuit.x(q[n-1])
    circuit.cx(q[n-1], aux)
    circuit.x(q[n-1])
    create_QFT(circuit,q,n,0)
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 0)

def ccphiADDmodN_inv(circuit, q, ctl1, ctl2, aux, a, N, n):
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 1)
    create_inverse_QFT(circuit, q, n, 0)
    circuit.x(q[n-1])
    circuit.cx(q[n-1],aux)
    circuit.x(q[n-1])
    create_QFT(circuit, q, n, 0)
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 0)
    cphiADD(circuit, q, aux, N, n, 1)
    create_inverse_QFT(circuit, q, n, 0)
    circuit.cx(q[n-1], aux)
    create_QFT(circuit, q, n, 0)
    phiADD(circuit, q, N, n, 0)
    ccphiADD(circuit, q, ctl1, ctl2, a, n, 1)

def cMULTmodN(circuit, ctl, q, aux, a, N, n):
    create_QFT(circuit,aux,n+1,0)
    for i in range(0, n):
        ccphiADDmodN(circuit, aux, q[i], ctl, aux[n+1], (2**i)*a % N, N, n+1)
    create_inverse_QFT(circuit, aux, n+1, 0)

    for i in range(0, n):
        circuit.cswap(ctl,q[i],aux[i])

    a_inv = modinv(a, N)
    create_QFT(circuit, aux, n+1, 0)
    i = n-1
    while i >= 0:
        ccphiADDmodN_inv(circuit, aux, q[i], ctl, aux[n+1], math.pow(2,i)*a_inv % N, N, n+1)
        i -= 1
    create_inverse_QFT(circuit, aux, n+1, 0)

#Main Program
if __name__ == '__main__':

    N = int(input('Please insert integer number N: '))

    print('input number was: {0}\n'.format(N))

    if N==1 or N==0: 
       print('Please put an N different from 0 and from 1')
       exit()
    

    if (N%2)==0:
        print('N is even, so does not make sense!')
        exit()
    
    if check_if_power(N)==True:
       exit()

    print('Not an easy case, using the quantum circuit is necessary\n')

    a = get_value_a(N)

    n = math.ceil(math.log(N,2))
    
    print('Total number of qubits used: {0}\n'.format(2*n+3))
    aux = QuantumRegister(n+2)
    up_reg = QuantumRegister(1)
    down_reg = QuantumRegister(n)
    up_classic = ClassicalRegister(2*n)
    c_aux = ClassicalRegister(1)
    circuit = QuantumCircuit(down_reg , up_reg , aux, up_classic, c_aux)
    circuit.x(down_reg[0])

    for i in range(0, 2*n):
        circuit.x(up_reg).c_if(c_aux, 1)
        circuit.h(up_reg)
        cMULTmodN(circuit, up_reg[0], down_reg, aux, a**(2**(2*n-1-i)), N, n)
        for j in range(0, 2**i):
            circuit.u1(getAngle(j, i), up_reg[0]).c_if(up_classic, j)
        circuit.h(up_reg)
        circuit.measure(up_reg[0], up_classic[i])
        circuit.measure(up_reg[0], c_aux[0])

    number_shots=int(input('Number of times to run the circuit: '))
    if number_shots < 1:
        print('Please run the circuit at least one time...')
        exit()

    print('Executing the circuit {0} times for N={1} and a={2}\n'.format(number_shots,N,a))

    simulation = execute(circuit, backend=BasicAer.get_backend('qasm_simulator'),shots=number_shots)

    sim_result=simulation.result()
    counts_result = sim_result.get_counts(circuit)

    print('Printing the various results followed by how many times they happened (out of the {} cases):\n'.format(number_shots))
    i=0
    while i < len(counts_result):
        print('Result \"{0}\" happened {1} times out of {2}'.format(list(sim_result.get_counts().keys())[i],list(sim_result.get_counts().values())[i],number_shots))
        i=i+1
    
    print(' ')
    prob_success=0

    i=0
    while i < len(counts_result):

        all_registers_output = list(sim_result.get_counts().keys())[i]
        output_desired = all_registers_output.split(" ")[1]
        x_value = int(output_desired, 2)
        prob_this_result = 100 * ( int( list(sim_result.get_counts().values())[i] ) ) / (number_shots)

        print("------> Analysing result {0}. This result happened in {1:.4f} % of all cases\n".format(output_desired,prob_this_result))

        print('In decimal, x_final value for this result is: {0}\n'.format(x_value))
 
        success = get_factors(int(x_value),int(2*n),int(N),int(a))

        if success==True:
            prob_success = prob_success + prob_this_result

        i=i+1

    print("\nUsing a={0}, found the factors of N={1} in {2:.4f} % of the cases\n".format(a,N,prob_success))
