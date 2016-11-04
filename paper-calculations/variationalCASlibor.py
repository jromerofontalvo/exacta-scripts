"""This module explicitly solves for the ghetto CC variational ansatz.
"""
import time
import numpy
import scipy
import warnings
import commutators
import scipy.sparse
import scipy.optimize
import resource
import scipy.sparse.linalg
from itertools import combinations
from sys import argv
import os

# Return matrix representation of creation or annihilation operator.

def JordanWignerTerm(index, n_orbitals):
  """Make a matrix representation of a fermion operator.

  Args:
    index: This is a nonzero integer. The integer indicates the tensor
      factor and the sign indicates raising or lowering.
    n_orbitals: This int gives the number of orbitals.

  Returns: The corresponding fermion operator.
  """
  # Define operators.
  I = scipy.sparse.csr_matrix([[1, 0], [0, 1]], dtype=float)
  Z = scipy.sparse.csr_matrix([[1, 0], [0, -1]], dtype=float)
  Q_raise = scipy.sparse.csr_matrix([[0, 0], [1, 0]], dtype=float)
  Q_lower = scipy.sparse.csr_matrix([[0, 1], [0, 0]], dtype=float)

  # Construct fermionic operator.
  assert index
  orbital = abs(index)
  operator = 1
  for tensor_factor in range(orbital - 1):
    operator = scipy.sparse.kron(operator, I, 'csr')
  if index > 0:
    operator = scipy.sparse.kron(operator, Q_raise, 'csr')
  else:
    operator = scipy.sparse.kron(operator, Q_lower, 'csr')
  for tensor_factor in range(n_orbitals - orbital):
    operator = scipy.sparse.kron(operator, Z, 'csr')
  assert scipy.sparse.isspmatrix_csr(operator)
  return operator

def JordanWignerTerm2(index, n_orbitals):
  """Make a matrix representation of a fermion operator.

  Args:
    index: This is a nonzero integer. The integer indicates the tensor
      factor and the sign indicates raising or lowering.
    n_orbitals: This int gives the number of orbitals.

  Returns: The corresponding fermion operator.
  """
  # Define operators.
  I = scipy.sparse.csr_matrix([[1, 0], [0, 1]], dtype=float)
  Z = scipy.sparse.csr_matrix([[1, 0], [0, -1]], dtype=float)
  Q_raise = scipy.sparse.csr_matrix([[0, 0], [1, 0]], dtype=float)
  Q_lower = scipy.sparse.csr_matrix([[0, 1], [0, 0]], dtype=float)

  # Construct fermionic operator.
  assert index
  orbital = abs(index)
  operator = 1
  for tensor_factor in range(orbital - 1):
    operator = scipy.sparse.kron(operator, I, 'csr')
  if index > 0:
    operator = scipy.sparse.kron(operator, Q_raise, 'csr')
  else:
    operator = scipy.sparse.kron(operator, Q_lower, 'csr')
  for tensor_factor in range(n_orbitals - orbital):
    operator = scipy.sparse.kron(operator, Z, 'csr')
  assert scipy.sparse.isspmatrix_csr(operator)
  return operator

def JordanWignerTermSym(index, n_orbitals):
  """Make a matrix representation of a fermion operator.

  Args:
    index: This is a nonzero integer. The integer indicates the tensor
      factor and the sign indicates raising or lowering.
    n_orbitals: This int gives the number of orbitals.

  Returns: The corresponding fermion operator.
  """

  # Construct fermionic operator.
  assert index
  orbital = abs(index)
  op=[]
  if index > 0:
    op += 'Q'
    op += str(orbital)
    op += '-'
  else:
    op += 'q'
    op += str(orbital)
    op += '-'
  for tensor_factor in range(n_orbitals - orbital):
    op += 'Z'
    op += repr(int(tensor_factor))
    op += '-'
  form = ''.join(op)
  return form

# Initialize Jordan-Wigner dictionary so terms aren't constantly reformed.
def GetJordanWignerTerms(n_orbitals):
    jw_terms = {}
    for index in range(-n_orbitals, n_orbitals + 1):
      if index:
        jw_terms[index] = JordanWignerTerm(index, n_orbitals)
    return jw_terms

def GetJordanWignerTermsSym(n_orbitals):
    jw_terms_sym = {}
    for index in range(-n_orbitals, n_orbitals + 1):
      if index:
        jw_terms_sym[index] = JordanWignerTermSym(index, n_orbitals)
    return jw_terms

# Return a particular operator given the term and coefficient.
def MatrixForm(coefficient, term, jw_terms,
    add_conjugates=False, anti_hermitian=False):
  operator = coefficient
  for index in term:
    operator = operator * jw_terms[index]
  if add_conjugates and commutators.GetConjugate(term):
    if anti_hermitian:
      #print 'Enter to T-T^dagger'
      operator = 1j * (operator - operator.getH())
    else:
      #print 'Enter to i(T+T^dagger')
      operator = operator + operator.getH()
  return operator

# Return a particular operator given the term and coefficient.
def MatrixFormNoZ(coefficient, term, n_orbitals, pauli, t1, t2):
  # Construct fermionic operator.
  operator1 = coefficient/2.0
  operator2 = -coefficient/2.0
  counter = -1
  for factor in range(1,n_orbitals+1):
    if factor in [abs(x) for x in term]:
      print("inside matrixFormNoZ:", term, factor)
      counter += 1
      operator1 = scipy.sparse.kron(operator1, pauli[t1[counter]], 'csr')
      operator2 = scipy.sparse.kron(operator2, pauli[t2[counter]], 'csr')
    else:    
      operator1 = scipy.sparse.kron(operator1, pauli['I'], 'csr')
      operator2 = scipy.sparse.kron(operator2, pauli['I'], 'csr')
  operator = operator1 + operator2
  assert scipy.sparse.isspmatrix_csr(operator)
  return operator


# Expectation.
def Expectation(operator, psi):
  psi = scipy.sparse.csr_matrix(psi)
  operator = scipy.sparse.csr_matrix(operator)
  expectation = psi.getH() * operator * psi
  assert expectation.get_shape() == (1, 1)
  return numpy.real(expectation[0, 0])

# Expectation.
def Overlap(psi1, psi2):
  psi1 = scipy.sparse.csr_matrix(psi1)
  psi2 = scipy.sparse.csr_matrix(psi2)
  overlap = psi1.getH() * psi2
  assert overlap.get_shape() == (1, 1)
  return numpy.real(overlap[0, 0])

# print state
def printState(state):
  print(state.shape)
  for n in range(0,state.shape[0]):
    if abs(state[n,0])>1E-05:
      print(n,state[n,0])

# read vector as state from MOLPRO file
def readMolproState(molproFile):
  state = numpy.loadtxt(molproFile,dtype=complex,skiprows=1)
  state = scipy.sparse.csr_matrix(state[::-1])
  print("within readmolpro:",state.shape)
  state = state.transpose(copy=True)
  print("within readmolpro:",state.shape)
  return state

# discard Term
def keepTerm(term, n_frozen, max_virtual):
  # verifies if excitation operator is within the 
  # active space
  value = True
  for a in term:
    print a
    if abs(a) <= n_frozen or abs(a)>max_virtual:
      value = False
      break 
  return value

# Get number operator.
def NumberOperator(n_orbitals):
  number_operator = 0
  jw_terms = GetJordanWignerTerms(n_orbitals)
  for tensor_factor in range(1, n_orbitals + 1):
    term = [tensor_factor, -tensor_factor]
    operator = MatrixForm(1, term, jw_terms)
    number_operator = number_operator + operator
  return number_operator


# Initialize Hartree-Fock state.
def HartreeFockState(n_electrons, n_orbitals):
  occupied = scipy.sparse.csr_matrix([[0], [1]], dtype=float)
  psi = 1.
  unoccupied = scipy.sparse.csr_matrix([[1], [0]], dtype=float)
  for orbital in range(n_electrons):
    psi = scipy.sparse.kron(psi, occupied, 'csr')
  for orbital in range(n_orbitals - n_electrons):
    psi = scipy.sparse.kron(psi, unoccupied, 'csr')
  return psi

# Initialize Hartree-Fock state.
def HartreeFockState2(state):
  occupied = scipy.sparse.csr_matrix([[0], [1]], dtype=float)
  psi = 1.
  unoccupied = scipy.sparse.csr_matrix([[1], [0]], dtype=float)
  for x in state:
    print(x)
    if x=='1':
      psi = scipy.sparse.kron(psi, occupied, 'csr')
    elif x=='0':
      psi = scipy.sparse.kron(psi, unoccupied, 'csr')
  return psi

# Test if matrix is Hermitian.
def IsHermitian(matrix):
  conjugate = matrix.getH()
  difference = matrix - conjugate
  if difference.nnz:
    discrepancy = max(map(abs, difference.data))
    if discrepancy > 1e-9:
      print 'Hermitian discrepancy = %s.' % repr(discrepancy)
      return False
  return True


# Construct an operator from terms and coefficients.
def MakeOperator(coefficients, terms, verbose=False):

  # Initialize.
  n_terms = len(coefficients)
  one_percent = numpy.rint(numpy.ceil(n_terms / 100.))
  start = time.clock()
  n_orbitals = commutators.OrbitalCount(terms)
  jw_terms = GetJordanWignerTerms(n_orbitals)

  # Loop over terms.
  operator = 0
  for i, (coefficient, term) in enumerate(zip(coefficients, terms)):
    operator = operator + MatrixForm(coefficient, term, jw_terms)

    # Report progress.
    if verbose and not (i + 1) % one_percent:
      percent_complete = numpy.rint(100. * (i + 1) / n_terms)
      elapsed = time.clock() - start
      rate = elapsed / percent_complete
      eta = rate * (100 - percent_complete)
      print('%s. Computation %i%% complete. Approximately %i '
            'minute(s) remaining.' % (time.strftime(
                '%B %d at %H:%M:%S', time.localtime()),
            percent_complete, round(eta / 60)))

  assert IsHermitian(operator)
  return operator


# Save a sparse matrix.
def SaveSparse(name, operator):
  numpy.savez(name, data=operator.data, indices=operator.indices,
              indptr=operator.indptr, shape=operator.shape)


# Load a sparse matrix.
def LoadSparse(name):
  loader = numpy.load(name)
  operator = scipy.sparse.csr_matrix((loader['data'], loader['indices'],
                                      loader['indptr']), shape=loader['shape'])
  return operator


# Get the relevant Slater determinants.
def GetDeterminants(n_orbitals, n_electrons, n_excitations):

  # Initialize reference state.
  reference = numpy.zeros(n_orbitals, int)
  reference[:n_electrons] = 1
  states = []

  # Loop over excitation level.
  for m in range(n_excitations + 1):
    for occupied in combinations(range(n_electrons), r=m):
      for unoccupied in combinations(range(n_electrons, n_orbitals), r=m):
        state = numpy.copy(reference)
        state[list(unoccupied)] = 1
        state[list(occupied)] = 0
        states += [state]
  return numpy.array(states)


# Return a projector into a restricted Fock space.
def ConfigurationProjector(n_orbitals, n_electrons, n_excitations):

  # Name.
  if n_excitations == 'all':
    n_excitations = n_orbitals - n_electrons
  name = 'data/operators/projector_%i_%i_%i.npz'\
      % (n_orbitals, n_electrons, n_excitations)

  # Attempt to load.
  try:
    projector = LoadSparse(name)

  except:
    # Initialize projector computation.
    n_hilbert = 2 ** n_orbitals
    states = GetDeterminants(n_orbitals, n_electrons, n_excitations)
    unoccupied = scipy.sparse.csr_matrix([[1], [0]], dtype=int)
    occupied = scipy.sparse.csr_matrix([[0], [1]], dtype=int)
    number_operator = NumberOperator(n_orbitals)

    # Construct projector.
    projector = 0
    for state in states:
      
      # Construct computational basis state in Hilbert space.
      ket = 1
      for qubit in state:
        if qubit:
          ket = scipy.sparse.kron(ket, occupied, 'csr')
        else:
          ket = scipy.sparse.kron(ket, unoccupied, 'csr')

      # Add to projector.
      density = ket * ket.getH()
      projector = projector + density

    # Save and return projector.
    SaveSparse(name, projector)
  return projector


# Return an operator.
def GetHamiltonian2(molecule, basis, n_excitations='FCI'):
  recompute = 1
  # Try to load it.
  name = 'data/operators/%s_%s_hamiltonian_%s.npz'\
      % (molecule, basis, str(n_excitations))
  try:
    assert not recompute
    hamiltonian = LoadSparse(name)

  # Compute pre-projected Hamiltonian.
  except:
    repulsion, coefficients, terms = commutators.GetHamiltonianTerms(molecule, basis, add_conjugates=True)
    hamiltonian = MakeOperator(coefficients, terms)

    # Project, if necessary.
    if n_excitations != 'FCI':
      n_hilbert = hamiltonian.shape[0]
      n_electrons = commutators.ElectronCount(molecule)
      n_orbitals = int(numpy.rint(numpy.log2(n_hilbert)))
      if n_orbitals - n_electrons > n_excitations:
        projector = ConfigurationProjector(
            n_orbitals, n_electrons, n_excitations)
        hamiltonian = projector * hamiltonian * projector.getH()

    # Save and return.
    SaveSparse(name, hamiltonian)
  return hamiltonian


# Return an operator.
def GetHamiltonian(molecule, basis, n_excitations='FCI', n_orbitals=None, path=None):
  # Try to load it.
  name = 'data/operators/%s_%s_hamiltonian_%s.npz'\
      % (molecule, basis, str(n_excitations))
  # Compute pre-projected Hamiltonian.
  print "Inside GetHamiltonian: path", path
  fci_energy, repulsion, coefficients, terms = commutators.GetHamiltonianTerms(molecule, basis, add_conjugates=True, n_orbitals=n_orbitals, path=path)
  print terms
  hamiltonian = MakeOperator(coefficients, terms)
  # n_orbitals = commutators.OrbitalCount(terms)
  print "orbitals inside getHamiltonian:", n_orbitals
  # Project, if necessary.
  if n_excitations != 'FCI':
    n_hilbert = hamiltonian.shape[0]
    n_electrons = commutators.ElectronCount(molecule)
    n_orbitals = int(numpy.rint(numpy.log2(n_hilbert)))
    if n_orbitals - n_electrons > n_excitations:
      projector = ConfigurationProjector(
        n_orbitals, n_electrons, n_excitations)
      hamiltonian = projector * hamiltonian * projector.getH()
      
    # Save and return.
  print "orbitals inside getHamiltonian:", n_orbitals
  SaveSparse(name, hamiltonian)
  return hamiltonian, n_orbitals


# Compute and save information about lowest or highest operator eigenvalue.
def SparseDiagonalize(operator, which='SA'):
  max_steps = 1e10
  values, vectors = scipy.sparse.linalg.eigsh(
      operator, 1, which=which, maxiter=max_steps)
  eigenstate = vectors[:, 0]
  print "the minimum energy in the energy vector is:",min(values),values[0]
  eigenstate = scipy.sparse.csr_matrix(eigenstate)
  return scipy.sparse.csr_matrix.transpose(eigenstate)
  
# Return the unitary corresponding to evolution under hamiltonian.
def GetUnitary(time, term):
  exponent = scipy.sparse.csc_matrix(-1j * time * term)
  assert IsHermitian(1j * exponent)
  with warnings.catch_warnings(record=False):
    warnings.simplefilter('ignore')
  unitary = scipy.sparse.linalg.expm(exponent)
  unitary.eliminate_zeros()
  return unitary.tocsr()


# Make a function object to return energy of variational ansatz.
class Objective:

  # Initialize function object.
  def __init__(self, terms, reference_state, hamiltonian, cc_form=False, n_orbitals = None, state_prep_method='exact', nTrotterSteps=1):

    # Initialize dictionary of unitaries.
    self.terms = terms
    self.hamiltonian = hamiltonian
    self.reference_state = reference_state
    if n_orbitals==None:
      n_orbitals = commutators.OrbitalCount(terms)
    jw_terms = GetJordanWignerTerms(n_orbitals)
    self.matrices = {}
    I = scipy.sparse.csr_matrix([[1, 0], [0, 1]], dtype=complex)
    Z = scipy.sparse.csr_matrix([[1, 0], [0, -1]], dtype=complex)
    X = scipy.sparse.csr_matrix([[0, 1], [1, 0]], dtype=complex)
    Y = scipy.sparse.csr_matrix([[0, -1j], [1j, 0]], dtype=complex)
    pauli={'I':I,'X':X,'Y':Y,'Z':Z}
    for term in terms:
      if len(term)==2:
        t1='YX'
        t2='XY'
      if len(term)==4:
        t1='XXXY'
        t2='YYYX'
      self.matrices[tuple(term)] = MatrixForm(
         1., term, jw_terms, add_conjugates=True, anti_hermitian=cc_form)
      # self.matrices[tuple(term)] = MatrixFormNoZ(1., term, n_orbitals, pauli, t1, t2)
      assert IsHermitian(self.matrices[tuple(term)])
    self.min_energy = numpy.inf
    self.calls = 0
    self.state_prep_method = state_prep_method
    self.trotterSteps = nTrotterSteps

  # Method to return unitary which prepares a state.
  def GetStatePrep(self, coefficients):
    unitary = 1.
    total_unitary = 1.
    counter=0
    for coefficient, term in zip(coefficients, self.terms):
      counter=counter+1
      # print 'TERM %i' % (counter)
      # print '%s %s' % (term,repr(float(coefficient)))
      #sterm=SymbolicJordanWigner2(term,coefficient)
      #for i in range(len(sterm)):
      #  print sterm[i]
      term_matrix = self.matrices[tuple(term)].copy()
      term_unit = GetUnitary(coefficient/float(self.trotterSteps), term_matrix)
      # print 'Unitary'
      # print term_unit
      unitary =  term_unit * unitary
    for n in range(0,self.trotterSteps):
      total_unitary = total_unitary * unitary
    return total_unitary

  # Method to return unitary which prepares a state using single unitary.
  def ModStatePrep(self, coefficients):
    hamiltonian = 0.
    for coefficient, term in zip(coefficients, self.terms):
      # if abs(float(coefficient)) > 1E-06:
      #   print '%s %s' % (term,repr(float(coefficient)))
      hamiltonian = hamiltonian + coefficient * self.matrices[tuple(term)].copy()
    unitary = GetUnitary(1., hamiltonian)
    return unitary

  # Write paren method which returns energy given variables.
  def returnOverlap(self, variables, state2):
    state_prep = self.ModStatePrep(variables)
    #state_prep = self.GetStatePrep(variables)
    state = state_prep * self.reference_state
    value = Overlap(state, state2)
    return value

  # Write paren method which returns energy given variables.
  def __call__(self, variables):
    if self.state_prep_method == 'exact':
      state_prep = self.ModStatePrep(variables)
    elif self.state_prep_method == 'trotter':
      state_prep = self.GetStatePrep(variables)
    state = state_prep * self.reference_state
    energy = Expectation(self.hamiltonian, state)
    self.calls += 1
    print 'CALL NUMBER: %i' %(self.calls)
    if energy < self.min_energy:
      print "current energy:", energy
      print 'Energy difference: %e' % (self.min_energy-energy)
      self.min_energy = energy
#      print 'Lowest energy after %i calls is %f.' % (self.calls, energy)
    print("Memory usage:",resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)
    return energy


# Variationally determine optimal coefficients for state prep.
def VariationallyMinimize(terms, reference_state, initial_guess, hamiltonian, cc_form=False, n_orbitals=None, optimization_method=None, state_prep_method='exact', nTrotterSteps=1):
  ObjectiveFunction = Objective(terms, reference_state, hamiltonian, cc_form, n_orbitals, state_prep_method, nTrotterSteps)
  if optimization_method == None or optimization_method == 'LBFGS':
    sol = scipy.optimize.fmin_l_bfgs_b(ObjectiveFunction, initial_guess, fprime=None, approx_grad=1, pgtol=1e-05, epsilon=1e-08, iprint=-1, maxfun=20000, maxiter=20000)
    print sol
    solution = sol[0]
  elif optimization_method == 'COBYLA':
    method = 'COBYLA'
    solution = scipy.optimize.minimize(
      ObjectiveFunction, initial_guess, args=(), method=method, constraints=(), tol=None, callback=None, options={'iprint': 1, 'disp': True, 'maxiter': 20000, 'tol': 1e-04, 'rhobeg': 0.005}).x
  elif optimization_method == 'POWELL':
    method = 'Powell'
    solution = scipy.optimize.minimize(
      ObjectiveFunction, initial_guess, method=method, options={'ftol':1e-05,'xtol':1e-04,'maxfev':20000,'maxiter':20000,}).x
  elif optimization_method == 'ND':
    method = 'Nelder-Mead'
    solution = scipy.optimize.minimize(
      ObjectiveFunction, initial_guess, method=method, options={'ftol':1e-05,'xtol':1e-04,'maxiter':20000,'maxfev':20000}).x
  energy = ObjectiveFunction(solution)
  ncalls = ObjectiveFunction.calls
  return ncalls, solution, energy, ObjectiveFunction

def SymbolicJordanWigner(term,coefficient):
  termlen=len(term)
  if termlen==2:
    if abs(term[1])>abs(term[0]):
      p = abs(term[1])
      q = abs(term[0])
    else:
      p = abs(term[0])
      q = abs(term[1])
    sterm=[]
    coeff=0.5*coefficient
    for i in range(q+1,p):
      sterm += 'Z'
      sterm += str(i)
    sterm += ' '
    sterm += str(coeff)
    sterm += '('
    sterm += 'X'
    sterm += str(q)
    sterm += 'X'
    sterm += str(p)
    sterm += '+'
    sterm += 'Y'
    sterm += str(q)
    sterm += 'Y'
    sterm += str(p)
    sterm += ')'
  elif termlen==4:
    if abs(term[0])==abs(term[2]): 
      q=abs(term[0])
      p=abs(term[1])
      r=abs(term[3])
      sign=-1
      opt=1
    elif abs(term[0])==abs(term[3]):
      q=abs(term[0])
      p=abs(term[1])
      r=abs(term[2])
      sign=1
      opt=1
    elif abs(term[1])==abs(term[2]): 
      q=abs(term[1])
      p=abs(term[0])
      r=abs(term[3])
      sign=1
      opt=1
    elif abs(term[1])==abs(term[3]):
      q=abs(term[1])
      p=abs(term[0])
      r=abs(term[2])
      sign=-1
      opt=1
    else:
      opt=2
      sign=1
      if (abs(term[0])>abs(term[2]) and abs(term[1])>abs(term[3])) or (abs(term[0])>abs(term[3]) and abs(term[1])>abs(term[2])):
        if abs(term[1])>abs(term[0]):
          p=abs(term[1])
          q=abs(term[0])
          sign=-1*sign
        else:
          p=abs(term[0])
          q=abs(term[1])
          sign=1*sign
        if abs(term[3])>abs(term[2]):
          r=abs(term[3])
          s=abs(term[2])
          sign=-1*sign
        else:
          r=abs(term[2])
          s=abs(term[3])
          sign=1*sign
      else:
        if abs(term[1])>abs(term[0]):
          r=abs(term[1])
          s=abs(term[0])
          sign=-1*sign
        else:
          r=abs(term[0])
          s=abs(term[1])
          sign=1*sign
        if abs(term[3])>abs(term[2]):
          p=abs(term[3])
          q=abs(term[2])
          sign=-1*sign
        else:
          p=abs(term[2])
          q=abs(term[3])
          sign=1*sign
    if opt==1:
      if r>p:
        aux=r
        r=p
        p=aux
      sterm=[]
      coeff=sign*coefficient
      for i in range(r+1,p):
        sterm += 'Z'
        sterm += str(i)
      sterm += '['
      sterm += str(coeff)
      sterm += '('
      sterm += 'X'
      sterm += str(r)
      sterm += 'X'
      sterm += str(p)
      sterm += '+'
      sterm += 'Y'
      sterm += str(r)
      sterm += 'Y'
      sterm += str(p)
      sterm += ')'
      sterm += '-'
      sterm += 'Z'
      sterm += str(q)
      sterm += ' '
      sterm += str(coeff)
      sterm += '('
      sterm += 'X'
      sterm += str(r)
      sterm += 'X'
      sterm += str(p)
      sterm += '+'
      sterm += 'Y'
      sterm += str(r)
      sterm += 'Y'
      sterm += str(p)
      sterm += ')'
      sterm += ']'
      coeff=coefficient
    elif opt==2:
      sterm=[]
      coeff=sign*coefficient/8.0
      for i in range(s+1,r):
        sterm += 'Z'
        sterm += str(i)
      for i in range(q+1,p):
        sterm += 'Z'
        sterm += str(i)
      sterm += ' '
      sterm += str(coeff)
      sterm += '('
      sterm += 'X'
      sterm += str(s)
      sterm += 'X'
      sterm += str(r)
      sterm += 'X'
      sterm += str(q)
      sterm += 'X'
      sterm += str(p)
      sterm += '-'
      sterm += 'X'
      sterm += str(s)
      sterm += 'X'
      sterm += str(r)
      sterm += 'Y'
      sterm += str(q)
      sterm += 'Y'
      sterm += str(p)
      sterm += '+'
      sterm += 'X'
      sterm += str(s)
      sterm += 'Y'
      sterm += str(r)
      sterm += 'X'
      sterm += str(q)
      sterm += 'Y'
      sterm += str(p)
      sterm += '+'
      sterm += 'Y'
      sterm += str(s)
      sterm += 'X'
      sterm += str(r)
      sterm += 'X'
      sterm += str(q)
      sterm += 'Y'
      sterm += str(p)
      sterm += '+'
      sterm += 'Y'
      sterm += str(s)
      sterm += 'X'
      sterm += str(r)
      sterm += 'Y'
      sterm += str(q)
      sterm += 'X'
      sterm += str(p)
      sterm += '-'
      sterm += 'Y'
      sterm += str(s)
      sterm += 'Y'
      sterm += str(r)
      sterm += 'X'
      sterm += str(q)
      sterm += 'X'
      sterm += str(p)
      sterm += '+'
      sterm += 'X'
      sterm += str(s)
      sterm += 'Y'
      sterm += str(r)
      sterm += 'Y'
      sterm += str(q)
      sterm += 'X'
      sterm += str(p)
      sterm += '+'
      sterm += 'Y'
      sterm += str(s)
      sterm += 'Y'
      sterm += str(r)
      sterm += 'Y'
      sterm += str(q)
      sterm += 'Y'
      sterm += str(p)
      sterm += ')'
  symterm=''.join(sterm)
  return symterm

def SymbolicJordanWigner2(term,coefficient):
  termlen=len(term)
  if termlen==2:
    if abs(term[1])>abs(term[0]):
      p = abs(term[1])
      q = abs(term[0])
    else:
      p = abs(term[0])
      q = abs(term[1])
    term1=[]
    term2=[]
    aux=range(2)
    coeff=''.join(str(0.5*coefficient))
    # term 1
    for i in range(q+1,p):
      aux[0] = 'Z'
      aux[1] = str(i)
      part = ''.join(aux)
      term1.append(part) 
    aux[0] = 'X'
    aux[1] = str(q)
    part = ''.join(aux)
    term1.append(part)
    aux[0] = 'X'
    aux[1] = str(p)
    part = ''.join(aux)
    term1.append(part)
    term1.append(str(coeff))
    # term 2
    for i in range(q+1,p):
      aux[0] = 'Z'
      aux[1] = str(i)
      part = ''.join(aux)
      term2.append(part) 
    aux[0] = 'Y'
    aux[1] = str(q)
    part = ''.join(aux)
    term2.append(part)
    aux[0] = 'Y'
    aux[1] = str(p)
    part = ''.join(aux)
    term2.append(part)
    term2.append(str(coeff))
    sterm=range(2)
    term1=' '.join(term1)
    term2=' '.join(term2)
    sterm[0]=term1
    sterm[1]=term2
  elif termlen==4:
    if abs(term[0])==abs(term[2]): 
      q=abs(term[0])
      p=abs(term[1])
      r=abs(term[3])
      sign=-1
      opt=1
    elif abs(term[0])==abs(term[3]):
      q=abs(term[0])
      p=abs(term[1])
      r=abs(term[2])
      sign=1
      opt=1
    elif abs(term[1])==abs(term[2]): 
      q=abs(term[1])
      p=abs(term[0])
      r=abs(term[3])
      sign=1
      opt=1
    elif abs(term[1])==abs(term[3]):
      q=abs(term[1])
      p=abs(term[0])
      r=abs(term[2])
      sign=-1
      opt=1
    else:
      opt=2
      sign=1
      if (abs(term[0])>abs(term[2]) and abs(term[1])>abs(term[3])) or (abs(term[0])>abs(term[3]) and abs(term[1])>abs(term[2])):
        if abs(term[1])>abs(term[0]):
          p=abs(term[1])
          q=abs(term[0])
          sign=-1*sign
        else:
          p=abs(term[0])
          q=abs(term[1])
          sign=1*sign
        if abs(term[3])>abs(term[2]):
          r=abs(term[3])
          s=abs(term[2])
          sign=-1*sign
        else:
          r=abs(term[2])
          s=abs(term[3])
          sign=1*sign
      else:
        if abs(term[1])>abs(term[0]):
          r=abs(term[1])
          s=abs(term[0])
          sign=-1*sign
        else:
          r=abs(term[0])
          s=abs(term[1])
          sign=1*sign
        if abs(term[3])>abs(term[2]):
          p=abs(term[3])
          q=abs(term[2])
          sign=-1*sign
        else:
          p=abs(term[2])
          q=abs(term[3])
          sign=1*sign
    if opt==1:
      if r>p:
        aux=r
        r=p
        p=aux
      term1=[]
      term2=[]
      term3=[]
      term4=[]
      aux=range(2)
      coeff=''.join(str(sign*coefficient))
      # term 1
      for i in range(r+1,p):
        aux[0] = 'Z'
        aux[1] = str(i)
        part = ''.join(aux)
        term1.append(part) 
      aux[0] = 'X'
      aux[1] = str(r)
      part = ''.join(aux)
      term1.append(part)
      aux[0] = 'X'
      aux[1] = str(p)
      part = ''.join(aux)
      term1.append(part)
      term1.append(str(coeff))
      # term 2
      for i in range(r+1,p):
        aux[0] = 'Z'
        aux[1] = str(i)
        part = ''.join(aux)
        term2.append(part) 
      aux[0] = 'Y'
      aux[1] = str(r)
      part = ''.join(aux)
      term2.append(part)
      aux[0] = 'Y'
      aux[1] = str(p)
      part = ''.join(aux)
      term2.append(part)
      term2.append(str(coeff))
      # term 3
      coeff=''.join(str(-sign*coefficient))
      for i in range(r+1,p):
        aux[0] = 'Z'
        aux[1] = str(i)
        part = ''.join(aux)
        term3.append(part) 
      aux[0] = 'Z'
      aux[1] = str(q)
      part = ''.join(aux)
      term3.append(part)
      aux[0] = 'X'
      aux[1] = str(r)
      part = ''.join(aux)
      term3.append(part)
      aux[0] = 'X'
      aux[1] = str(p)
      part = ''.join(aux)
      term3.append(part)
      term3.append(str(coeff))
      # term 4
      for i in range(r+1,p):
        aux[0] = 'Z'
        aux[1] = str(i)
        part = ''.join(aux)
        term4.append(part) 
      aux[0] = 'Z'
      aux[1] = str(q)
      part = ''.join(aux)
      term4.append(part)
      aux[0] = 'Y'
      aux[1] = str(r)
      part = ''.join(aux)
      term4.append(part)
      aux[0] = 'Y'
      aux[1] = str(p)
      part = ''.join(aux)
      term4.append(part)
      term4.append(str(coeff))
      sterm=range(4)
      term1=' '.join(term1)
      term2=' '.join(term2)
      term3=' '.join(term3)
      term4=' '.join(term4)
      sterm[0]=term1
      sterm[1]=term2
      sterm[2]=term3
      sterm[3]=term4
    elif opt==2:
      aux=range(2)
      sterm=range(8)
      XY=['XXXX','XXYY','XYXY','YXXY','YXYX','YYXX','XYYX','YYYY']
      XYsigns=[1,-1,1,1,1,-1,1,1]
      for j in range(8):
        term1=[]
        coeff=''.join(str(XYsigns[j]*sign*coefficient/8.0))
        # term1
        for i in range(s+1,r):
          aux[0] = 'Z'
          aux[1] = str(i)
          part = ''.join(aux)
          term1.append(part) 
        for i in range(q+1,p):
          aux[0] = 'Z'
          aux[1] = str(i)
          part = ''.join(aux)
          term1.append(part) 
        aux[0] = XY[j][0]
        aux[1] = str(s)
        part = ''.join(aux)
        term1.append(part)
        aux[0] = XY[j][1]
        aux[1] = str(r)
        part = ''.join(aux)
        term1.append(part)
        aux[0] = XY[j][2]
        aux[1] = str(q)
        part = ''.join(aux)
        term1.append(part)
        aux[0] = XY[j][3]
        aux[1] = str(p)
        part = ''.join(aux)
        term1.append(part)
        term1.append(str(coeff))
        term1=' '.join(term1)
        sterm[j]=term1
  symterm=sterm
  return symterm

# Unit tests.
def main():

  # Test parameters.
  molecule = str(argv[1])
  basis = str(argv[2])
  try:
    cc_form = bool(int(argv[3]))
  except:
    cc_form = False
  try:
    run_type = str(argv[4])
  except:
    run_type = 'ucc'
  try:
    frozen_core = int(argv[5])
  except:
    frozen_core = None
  try:
    active_vir = int(argv[6])
  except:
    active_vir = None
  try:
    dirpath = str(argv[7])
  except:
    dirpath = None
  try:
    guesspath = str(argv[8])
  except:
    guesspath = None
  try:
    optimization_method = str(argv[9])
  except:
    optimization_method = None
  try:
    state_prep_method = str(argv[10])
  except:
    state_prep_method = 'exact'
  try:
    nTrotterSteps = int(argv[11])
  except:
    nTrotterSteps = 1
  try:
    guess = str(argv[12])
  except:
    guess = 'read'
  try:
    useHF = bool(int(argv[13]))
  except:
    useHF = True
  try:
    HFState = str(argv[14])
  except:
    HFState = '00110011'
  try:
    molproFile = str(argv[15])
  except:
    HFState = '00110011'

  print "active_vir and frozen_core:", frozen_core, active_vir

  # Get Hamiltonian
  print 'dirpath',dirpath
  # this code SHOULD make a reduction in the number of qubits, just the number of operators
  n_electrons = commutators.ElectronCount(molecule)
  n_orbitals = n_electrons+active_vir
  #cisd_hamiltonian, norbitals = GetHamiltonian(molecule, basis, n_excitations=2, n_orbitals=n_orbitals, path=dirpath)
  fci_hamiltonian, norbitals = GetHamiltonian(molecule, basis, n_excitations='FCI', n_orbitals=n_orbitals, path=dirpath)
  fci_energy, repulsion, coefficients, ham_terms = commutators.GetHamiltonianTerms(molecule, basis, n_orbitals=n_orbitals, path=guesspath)
  print 'Number of terms %i' % len(ham_terms)
  print "terms:", ham_terms
  print "n orbitals:", n_orbitals
  print "n electrons:", n_electrons
  basis_parsed = basis.split('-')

  print basis_parsed
  print basis_parsed[0]
  energies=[]

  # Get states and energies.
  fci_state = SparseDiagonalize(fci_hamiltonian)
  #print fci_state
  # cisd_state = SparseDiagonalize(cisd_hamiltonian)
  #print cisd_state
  # states=['00','01','10','11','11110000','00001111','10101010','01010101','11001100','00110011']
  # for state in states:
  #   print("initial state",state)
  # hf_state = HartreeFockState2(state)  
    # for n in range(0,2**(len(state))):
    #   if abs(hf_state[n,0])>0.00001:
    #     print(n,hf_state[n,0])
  if useHF == True:
    print("Use HF:",useHF)
    hf_state = HartreeFockState2(HFState)
  elif useHF == False:
    print("Use HF:",useHF)
    hf_state =  readMolproState(molproFile)
  print("HF state:")
  printState(hf_state)
  print("FCI state:")
  printState(fci_state)
  fci_energy_2 = Expectation(fci_hamiltonian, fci_state)
  #cisd_energy_2 = Expectation(fci_hamiltonian, cisd_state)
  hf_energy = Expectation(fci_hamiltonian, hf_state)
  print basis_parsed
  energies=[basis_parsed[0],str(hf_energy+repulsion),str(fci_energy+repulsion)]

  # Get terms
  initial_guess = []
  terms = []
  # JRF
  jw_terms = GetJordanWignerTerms(n_orbitals)
  counter=0
  coeffs=[]
  thres=1E-08

  if run_type=='troyer':
    
    for coefficient, term in zip(coefficients, ham_terms):
      # JRF
      counter=counter+1
      hterm=MatrixForm(1., term, jw_terms, add_conjugates=True, anti_hermitian=cc_form)
      hcontribution=Expectation(hterm, cisd_state)
      initial_guess += [hcontribution]
      terms += [term]

  elif run_type=='troyer-wc':

    # JRF
    for coefficient, term in zip(coefficients, ham_terms):
      # JRF
      hterm=MatrixForm(1., term, jw_terms, add_conjugates=True, anti_hermitian=cc_form)
      hcontribution=Expectation(hterm, cisd_state)
      if commutators.GetConjugate(term):
        counter=counter+1
        initial_guess += [hcontribution]
        terms += [term]

  elif run_type=='ucc':

    print "List of guesses: \n"
    print ham_terms
    print "\n"
    
    # I = scipy.sparse.csr_matrix([[1, 0], [0, 1]], dtype=complex)
    # Z = scipy.sparse.csr_matrix([[1, 0], [0, -1]], dtype=complex)
    # X = scipy.sparse.csr_matrix([[0, 1], [1, 0]], dtype=complex)
    # Y = scipy.sparse.csr_matrix([[0, -1j], [1j, 0]], dtype=complex)
    # pauli={'I':I,'X':X,'Y':Y,'Z':Z}

    for coefficient, term in zip(coefficients, ham_terms):
      # JRF
      value = keepTerm(term,frozen_core,n_orbitals)
      print term, value

      # if len(term)==2:
      #   t1='YX'
      #   t2='XY'
      # if len(term)==4:
      #   t1='XXXY'
      #   t2='YYYX'
      if value==True:
        # modify here to try ansatz without Z operators
        # hterm = MatrixFormNoZ(1., term, n_orbitals, pauli, t1, t2)
        hterm=MatrixForm(1., term, jw_terms, add_conjugates=True, anti_hermitian=cc_form)
        if guess=='zeros':
          hcontribution=0.0
        elif guess=='read':
          hcontribution=2.0*coefficient
        elif guess=='random':
          hcontribution=0.5*numpy.random.random()-0.25
          
        # print term, hcontribution
    # for term in operators:
    #   # JRF
    #   hterm=MatrixForm(1., term, jw_terms, add_conjugates=True, anti_hermitian=True)
    #   # Initial guess
    #   #hcontribution=coefficient
    #   # hcontribution=Expectation(hterm, cisd_state)
    #   #hcontribution=0.0

    #   # print "hamiltonian terms"
    #   # print ham_terms
    #   try:
    #     location = ham_terms.index(term)
    #     hcontribution = coefficients[location]
    #   except:
    #     hcontribution = 0.0

        counter=counter+1
        initial_guess += [hcontribution]
        terms += [term]
    
  print initial_guess
  print "Number of parameters:",len(initial_guess)

  filename1=molecule+basis+'.out'
  output1=open(filename1,'a')
  # output2=open('OptimalParameters_H3pd__CCD.out','a')
  print "repulsion:", repulsion
  print 'This is the number of terms for VQE:',counter
  # Obtain variational solution.
  print 'Analyzing %s in the %s basis.\n' % (molecule, basis)
  print 'Hartree-Fock energy is %s.' % repr(float(hf_energy+repulsion))
  print '\nExact (FCI) energy is %s.' % repr(float(fci_energy))
  ncalls, solution, variational_energy, myObjective = VariationallyMinimize(
     terms, hf_state, initial_guess, fci_hamiltonian, cc_form, n_orbitals, optimization_method, state_prep_method, nTrotterSteps=nTrotterSteps)
  print 'Variationally optimal energy is %s.\n' % repr(float(variational_energy+repulsion))
  #print "cisd and fci energies computed with this code:", fci_energy_2+repulsion, cisd_energy_2+repulsion
  print "fci energy:", fci_energy_2+repulsion
  overlap_fci = myObjective.returnOverlap(solution,fci_state)
  # overlap_cisd = myObjective.returnOverlap(solution,cisd_state)
  overlap_fcii = myObjective.returnOverlap(initial_guess,fci_state)
  print "the overlap with the fci wavefunction is:", overlap_fci
  print "initial overlap with the fci wavefunction is:", overlap_fcii
  # print "the overlap with the cisd wavefunction is:", overlap_cisd
  del myObjective
  diff=variational_energy+repulsion-fci_energy
  print 'Data %f %s %i %s.\n' % (thres,repr(float(variational_energy+repulsion)),counter,repr(float(diff)))
  print 'Optimal parameters \n'
  for t, ov in zip(terms,solution):
    if abs(ov)>1E-6:
      print t,ov
  print energies
  # print 'CISD energy is %s.' % repr(float(cisd_energy))


  # writing results
  energies.append(str(variational_energy+repulsion))
  energies.append(str(diff))
  energies.append(ncalls)
  energies.append(counter)
  output1.write((" ".join(str(i) for i in energies)+'\n'))
  output1.close()
  # output2.write((basis_parsed[0]+" ".join(str(i) for i in solution )+'\n'))
  # output2.close()

# Run.
if __name__ == '__main__':
  main()
