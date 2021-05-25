import random
import argparse
import time
import argparse
from collections import defaultdict
from time import process_time_ns
from attack_parameters import *
import sage.parallel.multiprocessing_sage

proof.arithmetic(false)
random.seed(2)
set_random_seed(2)

PARALLEL_CPUS = 4

def explore_dfs_optimized(E, basis, ell, f, c, returnCurve=False):
    """Generator for the DFS-exploration of isogenies
    Args:
        E (EllipticCurve): starting curve
        basis: tuple of points generating the ell^f-torsion of E
        ell (int): prime number
        f (int): exponent
        c (int): exponent (f + c = e)
    Yields:
        j, (phi, (m, n)), keys
        j: finite field element with the j-invariant of the codomain of phi
        phi: isogeny between E and an elliptic curve
(m, n): pair of integers generating phi
        keys: tqdm generator (for clean console-output purposes)
    """
    stack = []

    # Push initial points
    P, Q = basis[0], basis[1]

    keys = [(1, 0), (0, 1), (1, 1)]
    # shuffle(keys)
    for (m, n) in keys:
        R = m*P + n*Q
        S = P if m == 0 else Q
        phi = E.isogeny(ell ** (f - 1) * R)
        stack.append((phi, (m, n), 1, phi(R), phi(S)))
        
    coefs = list(range(ell))
    while len(stack) > 0:
        # Pop stack
        phi, (r, s), d, T, S = stack.pop()

        # Push new points, except the backtracking point
        if d < f:
            S *= ell
            shuffle(coefs)
            for x in coefs:
                m, n = r, s
                if m == 1:
                    n += x  * (ell ** d)
                else:
                    m += x  * (ell ** d)

                R = T + x * S
                psi = phi.codomain().isogeny(ell ** (f - d - 1) * R)
                stack.append((psi * phi, (m, n), d + 1, psi(R), psi(S)))

        # If depth f has been reached, yield
        else:
            if returnCurve:
                yield phi.codomain(), phi
            else:
                yield phi.codomain().j_invariant() #, (phi, (r, s))

def claw_with_dfs(E0, EA, degree):

    degree = Integer(degree)

    P, Q = get_points(E0, lM, degree, int(lc*lK^eK*lM^(eM - degree)))
    curves = set(explore_dfs_optimized(E0, (P, Q), lM, degree, degree))
    P1, Q1 = get_points(EA, Integer(lM), Integer(degree), int(lc*lK^eK*lM^(eM - degree)))
   
    for Et, phi_list in explore_dfs_optimized(EA, (P1, Q1), lM, degree, degree, True):
        if Et.j_invariant() in curves:
            return Et, list(phi_list)

    raise Exception("No collision found :(")


def claw(E0, EA, degree):
    degree = Integer(degree)

    t0 = time.time()
    P0, Q0 = get_points(E0, Integer(lM), Integer(degree), int(lc*lK^eK*lM^(eM - degree)))
    t1 = time.time()
    if verbose:
        print("claw get points took %0.2fs" % (t1 - t0))

    curves = set()

    t0 = time.time()
    for i in range(lM^degree):
        curves.add(isogeny(lM, degree, E0, P0 + i*Q0)[0].codomain().j_invariant())
    t1 = time.time()
    if verbose:
        print("claw collecting curves from E0 part 1 took %0.2fs" % (t1 - t0))

    t0 = time.time()
    for i in range(0, lM^degree, 2):
        curves.add(isogeny(lM, degree, E0, i*P0 + Q0)[0].codomain().j_invariant())
    t1 = time.time()
    if verbose:
        print("claw collecting curves from E0 part 2 took %0.2fs" % (t1 - t0))

    if verbose:
        print("Got the list of curves from E0")

    t0 = time.time()
    P1, Q1 = get_points(EA, Integer(lM), Integer(degree), int(lc*lK^eK*lM^(eM - degree)))
    t1 = time.time()
    if verbose:
        print("claw get points part 2 took %0.2fs" % (t1 - t0))

    for i in range(lM^degree):
        phi, phi_list = isogeny(lM, degree, EA, P1 + i*Q1)
        Etarget = phi.codomain()
        if Etarget.j_invariant() in curves:
            print("Found collision of degree", degree)
            return Etarget, phi_list

    for i in range(0, lM^degree, 2):
        phi, phi_list = isogeny(lM, degree, EA, i*P1 + Q1)
        Etarget = phi.codomain()
        if Etarget.j_invariant() in curves:
            print("Found collision of degree", degree)
            return Etarget, phi_list

    raise Exception("No collision found :(")

def claw_return_both(E0, EA, degree):

    degree = Integer(degree)

    t0 = time.time()
    P0, Q0 = get_points(E0, Integer(lM), Integer(degree), int(lc*lK^eK*lM^(eM - degree)))
    t1 = time.time()
    if verbose:
        print("claw get points took %0.2fs" % (t1 - t0))

    curves = {}

    t0 = time.time()
    for i in range(lM^degree):
        phi, phi_list = isogeny(lM, degree, E0, P0 + i*Q0)
        Et = phi.codomain()
        curves[hash(Et.j_invariant())] = phi_list
    t1 = time.time()
    if verbose:
        print("claw collecting curves from E0 part 1 took %0.2fs" % (t1 - t0))

    t0 = time.time()
    for i in range(0, lM^degree, 2):
        phi, phi_list = isogeny(lM, degree, E0, i*P0 + Q0)
        Et = phi.codomain()
        curves[hash(Et.j_invariant())] = phi_list
    t1 = time.time()
    if verbose:
        print("claw collecting curves from E0 part 1 took %0.2fs" % (t1 - t0))

    if verbose:
        print("Got the list of curves from E0")

    t0 = time.time()
    P1, Q1 = get_points(EA, Integer(lM), Integer(degree), int(lc*lK^eK*lM^(eM - degree)))
    t1 = time.time()
    if verbose:
        print("claw get points  part 2 took %0.2fs" % (t1 - t0))

    for i in range(lM^degree):
        phi, phi_list = isogeny(lM, degree, EA, P1 + i*Q1)
        Etarget = phi.codomain()
        if hash(Etarget.j_invariant()) in curves.keys():
            print("Found collision of degree", degree)
            return Etarget, curves[hash(Etarget.j_invariant())], phi_list

    for i in range(0, lM^degree, 2):
        phi, phi_list = isogeny(lM, degree, EA, i*P1 + Q1)
        Etarget = phi.codomain()
        if hash(Etarget.j_invariant()) in curves.keys():
            print("Found collision of degree", degree)
            return Etarget, curves[hash(Etarget.j_invariant())], phi_list

    raise Exception("No collision found :(")

def DL(bs,l,c,g):
    #costs sqrt(ell) operations in the finite field.
    for i in [0..(c-1)]:
        j= bs.index(l)
        if j != 0:         
            return i*c+j
        l= l* g;
    return 0

def DLx(bs,l,c,ell,e,L,g):
    n=0;
    for s in [0..(e-1)]:
        t=DL(bs,l^(ell^(e-1-s)),c,g)
        l =l /(L[s]^t)
        n =n+ t * ell^s
    return n;

def find_basis(E,P,Q,R,ell,e):
    w = P.weil_pairing(Q,ell^e)
    L = [w];
    wx=w
    for i in [1..(e-1)]:
        wx=wx^ell
        L.append(wx)
    
    bs = []
    g = wx.parent()(1)
    c=ceil(sqrt(ell))
    for i in [0..(c-1)]:
        bs.append(g)
        bs = sorted(bs)
        g  = g*wx 

    g = g^-1;
    l = P.weil_pairing(R,ell^e)
    nQ=DLx(bs,l,c,ell,e,L,g);
    assert l == w^nQ;
   
    l = R.weil_pairing(Q,ell^e)
    nP=DLx(bs,l,c,ell,e,L,g);
    assert l == w^nP;
    assert R == nP*P+nQ*Q;
    return nP,nQ;

def find_basis_identity(P,Q):
    if P == 0:
        return 1,0

    if Q == 0:
        return 0,1

    try:
        x = discrete_log(Q,P,lM^eM,operation='+')
        return -x, 1
    except Exception as e:
        x = discrete_log(P,Q,lM^eM,operation='+')
        return 1, -x

def two_isogeny_optimal(base, degree, E, R):
    if degree == 0:
        return

    if degree == 1:
        phi = E.isogeny(R)
        return phi
    else:
        right = int(degree/base)
        left = int(degree - right)
        Q = int(base^right)*R
        phi = two_isogeny_optimal(base, right, E, Q)
        return two_isogeny_optimal(base, left, phi.codomain(), phi(R)) * phi

def isogeny(base, degree, E, R):    
    Ea = E
    Ra = R
    l = []
    while degree > 0:
        phi = Ea.isogeny(base^(degree-1)*Ra)
        Ea = phi.codomain()
        Ra = phi(Ra)
        degree -= 1
        l.append(phi)
    return phi, l

def isogeny2(base, degree, E, K, P,Q,R):
    Ea = E
    Ka = K
    Ra = R
    Pa = P
    Qa = Q
    while degree > 0:
        phi = Ea.isogeny(base^(degree-1)*Ka)
        Ea = phi.codomain()
        Ka = phi(Ka)
        Pa = phi(Pa)
        Qa = phi(Qa)
        Ra = phi(Ra)
        degree -= 1
    return phi.codomain(), Pa, Qa, Ra

def isogenies(bases, degrees, E, R, P, Q):
    t0 = time.time()
    Ra = R
    Pa = P
    Qa = Q

    length = len(bases)
    i = 0
    while i < length:
        j = 0
        K =Ra
        while j < length:
            if i!=j:
                K = (bases[j]^degrees[j])*K
            j = j+1
        E,Pa,Qa,Ra = isogeny2(bases[i], degrees[i], E, K, Pa, Qa, Ra)
        i = i+1

    t1 = time.time()
    # print("`isogenies' with bases", bases, "and degrees", degrees, "took %0.2fs" % (t1 - t0))

    return E,Pa,Qa 

def get_points(E, base, exp, cofactor, P = -1):
    if base == lc:
        base = [l1, l2, l3]
        exp = [e1, e2, e3]
    else:
        base = [base]
        exp = [exp]

    order = 1
    for i in range(len(base)):
        order *= base[i]^(exp[i])

    if P == -1:
        is_p_wrong_order = True
        while is_p_wrong_order:
            P = E.random_point()
            P = cofactor * P
            is_p_wrong_order = False

            for b in base:
                is_p_wrong_order = is_p_wrong_order or (order//b)*P == 0

    is_q_wrong_order = True
    while is_q_wrong_order:
        Q = E.random_point()
        Q = cofactor * Q
        is_q_wrong_order = False

        for b in base:
            is_q_wrong_order = is_q_wrong_order or (order//b)*Q == 0

        if not is_q_wrong_order:
            wp = P.weil_pairing(Q, order)

            for b in base:
                wp_power = wp^(order//b)

                is_q_wrong_order = is_q_wrong_order or (wp_power == 1 or wp_power^b != 1)

    return P,Q

def solve_query(V, i):
    assert V.order() == lM^eM * lc

    t0 = time.time()
    Ev, Pv, Qv    = isogenies([lM,l1,l2,l3],[eM,e1,e2,e3],E,V,P,Q)
    Ev, PBv, QBv  = isogenies([lM,l1,l2,l3],[eM,e1,e2,e3],E,V,PK,QK)
 
    Evk, Pvk, Qvk, _ = isogeny2(lK, eK, Ev, PBv + beta*QBv, Pv, Qv, Ev(0))

    t1 = time.time()
    print("Solve query #%d took %0.2fs" % (i, t1 - t0))

    return (Evk, Pvk, Qvk)


def bit_reverse(n, width):
    b = '{:0{width}b}'.format(n, width=width)
    return int(b[::-1], 2)

def normalize_m_n(m, n, N):
    """
    Given the pair (m,n) modulo N, with one of m or n coprime to N, invert
    the  invertible component to create a pair (1, n') or (m', 1).

    Args:
        m (int)
        n (int)
        N (int)
    """

    if gcd(m, N) == 1:
        inv = inverse_mod(m, N)
    else:
        inv = inverse_mod(n, N)

    m, n = m * inv, n * inv
    return (m % N, n % N)

query_history = {}
curve_history = {}

def get_leaf_curves(PV, QV, q, use_Qvki, Pbase, Qbase):
    global query_history

    curves = []
    points  =[]

    coeffs = [0, 1]
    degree = eM//q

    for i in range(1, q):
        coeffs2 = []
        for c in coeffs:
            coeffs2.append(c)
            coeffs2.append(c + 2**(i*degree))
        coeffs = coeffs2

    coeffs = coeffs[:len(coeffs)/4]

    curves = []

    inputs = []

    for i, c in enumerate(coeffs):
        inputs.append((c, PV, QV, use_Qvki, i))
        # curves.append(compute_leaf(c, PV, QV, use_Qvki))

    results = list(compute_leaf(inputs))

    result_dict = {}

    for result in results:
        (f_input, f_output) = result
        result_dict[f_input[0][0]] = f_output

    for c in reversed(coeffs):
        curves.append(result_dict[c])

    return curves

@parallel(PARALLEL_CPUS)
def compute_leaf(coeff, PV, QV, use_Qvki, i):

    Vi = PV + coeff*QV
    (Evki, Pvki, Qvki) = solve_query(Vi, i)

    kernel = Qvki if use_Qvki else Pvki
    Evkvi, _, _ = \
        isogenies([l1, l2, l3], [e1, e2, e3], Evki, kernel, Evki(0), Evki(0))
    return Evkvi

@parallel(PARALLEL_CPUS)
def find_common_parent(E1,E2,degree, i):
    t0 = time.time()
    if algorithm == "DFS":
        c = claw_with_dfs(E1, E2, degree)
    else :
        c = claw(E1, E2, degree)
    t1 = time.time()
    print("Claw #%d took %0.2fs" % (i, t1 - t0))
    return c

def traverse(curves, backtrack, q):
    p = []
    length = len(curves)
    i = 0

    degree = int(eM/q)

    inputs = []

    while i < length:
        E1 = curves[i]
        i = i+1
        E2 = curves[i]

        inputs.append((E1, E2, degree, (i+1)//2))
        i = i+1
  

    result_dict = {}
    results = list(find_common_parent(inputs))

    for result in results:
        (f_input, f_output) = result
        result_dict[hash(f_input[0][0].j_invariant())] = f_output

    for i in range(0, len(curves), 2):
        j_E1 = hash(curves[i].j_invariant())
        p.append(result_dict[j_E1][0])

    last_E1 = hash(curves[-2].j_invariant())
    backtrack += result_dict[last_E1][1]

    return p, backtrack

def do_backtrack(phi, PA, QA):
    if isinstance(phi, sage.categories.map.FormalCompositeMap):
        first = phi.first()
        PA,QA = (do_backtrack(phi.then(), PA, QA))
        phihat = first.dual()
        PA = phihat(PA)
        QA = phihat(QA)
        return PA,QA
    else:
        phihat = phi.dual()
        PA = phihat(PA)
        QA = phihat(QA)
        return PA,QA

def apply_isogeny_list(isogenies, P, dual=False):
    PP = P
    isogenies = reversed(isogenies) if dual else isogenies
    for phi in isogenies:
        if dual:
            if phi.degree() == 1:
                continue
            phi = phi.dual()

        PP = phi(PP)

    return PP

def find_kernel(PV, QV, q, use_Qvki=False, Pbase=(-1, -1), Qbase=(-1, -1)):
    
    backtrack = []
    curves =  get_leaf_curves(PV, QV, q, use_Qvki, Pbase, Qbase)

    print("start traversing")
    while len(curves) > 1:
        curves, backtrack = traverse(curves, backtrack, q)

    Em, phi0, phi1 = claw_return_both(EK, curves[0], int(eM/q))

    print("start backtracking")
    (PA,QA) = get_points(EK, lM, eM, lc^ec*lK^eK)
    PA_K = PA
    QA_K = QA

    PA = apply_isogeny_list(phi0, PA)
    QA = apply_isogeny_list(phi0, QA)

    iso = PA.curve().isomorphism_to(phi1[-1].codomain())

    PA = iso(PA)
    QA = iso(QA)

    PA = apply_isogeny_list(phi1, PA, dual=True)
    QA = apply_isogeny_list(phi1, QA, dual=True)

    for phi in reversed(backtrack):
        PA,QA = do_backtrack(phi,PA,QA)

    print("Done backtrack")

    t0 = time.time()
    (k1,k2) = find_basis_identity(PA,QA)
 
    T0 = k1*PA_K  + k2*QA_K
    t1 = time.time()
    print("Found kernel, it took %0.2fs" % (t1 - t0))

    assert isogeny(lM, eM,EK, T0)[0].codomain().j_invariant() ==  backtrack[0].domain().j_invariant()
    return T0

def compute_V_basis(P, Q):

    FFF = Integers(lc)
    alpha = FFF(5*7*11)

    A = Matrix(FFF, [[alpha, 1], [1, alpha]])
    A = A.inverse()

    PM, QM = get_points(E, lM, eM, lK^eK * lc)

    PV = Integer(FFF(1/lM^eM)) * (Integer(A[0][0])*P + Integer(A[0][1])*Q) + PM
    QV = Integer(FFF(1/lM^eM)) * (Integer(A[0][1])*P + Integer(A[0][0])*Q) + QM

    assert (PV.order() == lc * lM ^eM and QV.order() == lc * lM ^eM)
    assert (P == lM^eM*(5*7*11*PV + QV))
    assert (Q == lM^eM*(PV + 5*7*11*QV))

    return PV, QV

# ###### PARAMETERS DEFINITION ######

# Console arguments parsing
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--verbose", help="Enable Verbose Mode", type=bool, default=false)
parser.add_argument("-e", "--expnum", help="Experiment to perform (see paramsets), defaults to 0", type=str, default=3)
parser.add_argument("-q", "--queries", help="Number of queries to perform, defaults to 4", type=int, default=4)
parser.add_argument("-a", "--algorithm", help="Choose between the MN and DFS algorithms to explore the isogeny graph (default is DFS)", type=str, default="DFS")
args = parser.parse_args()

verbose = args.verbose

lM = 2
lK = 3
l1 = 5
l2 = 7
l3 = 11

# Load parameters
load(f"paramsets/{args.expnum}.sage")
params = paramset["params"]

eM = params.eM
eK = params.eK
e1 = params.e1
e2 = params.e2
e3 = params.e3
ec = params.ec

lc = l1^e1*l2^e2*l3^e3
f = 1
queries = args.queries # the number of queries is 2^queries
algorithm = args.algorithm

print("Setting everything up...")
######  ELLIPTIC CURVE DEFINITION ######
p = lM^eM*lK^eK*lc^ec*f - 1
F.<i> = GF(p^2)
E = EllipticCurve(F, [1, 0])

t0 = time.time()
(P,Q) = get_points(E, lc, ec, lM^eM * lK^eK)
## P, Q are the torsion point whose image are revealed
t1 = time.time()
print("Generated P, Q took %0.2fs" % (t1 - t0))

t0 = time.time()
PK, QK = get_points(E, lK, eK, lc^ec * lM^eM)
# These torsion points are only used to compute K but are not public
t1 = time.time()
print("get_points took %0.2fs" % (t1 - t0))

t0 = time.time()
beta = random.randint(1, lK^eK)
K = PK + beta*QK
EK = isogeny(lK,eK,E,K)[0].codomain()

t1 = time.time()
print("Private E_K computations took %0.2fs" % (t1 - t0))

print("j(Ek) =", EK.j_invariant())

print("Complexity: lM = %d, eM = %d, security = %d, queries = %d" % (lM, eM, eM // 5 * 2, queries))
print("STARTING ATTACK")
# Section 1.1

t_start_attack = time.time()
t0 = time.time()

PV, QV = compute_V_basis(P, Q)
PVQV = PV + QV

t1 = time.time()
print("Generating (PV,QV) took %0.2fs" % (t1 - t0))

print("recovering PV part")
t_PV_start = time.time()
T0 = find_kernel(PV,QV, queries, False, (1, 0), (0, 1))
t_PV_end = time.time()
print("recovered PV part ==> OK took %0.2fs" % (t_PV_end - t_PV_start))

print("recovering QV part")
t_QV_start = time.time()
T1 = find_kernel(QV,PV, queries, True, (0, 1), (1, 0))
t_QV_end = time.time()
print("recovered QV part ==> OK took %0.2fs" % (t_QV_end - t_QV_start))

print("recovering PV+QV part")
t_PQV_start = time.time()
T2 = find_kernel(PVQV, QV, queries, False, (1, 1), (0, 1))
t_PQV_end = time.time()
print("recovered PV+QV part ==> OK took %0.2fs" % (t_PQV_end - t_PQV_start))

# Section 1.2

PVM = lc*PV
QVM = lc*QV

print("second part")
t_second_part_start = time.time()
e1 = T2.weil_pairing(T0, lM^eM)
e2 = T2.weil_pairing(T1, lM^eM)

e1_ = (PVM + QVM).weil_pairing(PVM, lM^eM)
e2_ = (PVM + QVM).weil_pairing(QVM, lM^eM)

zx = e1.log(e1_)
zy = e2.log(e2_)

zx = e1.log(e1_)
zy = e2.log(e2_)

M,_ = get_points(E, lM, eM, lc^ec*lK^eK)

print("challenge M ", M)
Z = Integers(lM^eM)
xdivy = int(Z(zx)/Z(zy))
k1,k2 = find_basis(E,PVM,QVM,M,lM,eM)
R = k1*T0 + k2*xdivy*T1

_, EMK, _, _ = isogeny2(lK, eK, E, K, M, E(0), E(0))

j2 = isogeny(lM, eM, EK, EMK)[0].codomain().j_invariant()
j1 = isogeny(lM, eM, EK, R)[0].codomain().j_invariant()

assert j1 == j2
t_second_part_end = time.time()
t_end_attack = time.time()
print("attack succedeed!!! took %0.2fs" % (t_end_attack - t_start_attack))
print("====SUMMARY====")
print("Complexity: lM = %d, eM = %d, security = %d, queries = %d" % (lM, eM, eM // 5 * 2, queries))
print("recovered PV part ==> OK took %0.2fs" % (t_PV_end - t_PV_start))
print("recovered QV part ==> OK took %0.2fs" % (t_QV_end - t_QV_start))
print("recovered PV+QV part ==> OK took %0.2fs" % (t_PQV_end - t_PQV_start))
print("second part took %0.2fs" % (t_second_part_end - t_second_part_start))
