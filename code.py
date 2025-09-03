import numpy as np

# System of nonlinear equations (Mass + Energy conservation)
def equations(Q, R):
    Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8 = Q
    
    # ---- Mass Conservation ----
    eq1 = Q1 + Q2 + Q5 - 5              # Node 1
    eq2 = Q5 + Q7 - Q8 - Q6             # Node 2
    eq3 = Q3 + Q4 + Q6 - 3              # Node 3
    eq4 = Q2 + Q8 - Q3 - 2              # Node 4
    
    # ---- Energy Conservation (Loop equations) ----
    eq5 = R[0]*Q1*abs(Q1) - R[4]*Q5*abs(Q5) + R[6]*Q7*abs(Q7)
    eq6 = R[3]*Q4*abs(Q4) - R[5]*Q6*abs(Q6) - R[6]*Q7*abs(Q7)
    eq7 = R[4]*Q5*abs(Q5) - R[1]*Q2*abs(Q2) + R[7]*Q8*abs(Q8)
    eq8 = R[5]*Q6*abs(Q6) - R[7]*Q8*abs(Q8) - R[2]*Q3*abs(Q3)
    
    return np.array([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8])

# Jacobian matrix
def jacobian(Q, R):
    Q1, Q2, Q3, Q4, Q5, Q6, Q7, Q8 = Q
    d = lambda Ri, Qi: 2 * Ri * abs(Qi)   # derivative of R*Q*|Q|
    
    J = np.array([
        [1, 1, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, -1, 1, -1],
        [0, 0, 1, 1, 0, 1, 0, 0],
        [0, 1, -1, 0, 0, 0, 0, 1],
        
        [d(R[0],Q1), 0, 0, 0, -d(R[4],Q5), 0, d(R[6],Q7), 0],
        [0, 0, 0, d(R[3],Q4), 0, -d(R[5],Q6), -d(R[6],Q7), 0],
        [0, -d(R[1],Q2), 0, 0, d(R[4],Q5), 0, 0, d(R[7],Q8)],
        [0, 0, -d(R[2],Q3), 0, 0, d(R[5],Q6), 0, -d(R[7],Q8)]
    ])
    return J

# Newton-Raphson solver
def newton_raphson(Q, R, tol=1e-6, max_iter=200):
    for i in range(max_iter):
        F = equations(Q, R)
        if np.linalg.norm(F, ord=2) < tol:   # convergence check
            print(f"Converged in {i} iterations")
            return Q
        
        J = jacobian(Q, R)
        delta_Q = np.linalg.solve(J, -F)
        Q += delta_Q
        
    print("Did not converge within max iterations")
    return Q

# ------------------- MAIN -------------------
if __name__ == "__main__":
    # Value of p
    p = 40
    
    # Pipe resistances
    resistances = np.array([
        120 + p,  # R1
        200 + p,  # R2
        150 + p,  # R3
        300 + p,  # R4
        400 + p,  # R5
        400 + p,  # R6
        400 + p,  # R7
        400 + p   # R8
    ])
    
    # Initial guess for flows (closer to inflow/outflow balance)
    Q_initial = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    
    # Solve system
    flows = newton_raphson(Q_initial, resistances)
    print("Flow in each pipe:", np.round(flows, 6))
