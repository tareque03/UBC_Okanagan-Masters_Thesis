"""Figure 4.2: structured SIR model, slow-fast decomposition with a fixed rho2 = 0.2
over a 60-day period (rho1 is chosen at each step so that the Perron root of Q is 0)."""
import numpy as np
import matplotlib.pyplot as plt

# Parameters
epsilon = 0.01
S_i = np.array([0.6, 0.4])
S0 = (1 - epsilon) * S_i
I0 = epsilon * np.array([0.6, 0.4])
C = np.array([[6.03, 3], [3, 4.08]])
p = np.array([0.7, 0.3])
R0 = np.array([0, 0])
beta = np.matmul(C, np.diag(p))
gamma = np.array([[0.4, 0], [0, 0.2]])
rho2 = 0.2

T = 60      # duration of simulation
h = 0.05    # time step
time = np.arange(0, T, h)
num_iterations = len(time)


def simulate_SI_model(S0, I0, R0, beta, gamma, num_iterations, h, rho2):
    S = np.zeros((num_iterations, len(S0)))
    I = np.zeros((num_iterations, len(I0)))
    R = np.zeros((num_iterations, len(R0)))
    P = np.zeros((num_iterations, 2, 2))
    Q = np.zeros((num_iterations, 2, 2))
    S[0] = S0
    I[0] = I0
    R[0] = R0
    for n in range(0, num_iterations - 1):
        rho1 = (((beta[0][1]*beta[1][0]*S[n][0]*S[n][1]*rho2*rho2) / (rho2*beta[1][1]*S[n][1]-gamma[1][1]))+gamma[0][0]) / (beta[0][0]*S[n][0])
        if not (0 <= rho1 <= 1):
            break
        beta_new = np.array([[beta[0][0] * rho1, beta[0][1] * rho2],
                             [beta[1][0] * rho2, beta[1][1] * rho2]])
        Q[n] = beta_new * S[n] - gamma
        eigenvalues1, eigenvectors1 = np.linalg.eig(Q[n])
        dominant_eigenvector1 = eigenvectors1[:, np.argmax(eigenvalues1)]
        eigenvalues2, eigenvectors2 = np.linalg.eig(Q[n].T)
        dominant_eigenvector2 = eigenvectors2[:, np.argmax(eigenvalues2)]
        P[n] = np.outer(dominant_eigenvector1, dominant_eigenvector2.conj()) / np.matmul(dominant_eigenvector2.conj(), dominant_eigenvector1)
        S[n + 1] = S[n] + h * (-np.matmul(I[n], beta_new * S[n]))
        I[n + 1] = np.matmul(P[n], I[n])
        R[n + 1] = R[n] + h * (np.matmul(I[n], gamma))
    return S, I, R


S, I, R = simulate_SI_model(S0, I0, R0, beta, gamma, num_iterations, h, rho2)

fig = plt.figure(figsize=(8, 8))

plt.subplot(3, 1, 1)
plt.plot(time, S[:, 0], label='S-non-seniors')
plt.plot(time, S[:, 1], label='S-seniors')
plt.xlabel('Time'); plt.ylabel('S'); plt.title('SIR Model Simulation - S'); plt.legend()

plt.subplot(3, 1, 2)
plt.plot(time, I[:, 0], label='I-non-seniors')
plt.plot(time, I[:, 1], label='I-seniors')
plt.xlabel('Time'); plt.ylabel('I'); plt.title('SIR Model Simulation - I'); plt.legend()

plt.subplot(3, 1, 3)
plt.plot(time, R[:, 0], label='R-non-seniors')
plt.plot(time, R[:, 1], label='R-seniors')
plt.xlabel('Time'); plt.ylabel('R'); plt.title('SIR Model Simulation - R'); plt.legend()

plt.tight_layout()
plt.savefig('fig4_2_output.png', dpi=150)
plt.show()


# ---- Save each panel (S, I, R) as its own plot in ../plots ----
import os
os.makedirs('../plots', exist_ok=True)
for name, data, ylabel in [('S', S, 'S'), ('I', I, 'I'), ('R', R, 'R')]:
    plt.figure(figsize=(8, 3.5))
    plt.plot(time, data[:, 0], label=f'{name}-non-seniors')
    plt.plot(time, data[:, 1], label=f'{name}-seniors')
    plt.xlabel('Time'); plt.ylabel(ylabel); plt.title(f'SIR Model Simulation - {name}'); plt.legend()
    plt.tight_layout()
    plt.savefig(f'../plots/fig4_2_{name}.png', dpi=150)
    plt.close()
