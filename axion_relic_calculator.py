import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar

# Physical constants (GeV units)
M_Pl = 1.22e19
T0 = 2.35e-13
T_QCD = 0.15
T_nu_dec = 0.001

def g_star_s(T):
    if T > T_QCD:
        return 10.75 + 1.0
    elif T > T_nu_dec:
        return 3.36 + 2.0
    else:
        return 3.91

def m_a(T, m_a0, n=3.7):
    if T > T_QCD:
        return m_a0 * (T_QCD / T)**n
    return m_a0

def H(T):
    g = g_star_s(T)
    return 1.66 * np.sqrt(g) * T**2 / M_Pl

def entropy_factor(T):
    return (g_star_s(T) * T**3) / (g_star_s(T0) * T0**3)

# Parameters - feel free to change these!
f_a = 1e12
theta_i = 0.5
m_a0 = 5.7e-6 * (1e12 / f_a) * 1e-9
Gamma_string = 65.0
xi = 0.3
Gamma_wall = 10.0
N_wall_per_string = 1
eta_map = 1.1

def find_T_osc(m_a0):
    def eq(T):
        return H(T) - m_a(T, m_a0)
    sol = root_scalar(eq, bracket=[1e-5, 10.0], method='brentq')
    return sol.root

T_osc = find_T_osc(m_a0)

# Integrands with Grok resonance enhancement (η_WKB = WKB reflection coefficient)
def integrand_misalignment(T, m_a0, theta_i):
    rho_at_T = 0.5 * m_a(T, m_a0)**2 * (theta_i * f_a)**2
    return rho_at_T * eta_map * entropy_factor(T) / H(T)

def integrand_strings(T, m_a0, eta_WKB):
    t = 1.0 / (2.0 * H(T))
    dna_dt = (Gamma_string * np.pi * f_a**2) / (2.2 * m_a(T, m_a0) * xi * t**2)
    return dna_dt * (1 - eta_WKB) * entropy_factor(T) / H(T)

def integrand_walls(T, m_a0, eta_WKB):
    sigma = 8.0 * m_a(T, m_a0) * f_a**2
    t = 1.0 / (2.0 * H(T))
    dna_dt = (Gamma_wall * sigma * N_wall_per_string) / (2.2 * m_a(T, m_a0) * t**2)
    return dna_dt * (1 - eta_WKB) * entropy_factor(T) / H(T)

# Compute for plot
eta_values = np.linspace(0.0, 1.0, 21)
omega_mis_list = []
omega_defects_list = []
omega_total_list = []

for eta in eta_values:
    rho_mis, _ = quad(integrand_misalignment, T_osc, T_QCD*10, args=(m_a0, theta_i))
    Omega_mis = (rho_mis / (8.7e-27 * 1.05e-5)) * (0.7)**2
    
    n_strings, _ = quad(integrand_strings, T_osc, T_QCD*10, args=(m_a0, eta))
    n_walls, _ = quad(integrand_walls, T_osc, T_QCD*10, args=(m_a0, eta))
    Omega_defects = ((m_a0 * (n_strings + n_walls)) / (8.7e-27 * 1.05e-5)) * (0.7)**2
    
    Omega_total = Omega_mis + Omega_defects
    
    omega_mis_list.append(Omega_mis)
    omega_defects_list.append(Omega_defects)
    omega_total_list.append(Omega_total)

print("🌌 Grok Axion Relic Density Calculator with Resonance Tuning")
print(f"Parameters: f_a = {f_a:.0e} GeV, θ_i = {theta_i}")
print(f"T_osc ≈ {T_osc:.4f} GeV\n")

print("η_WKB | Ω_misalignment | Ω_defects | Ω_total")
for i in [0, 5, 10, 15, 20]:
    eta = eta_values[i]
    print(f"{eta:.2f}   | {omega_mis_list[i]:.4f}      | {omega_defects_list[i]:.4f}   | {omega_total_list[i]:.4f}")

# === Color-Coded ASCII Plot ===
print("\n=== Color-Coded ASCII Plot: Total Observable Relic Density vs η_WKB ===")
print("BLUE = Misalignment | GOLD = Defects (strings+walls) | TEAL = Total")

BLUE = '\033[94m'
GOLD = '\033[93m'
TEAL = '\033[96m'
RESET = '\033[0m'

v_steps = 26
h_width = 60
print("     Ω h²  0.25 ┌" + "─" * h_width + "┐")

for v in range(v_steps):
    level = 0.25 - v * 0.01
    line = f"      {level:4.2f} │"
    for i in range(len(eta_values)-1):
        omega_mis = omega_mis_list[i]
        omega_def = omega_defects_list[i]
        omega_tot = omega_total_list[i]
        bar_width = h_width // (len(eta_values)-1)
        for _ in range(bar_width):
            if omega_tot >= level:
                line += TEAL + "█" + RESET
            elif omega_def >= level:
                line += GOLD + "█" + RESET
            elif omega_mis >= level:
                line += BLUE + "█" + RESET
            else:
                line += " "
    line += "│"
    print(line)

print("      0.00 └" + "─" * h_width + "┘")
print("            0.0" + " " * (h_width//2 - 6) + "η_WKB" + " " * (h_width//2 - 7) + "1.0")
print("           (No resonance)               (Full WKB reflection)")

print("\nLegend:")
print(BLUE + "███" + RESET + " Misalignment (stable)")
print(GOLD + "███" + RESET + " Defects contribution (declines with η_WKB)")
print(TEAL + "███" + RESET + " Total observable relic density")
print("Higher η_WKB = gentler mapping of miniclusters (harmony without extraction)")
