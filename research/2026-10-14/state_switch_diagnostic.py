"""Demonstration of identifiability limits in an acne follicle state-reset hypothesis.
Synthetic states ONLY: not a fitted skin model or a clinical forecast.
The state variable x is dimensionless and MUST NOT be interpreted as lesions.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

OUT = Path(__file__).resolve().parent
SWITCH = 2.0/(3.0*np.sqrt(3))  # cubic saddle-node input magnitude

def equilibria(h):
    roots = np.roots([1, 0, -1, -h]) # x**3-x-h == 0
    return np.sort(np.real(roots[np.abs(np.imag(roots)) < 1e-9]))

def f(x, h, model, tau=15.0):
    return (h-x)/tau if model == 'slow_monostable' else x-x**3+h

def integrate_leg(model, start, end, duration, x0, tau=15.0, steps=401):
    def rhs(t, x):
        return [f(x[0], start+(end-start)*t/duration, model, tau)]
    grid = np.linspace(0, duration, steps)
    sol = solve_ivp(rhs, (0,duration), [x0], t_eval=grid, rtol=1e-9, atol=1e-10, max_step=1.)
    if not sol.success: raise RuntimeError(sol.message)
    return start+(end-start)*grid/duration, sol.y[0]

def ramp_area(model, duration, tau=15.0):
    h_low, h_high = -0.6, 0.6
    initial = h_low if model == 'slow_monostable' else equilibria(h_low)[0]
    h_up, x_up = integrate_leg(model, h_low, h_high, duration, initial, tau)
    _, x_down = integrate_leg(model, h_high, h_low, duration, x_up[-1], tau)
    area = float(abs(np.trapezoid(x_up-x_down[::-1], h_up)))
    return area, (h_up, x_up, x_down[::-1])

def post_washout_separation(model, t, tau=100.):
    # Different initial conditions, identical environmental input h=0, no treatment.
    if model == 'slow_monostable': return float(2*np.exp(-t/tau))
    return 2.0

def main():
    assert np.allclose(equilibria(0), [-1,0,1])
    assert len(equilibria(SWITCH+1e-5)) == 1
    assert len(equilibria(SWITCH-1e-5)) == 3
    cases=[]
    for model in ['slow_monostable', 'bistable_cubic']:
        for duration in [20, 80, 320, 1280]:
            area, curve = ramp_area(model,duration)
            cases.append(dict(model=model,ramp_duration_dimensionless=duration,hysteresis_loop_area_dimensionless=area))
    df = pd.DataFrame(cases)
    df.to_csv(OUT/'hysteresis_ramp_sweep.csv',index=False)
    follow=[]
    for t in [0, 10, 30, 100, 300, 1000]:
        follow.append(dict(post_washout_time_dimensionless=t,
                           monostable_tau100=post_washout_separation('slow_monostable',t),
                           monostable_tau1000=post_washout_separation('slow_monostable',t,tau=1000),
                           bistable=post_washout_separation('bistable_cubic',t)))
    pd.DataFrame(follow).to_csv(OUT/'washout_separation.csv',index=False)
    fig, ax = plt.subplots(figsize=(8,4.8))
    for model,sub in df.groupby('model'):
        ax.plot(sub['ramp_duration_dimensionless'],sub['hysteresis_loop_area_dimensionless'], marker='o',label=model)
    ax.set_xscale('log');ax.set_xlabel('Time per ramp (dimensionless)')
    ax.set_ylabel('Hysteresis loop area (arbitrary units)')
    ax.set_title('Synthetic state models: rate lag vs a bistable loop')
    ax.legend();fig.tight_layout();fig.savefig(OUT/'hysteresis_ramp_sweep.png',dpi=170);plt.close(fig)
    print('Saddle-node |input| threshold:',SWITCH)
    print(df.to_string(index=False))
    print('No result above is human biological data.')
    return df

if __name__ == '__main__': main()
