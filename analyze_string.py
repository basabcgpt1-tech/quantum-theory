"""String tension vs mass perturbation theory.
First order (exact):  sigma = S m (1 - cos 2 pi alpha),  S = e^gamma/(2 pi^1.5) g.
Second order (Adam 1997; constants quoted from memory, sign fixed by requiring a
negative second-order vacuum energy):  -0.358 m^2 at alpha = 1/4,  0 at alpha = 1/2."""
import json
import numpy as np

d = json.load(open("data/string.json"))
S = np.exp(np.euler_gamma) / (2 * np.pi ** 1.5)
E_plus = -8.9139
b2 = lambda al: -(S * np.sqrt(np.pi)) ** 2 / 4 * E_plus * (np.cos(4 * np.pi * al) - 1)
xs = sorted({float(k.split("|")[0]) for k in d})
ms = sorted({float(k.split("|")[1]) for k in d} - {0.0})
for al in (0.25, 0.5):
    print(f"\nalpha = {al}:   theory a = {S * (1 - np.cos(2 * np.pi * al)):.4f}, b = {b2(al):.3f}")
    for x in xs:
        if f"{x}|{ms[-1]}|{al}" not in d:
            continue
        s = np.array([d[f"{x}|{m}|{al}"] for m in ms])
        th = np.array([S * m * (1 - np.cos(2 * np.pi * al)) + b2(al) * m * m for m in ms])
        b, a = np.polyfit(ms[:4], s[:4] / np.array(ms[:4]), 1)
        print(f"  x={x}: residual at m=0: {d[f'{x}|0.0|{al}']:+.5f};  fit a={a:.4f} b={b:.3f}")
        for m, v, t in zip(ms, s, th):
            print(f"     m/g={m:5.3f}  lattice {v:.5f}  theory(1st+2nd) {t:.5f}  ratio {v / t:.3f}")
