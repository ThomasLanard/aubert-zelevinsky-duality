# Aubert–Zelevinsky Duality Algorithms

This repository provides a unified implementation of the **Aubert–Zelevinsky involution** for representations of `GL(n)` and the classical groups `Sp(2n)`,`SO(2n+1)`.

## 🔁 Unified Interface: `AD`

```python
AD(case=None, *args, verbose=False)
```

A function for computing the dual Langlands data, supporting the following cases:

- `"gl"` – for GL(n) (Mœglin–Waldspurger algorithm)
- `"good"` – classical group with ρ of good parity
- `"bad"` – classical group with ρ of bad parity
- `"ugly"` – classical group with ρ ugly (use ρ and ρᵛ (the contragredient))

If `case` is omitted, it is automatically inferred from the arguments when possible.

---

## 💡 Input Formats

### ▶️ GL(n)

```python
m = ([a1, b1], [a2, b2], ..., [ak, bk])
AD("gl", m)
# or simply:
AD(m)
```

Each segment `[a, b]` satisfies:
- `a ≤ b`
- `b - a ∈ Z`
- all `a_i` lie on the same affine line: `a_i - a_j ∈ Z`

---

### ▶️ Classical groups - Good Parity

```python
m = ([a1, b1], ..., [ak, bk])
T = ([x1, ε1], ..., [xr, εr])  # where εi ∈ {+1, -1}
AD("good", m, T)
```

Conditions:
- All `a`, `b`, and `x ∈ T` are in Z or (1/2)Z-Z.
- All segments satisfy `a + b < 0`.
- Each sign `εi` is either `+1` or `-1`.
- If `(xi, εi)` and `(xj, εj)` are in `T` with `xi = xj` then `εi = εj`

---

### ▶️ Classical groups - Bad Parity

```python
m = ([a1, b1], ..., [ak, bk])
T = (x1, x2, ..., xr)
AD("bad", m, T)
```

Conditions:
- All `a`, `b`, and `x ∈ T` are in Z or (1/2)Z-Z.
- Each `x ∈ T` appears an **even** number of times.
- All segments satisfy `a + b < 0`.

---

### ▶️ Classical groups - case Ugly

```python
m = ([a1, b1], ..., [ak, bk]) # multisegment supported in ℤρ
mcontr = ([c1, d1], ..., [cl, dl]) # multisegment supported in ℤρᵛ
T = (x1, x2, ..., xr)
AD("ugly", m, mcontr, T)
```

Conditions:
- All segments satisfy `a + b < 0` and `c + d < 0` .
- All the coefficients are in the same (Zρ u Zρᵛ) - line.

---

## 🧠 Case Inference

The case `"gl"`, `"bad"`, `"good"` or `"ugly"` can be omitted and inferred from the number and type of arguments:

```python
AD(m)                 # → inferred as "gl"
AD(m, T)              # → "bad" or "good" depending on structure of T (if T is not empty)
AD(m, mcontr, T)      # → "ugly"
```

If `verbose=True`, the inferred case will be printed before execution.

---

## 📢 Verbose Mode

Use `verbose=True` to print step-by-step output during the computation. Example:

```python
m=([-3, -3], [-2, -1]) 
T=([1, 1],)
AD("good",m,T,verbose=True)
```

Output:

```
Good Case — Step 0
  Symmetrisation
  Remaining: ([-3, -3], [-2, -1], [-1, 1], [3, 3], [1, 2])
  Sign updates: ε([-1, 1]) = 1

Good Case — Step 1
  Labeled segments: ([3, 3, 1], [1, 2, 1], [-1, 1, 0], [-2, -1, -1], [-3, -3, -1])
  Dual: ([1, 3], [-3, -1])
  Remaining: ([1, 1], [0, 0], [-1, -1])
  Sign updates: ε([0, 0]) = 1

Good Case — Step 2
  Labeled segments: ([1, 1, 1], [0, 0, 0], [-1, -1, -1])
  Dual: ([-1, 1],)
  Sign dual: ε([-1, 1]) = 1
  Remaining: ()
```

---

## 🐍 Python / Sage Usage

This project is compatible with both **Python** and **SageMath**. To use it:

```python
from aubert_dual import AD
```

---

## 📚 Bibliography

The implementation is based on:

- T. Lanard and A. Mínguez, *A Moeglin-Waldspurger algorithm for the Aubert-Zelevinsky duality*.

---

## 📁 Files

- `aubert_dual.py` – main implementation of the algorithms
- `README.md` – documentation

