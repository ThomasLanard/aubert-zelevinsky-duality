###########################################################
# Aubert–Zelevinsky Involution Algorithms
###########################################################


# =========================================================
# Imports
# =========================================================

from typing import List, Tuple, Sequence, Dict, Any, Union
from collections import defaultdict, Counter
from functools import wraps
from contextlib import contextmanager

# =========================================================
# Types
# =========================================================

Segment = List[Any]
Multisegment = Sequence[Segment]
LanglandsParameter = Tuple[Any, ...]
EnhancedParameter = Tuple[Tuple[Any, int], ...]
LanglandsData = Tuple[Multisegment, LanglandsParameter]
LabeledSegment = List[Any]

# =========================================================
# Global Validation Control
# =========================================================

VALIDATION_ENABLED = True

def set_validation(enabled: bool) -> None:
    """
    Enable or disable global validation logic.

    Args:
        enabled (bool): If True, enable validation. If False, disable it.
    """
    global VALIDATION_ENABLED
    VALIDATION_ENABLED = enabled

@contextmanager
def disable_validation():
    """
    Context manager to temporarily disable validation inside a with-block.
    """
    global VALIDATION_ENABLED
    old_value = VALIDATION_ENABLED
    VALIDATION_ENABLED = False
    try:
        yield
    finally:
        VALIDATION_ENABLED = old_value

# =========================================================
# Validation Utilities
# =========================================================

def validated(validator):
    """
    Decorator to apply a validator function to the arguments of an AD_* function.
    Automatically applies validation before the function runs.
    """
    def decorator(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            validator(*args)
            return f(*args, **kwargs)
        return wrapper
    return decorator

def is_integer(x):
    """
    Check if x is an integer
    """
    return x == int(x) if isinstance(x, (int, float)) else hasattr(x, 'is_integer') and x.is_integer()

def validate_data_gl(mult: Multisegment) -> None:
    """
    Validates that:
      1. Each segment [a, b] satisfies a ≤ b and (b - a) is an integer.
      2. All segment starting points lie on the same affine integer lattice:
         for all i, j, a_i - a_j ∈ ℤ.
    
    Raises:
        ValueError if any of these conditions are violated.
    """
    if not mult:
        return  # empty multisegment is valid

    for seg in mult:
        if len(seg) != 2:
            raise ValueError(f"Invalid segment format: {seg}")
        a, b = seg
        if a > b:
            raise ValueError(f"Invalid segment {seg}: beginning > end")
        if not is_integer(b - a):
            raise ValueError(f"Invalid segment {seg}: end - beginning is not an integer")

    base_beginning = mult[0][0]
    for seg in mult:
        a = seg[0]
        if not is_integer(a - base_beginning):
            raise ValueError(f"Segments are not in the same ℤ-line: {a} - {base_beginning} is not an integer")

def validate_langlands_data(m: Multisegment, T: Tuple[Any, ...]) -> None:
    """
    Validates Langlands data conditions:

    Conditions:
        1. m is a valid multisegment (via validate_data_gl)
        2. Each segment [a, b] satisfies a + b < 0
        3. All elements of T are in (1/2)ℤ
        4. All elements of T are ≥ 0
        5. For all x, y ∈ T: x - y ∈ ℤ
        6. For all x ∈ T and all [a, b] ∈ m: x - a ∈ ℤ

    T may consist of numbers or pairs (x, ±1).
    """
    def is_half_integer(x): return is_integer(2 * x)

    def extract_x(t_elem):
        return t_elem[0] if isinstance(t_elem, (tuple, list)) else t_elem

    # 1. m is valid multisegment
    validate_data_gl(m)

    # 2. a + b < 0
    for a, b in m:
        if a + b >= 0:
            raise ValueError(f"Segment [{a}, {b}] does not satisfy a + b < 0")

    x_list = [extract_x(t) for t in T]

    # 3. x ∈ (1/2)ℤ
    for x in x_list:
        if not is_half_integer(x):
            raise ValueError(f"Element {x} in T is not in (1/2)ℤ")

    # 4. x ≥ 0
    for x in x_list:
        if x < 0:
            raise ValueError(f"Element {x} in T is negative")

    # 5. x - y ∈ ℤ for all x, y ∈ T
    for i in range(len(x_list)):
        for j in range(i + 1, len(x_list)):
            if not is_integer(x_list[i] - x_list[j]):
                raise ValueError(f"Elements in T are not in the same ℤ-line: {x_list[i]} - {x_list[j]} is not an integer")

    # 6. x - a ∈ ℤ for one x ∈ T and one segment [a, b] ∈ m
    if m and x_list:
        a = m[0][0]
        x = x_list[0]
        if not is_integer(x - a):
            raise ValueError(f"T and m are not in the same ℤ-line: {x} - {a} is not an integer")
            
def validate_data_ugly(
    m: Multisegment,
    m_contr: Multisegment,
    T: LanglandsParameter,
) -> None:
    """
    Validation for the ugly case.

    Adds:
        For all [a, b] ∈ m, [c, d] ∈ m_contr: c + a ∈ ℤ
    """
    validate_langlands_data(m, T)
    validate_langlands_data(m_contr, T)

    # c + a ∈ ℤ for all [a, b] in m and [c, d] in m_contr
    if m and m_contr:
        a = m[0][0]
        c = m_contr[0][0]
        if not is_integer(a + c):
            raise ValueError(f"m and m_contr are not in the same ℤ-line: {a} + {c} is not an integer")

def validate_data_bad(m: Multisegment, T: Tuple[Any, ...]) -> None:
    """
    Validation for the bad parity case.

    Adds:
        1. All coefficients in m and T are in (1/2)ℤ
        2. Each element of T appears an even number of times.
    """
    validate_langlands_data(m, T)
        
    def is_half_integer(x):
        return is_integer(2 * x)
    
    # (1) All coefficients in (1/2)ℤ
    for a, b in m:
        if not is_half_integer(a):
            raise ValueError(f"Segment [{a}, {b}] has non-(1/2)-integer endpoints")

    # (2) Even multiplicity for each x in T
    counts = Counter(T)
    for x, c in counts.items():
        if c % 2 != 0:
            raise ValueError(f"Element {x} appears an odd number of times in T")

def validate_data_good(m: Multisegment, T: EnhancedParameter) -> None:
    """
    Validation for the good parity case.

    Adds:
        1. Each [x, e] ∈ T satisfies e ∈ {+1, -1}
        2. If [x1, e1], [x2, e2] ∈ T and x1 == x2, then e1 == e2.
    """
    validate_langlands_data(m, T)

    seen = {}
    for x, e in T:
        # (1) e ∈ {+1, -1}
        if e not in {+1, -1}:
            raise ValueError(f"Invalid sign {e} for {x}: expected +1 or -1")
        
        # (2) If [x1, e1], [x2, e2] ∈ T and x1 == x2, then e1 == e2.
        if x in seen and seen[x] != e:
            raise ValueError(f"Inconsistent signs for {x} in T: got both {seen[x]} and {e}")
        seen[x] = e

# =========================================================
# Sorting Utilities
# =========================================================

def index_max_end(m: List[Segment]) -> int:
    """
    Return the index of the segment with the maximal end value in the list m.

    Raises:
        ValueError: if m is empty.
    """
    if not m:
        raise ValueError("Empty multisegment provided")
    return max(range(len(m)), key=lambda i: m[i][1])

def multisegment_sort_key(seg: Segment) -> Tuple:
    """
    Sorting key for multisegments: sort by a, then by -b.
    """
    a, b = seg
    return (a, -b)

def sort_multisegment(m: List[Segment]) -> List[Segment]:
    """
    In-place sort of multisegments using multisegment_sort_key.
    """
    m.sort(key=multisegment_sort_key, reverse=True)

def labeled_segment_sort_key(seg: LabeledSegment) -> Tuple:
    """
    Return a key for sorting labeled segments [b, e, t] as follows:
      - First by type t (ascending)
      - Then:
          - if t == ±1: by b ascending, then e descending
          - if t == 0 : by b descending (equivalent to e ascending since b = -e)
    """
    b, e, t = seg
    if t in {1, -1}:
        return (t, b, -e)
    else:
        return (t, -b)
    
def sort_labeled_segments(segments: List[LabeledSegment]) -> None:
    """
    In-place sort of labeled segments using labeled_segment_sort_key.
    """
    segments.sort(key=labeled_segment_sort_key, reverse=True)

# =========================================================
# Verbose / Debug Utils
# =========================================================

def _print_step(
    tag: str,
    step: int,
    labeled_seg: Any = None,
    dual: Any = None,
    sign_dual : Any = None,
    remaining: Any = None,
    signs: Any = None,
    message: str = None
) -> None:
    """
    Unified and structured step printer for verbose output.

    Args:
        tag: algorithm name (e.g., "AD_gl", "AD_ugly", "AD_bad", "AD_good")
        step: step number
        labeled_seg: labeled segments (specific to AD_good)
        dual: newly added dual segment(s)
        sign_dual: signs associated to the new dual segment(s)
        remaining: remaining multisegment
        signs: updated sign dictionary (if applicable)
        message: optional custom message for the step
    """
    print(f"{tag} — Step {step}")
    if message:
        print(f"  {message}")
    if labeled_seg is not None:
        print(f"  Labeled segments: {labeled_seg}")
    if dual is not None:
        print(f"  Dual: {dual}")
    if sign_dual:
        print(f"  Sign dual: {sign_dual}")
    if remaining is not None:
        print(f"  Remaining: {remaining}")
    if signs:
        print(f"  Sign updates: {signs}")
    print()

def format_sign_dict(E):
    """
    Format a dictionary of signs indexed by doubled half-integers.

    In the 'good' case, signs are stored using integer keys representing 2x
    (e.g. x = 1/2 is stored as key 1). This function converts these keys back
    to half-integers and formats them in a readable string of the form:

        "ε(1/2) = -1,  ε(1) = +1,  ε(3/2) = -1"

    Args:
        E (dict[int, int]): Dictionary mapping int(2x) to ±1 signs.

    Returns:
        str: A human-readable string describing the sign function ε.
    """
    def format_key(k):
        if k % 2 == 0:
            return str(k // 2)
        else:
            return f"{k}/2"

    return ", ".join(f"ε([{format_key(-k)}, {format_key(k)}]) = {v}" for k, v in sorted(E.items()))

# =========================================================
# Case: GL(n)
# =========================================================

@validated(validate_data_gl)
def AD_gl(
    mult: Multisegment,
    verbose: bool = False,
    print_step: callable = None
) -> Tuple[Segment, ...]:
    """
    Compute the Aubert-Zelevinsky dual of a multisegment using the Moeglin-Waldspurger algorithm.

    Parameters:
        mult: sequence (tuple or list) of segments [beginning, end]
        verbose: if True, print each step of the algorithm
        print_step: optional function for formatted output. Signature:
            print_step(step: int, dual=[...], remaining=[...])

    Returns:
        Tuple of segments representing the dual.
    """
    # Copy input to mutable list of lists
    mutable_segments = [[a, b] for a, b in mult]

    # Initial sort
    sort_multisegment(mutable_segments)
    
    dual = []

    step = 1
    while mutable_segments:
        i = index_max_end(mutable_segments)
        e = mutable_segments[i][1]
        e0 = e
        mutable_segments[i][1] = e - 1

        # Create the linked chain of segments and modify mutable_segments
        for j in range(i + 1, len(mutable_segments)):
            if mutable_segments[j][1] == e - 1:
                mutable_segments[j][1] = e - 2
                e -= 1

        dual_seg = [e, e0]
        dual.append(dual_seg)

        # Filter out empty/invalid segments
        new_segments = []
        for seg in mutable_segments:
            if seg[1] >= seg[0]:
                new_segments.append(seg)
        mutable_segments = new_segments

        # Verbose output
        if print_step is not None:
            print_step(step, dual=[dual_seg], remaining=tuple(mutable_segments))
        elif verbose:
            _print_step("GL", step, dual=dual_seg, remaining=tuple(mutable_segments))
        step += 1

    return tuple(dual)

# =========================================================
# Case: ugly
# =========================================================

def symmetrisation_ugly(
    m: Multisegment,
    m_contr: Multisegment,
    T: LanglandsParameter,
) -> List[Segment]:
    """
    Build the symmetric Langlands data for the ugly case.

    Args:
        m: multisegment for ρ
        m_contr: multisegment for contragredient of ρ
        T: common Langlands parameter

    Returns:
        A symmetric multisegment (list of segments)
    """
    s = []
    s += [seg[:] for seg in m] 
    s += [[-x, x] for x in T]
    s += [[-seg[1], -seg[0]] for seg in m_contr]
    return s

def langlands_data_ugly(
    symmetric_multisegment: List[Segment]
) -> Tuple[LanglandsData, LanglandsData]:
    """
    Extract Langlands data (m, m_contr, T) from a symmetric multisegment.

    Args:
        symmetric_multisegment: list of segments [a, b] (a symmetric multisegment)

    Returns:
        (m, m_contr, T) - Langlands data
    """
    m = []
    m_contr = []
    T = []

    for start, end in symmetric_multisegment:
        center = start + end
        if center < 0:
            m.append([start, end])
        elif center == 0:
            T.append(end)
        else:
            m_contr.append([-end, -start])

    return (tuple(m), tuple(m_contr), tuple(T))

@validated(validate_data_ugly)
def AD_ugly(
    m: Multisegment,
    m_contr: Multisegment,
    T: LanglandsParameter,
    verbose: bool = False
) -> Tuple[LanglandsData, LanglandsData]:
    """
    Compute the Aubert–Zelevinsky dual in the ugly case.

    Args:
        m: multisegment for ρ
        m_contr: multisegment for contragredient of ρ
        T: common Langlands parameter
        verbose: if True, print intermediate steps
    Returns:
        (m_dual, m_dual_contr, T_dual) – dual Langlands data
    """
    #Symmetrisation of the datas
    sym = symmetrisation_ugly(m,m_contr,T)

    print_step = (lambda step, **kwargs: _print_step("Ugly Case", step, **kwargs)) if verbose else None

    if print_step:
        print_step(0, message="Symmetrisation", remaining=tuple(sym))

    #Applying the Moeglin-Waldspurger algorithm
    sym_dual = AD_gl(sym, print_step=print_step)
    
    return langlands_data_ugly(sym_dual)

# =========================================================
# Case: bad
# =========================================================

def symmetrisation_bad(m: Multisegment, T: LanglandsParameter) -> List[Segment]:
    """
    Build the symmetric multisegment for the bad parity case.

    Args:
        m: multisegment
        T: Langlands parameter (tuple of symmetric exponents)

    Returns:
        A symmetric multisegment (list of segments)
    """
    s = []
    s += [seg[:] for seg in m]
    s += [[-x, x] for x in T]
    s += [[-seg[1], -seg[0]] for seg in m]
    return s

def langlands_data_bad(symmetric_multisegment: List[Segment]) -> LanglandsData:
    """
    Extract (m, T) from a symmetric multisegment in the bad case.

    Args:
        symmetric_multisegment: list of segments

    Returns:
        Langlands data (m, T)
    """
    m = []
    T = []
    for start, end in symmetric_multisegment:
        center = start + end
        if center < 0:
            m.append([start, end])
        elif center == 0:
            T.append(end)
    return tuple(m), tuple(T)

def modify_segment_list(segments: List[Segment], seg: Segment) -> None:
    """
    Modify the list of segments in place:
    - Remove seg and its symmetric counterpart
    - Add their truncated versions if possible
    """
    segments.remove([seg[0], seg[1]])
    segments.remove([-seg[1], -seg[0]])
    if seg[0] <= seg[1] - 1:
        segments.append([seg[0], seg[1] - 1])
        segments.append([-seg[1] + 1, -seg[0]])

@validated(validate_data_bad)      
def AD_bad(m: Multisegment, T: LanglandsParameter, verbose: bool = False) -> LanglandsData:
    """
    Compute the Aubert–Zelevinsky dual in the bad case.

    Args:
        m: multisegment
        T: Langlands parameter
        verbose: if True, print each dual pair constructed

    Returns:
        Langlands data (m', T')
    """
    #Symmetrisation of the datas
    mutable_segments = symmetrisation_bad(m, T)

    #Ordering
    sort_multisegment(mutable_segments)

    if verbose:
        _print_step("Bad Case", 0, message="Symmetrisation", remaining=tuple(mutable_segments))

    dual = []

    step=1

    while mutable_segments:
        i = index_max_end(mutable_segments)
        e = mutable_segments[i][1]
        e0 = e

        #Computing the chain of segments
        chain_segments = [mutable_segments[i]]

        for j in range(i + 1, len(mutable_segments)):
            seg = mutable_segments[j]
            if seg[1] == e - 1:
                if ([-seg[1], -seg[0]] not in chain_segments) or (mutable_segments.count(seg) > 1):
                    chain_segments.append(seg)
                    e = seg[1]
        
        #Modification of mutable_segments with chain_segment
        for seg in chain_segments:
            modify_segment_list(mutable_segments, seg)

        dual_pair = [[e, e0], [-e0, -e]]
        dual.extend(dual_pair)
        
        #Reordering
        sort_multisegment(mutable_segments)

        if verbose:
            _print_step("Bad Case", step, dual=tuple(dual_pair), remaining=tuple(mutable_segments))
            step += 1
    return langlands_data_bad(dual)

# =========================================================
# Case: good
# =========================================================

def symmetrisation_good(m: Multisegment, T: EnhancedParameter) -> Tuple[List[Segment], Dict[Any, int]]:
    """
    Build the symmetric multisegment and sign dictionary for the good case.

    Args:
        m: multisegment
        T: enhanced Langlands parameter, i.e. tuple of [a, sign]

    Returns:
        - symmetric multisegment (list of segments)
        - dictionary of signs
    """
    s = []
    s += [seg[:] for seg in m]
    s += [[-x[0], x[0]] for x in T]
    s += [[-seg[1], -seg[0]] for seg in m]

    E = {int(2 * x) : e for x, e in T} #we use int(2 * x) for the keys to avoid type issues with half-integers

    return s, E

def langlands_data_good(s: List[Segment], E: Dict[Any, int]) -> LanglandsData:
    """
    Extract Langlands data (m, T) from symmetric Langlands datas.

    Args:
        s: symmetrical multisegment
        E: sign dictionary

    Returns:
        - m: multisegment
        - T: enhanced Langlands parameter ([a], sign)
    """
    m = []
    T = []
    for a, b in s:
        center = a + b
        if center < 0:
            m.append([a, b])
        elif center == 0:
            T.append([b, E[int(2*b)]])
    return tuple(m), tuple(T)

def labeled_segments(segments: List[Segment]) -> List[List]:
    """
    Label each segment [a, b] with a type:
        - 1 if a + b > 0
        - -1 if a + b < 0
        - 0 if a + b == 0

    For segments [a, b] with a + b == 0, assign balanced signs:
        - If there are an even number: half get 1, half get -1
        - If odd: one gets 0, others split evenly between 1 and -1

    Segments are then sorted for processing.
    
    Args:
        segments: list of segments [a, b]

    Returns:
        labeled: list of [a, b, t] where t is the label
    """
    # Initial labeling
    labeled = [[a, b, 1 if a + b > 0 else -1 if a + b < 0 else 0] for a, b in segments]
        
    # Group centered segments (a + b == 0) by their 'a' value
    centered_by_start = defaultdict(list)
    for idx, (a, b, t) in enumerate(labeled):
        if t == 0:
            centered_by_start[a].append(idx)

    # Assign balanced labels within each group
    for _, indices in centered_by_start.items():
        n = len(indices)
        if n % 2 == 1:
            mid = indices[n // 2]
            labeled[mid][2] = 0
            left_half = indices[:n // 2]
            right_half = indices[n // 2 + 1:]
        else:
            left_half = indices[:n // 2]
            right_half = indices[n // 2:]
        for idx in left_half:
            labeled[idx][2] = 1
        for idx in right_half:
            labeled[idx][2] = -1

    # Final sort
    sort_labeled_segments(labeled)

    return labeled

def is_stop_same_type(seg: LabeledSegment) -> bool:
    """
    Check the stop condition in the same parity case.
    """
    return seg == [0, 0, 1] or seg == [0, 0, 0]

def is_stop_opposite_type(seg: LabeledSegment, E: Dict[Any, int]) -> bool:
    """
    Check the stop condition in the opposite parity case.
    """
    return (
        seg == [1/2, 1/2, 1]
        or (seg == [-1/2, 1/2, 0] and E.get(int(1)) == -1)
        or (seg == [-1/2, 1/2, 1] and E.get(int(1)) == -1)
    )

def is_stop_segment(seg: LabeledSegment, E: Dict[Any, int]) -> bool:
    """
    Check the stop condition in all cases.

    Args:
        seg: a labeled segment [a, b, t]
        E: sign dictionary

    Returns:
        True if seg should stop the chain construction.
    """
    return is_stop_same_type(seg) or is_stop_opposite_type(seg, E)

def construct_dual_chains(
    s_labeled: List[List],
    E: Dict[Any, int]
) -> Tuple[List[int], List[int]]:
    """
    Construct the chain and its dual chain for the current iteration.

    Args:
        s_labeled: labeled symmetric multisegment [a, b, t]
        E: sign dictionary

    Returns:
        chain: list of indices of the segments in the initial sequence
        chain_dual: list of indices of the dual of these segments
    """
    #Dual of a labeled segment
    def dual_of_segment(seg: List) -> List:
        a, b, t = seg
        if a + b > 0:
            return [-b, -a, -1]
        elif a + b < 0:
            return [-b, -a, 1]
        elif t == -1:
            return [a, b, 1]
        else:
            return [a, b, t]

    #Construction of the chain
    chain = []
    i = index_max_end(s_labeled)
    e = s_labeled[i][1]
    chain.append(i)

    stop = is_stop_segment(s_labeled[i], E)

    while not stop:
        for j in range(i + 1, len(s_labeled)):
            seg_i = s_labeled[i]
            seg_j = s_labeled[j]
            if seg_j[1] == e - 1:
                if seg_i[0] + seg_i[1] == 0 and seg_j[0] + seg_j[1] == 0:
                    if E[int(2*seg_i[1])] * E[int(2*seg_j[1])] == -1:
                        chain.append(j)
                        i = j
                        e = seg_j[1]
                        stop = is_stop_segment(seg_j, E)
                        break
                else:
                    chain.append(j)
                    i = j
                    e = seg_j[1]
                    stop = is_stop_segment(seg_j, E)
                    break
        else:
            break

    # Precompute first occurrences of each segment
    segment_to_index = {}
    for idx, seg in enumerate(s_labeled):
        key = tuple(seg)
        if key not in segment_to_index:
            segment_to_index[key] = idx

    # Construct dual chain
    chain_dual = []
    for idx in chain:
        dual_seg = dual_of_segment(s_labeled[idx])
        dual_idx = segment_to_index.get(tuple(dual_seg))
        if dual_idx is not None:
            chain_dual.append(dual_idx)

    return chain, chain_dual

def compute_first_dual_segment(m: List[List], E: Dict[Any, int], s: List[int]) -> Tuple[List[Segment], Dict[Any, int]]:
    """
    Compute the first segment(s) of the dual and their associated sign(s).

    Returns:
        - dual_segments: list of segments to add to the dual
        - updated_signs: signs of these segments
    """
    e = m[s[0]][1]
    last_seg = m[s[-1]]

    if is_stop_same_type(last_seg):
        nb_centered = sum(1 for seg in m if seg[0] + seg[1] == 0)
        sign = (-1) ** (nb_centered + 1)
        base_sign = E.get(int(0), 1)
        return [[-e, e]], {int(2*e): sign * base_sign}

    elif is_stop_opposite_type(last_seg,E):
        nb_centered = sum(1 for seg in m if seg[0] + seg[1] == 0)
        sign = (-1) ** nb_centered
        return [[-e, e]], {int(2*e): sign}

    else:
        b = m[s[-1]][1]
        return [[b, e], [-e, -b]], {}

def update_segments_and_signs(
    s_labeled: List[LabeledSegment],
    E: Dict[Any, int],
    chain: List[int],
    chain_dual: List[int]
) -> Tuple[List[Segment], Dict[Any, int]]:
    """
    Update the symmetric multisegment and sign dictionary after one step of the algorithm.

    This function modifies:
      - the right endpoints of segments in the chain (decremented),
      - the left endpoints of dual segments in the chain_dual (incremented),
      - and removes degenerate segments [a, b] with a > b.

    Args:
        s_labeled: list of labeled segments [a, b, t]
        E: dictionary of signs for the symmetric parameter
        chain: indices of segments contributing to the dual segment
        chain_dual: indices of their dual counterparts

    Returns:
        - Updated multisegment as a list of [a, b] (no labels)
        - Updated sign dictionary
    """
    s_updated = [seg[:] for seg in s_labeled]  # safe copy

    for i in chain:
        s_updated[i][1] -= 1
    for i in chain_dual:
        s_updated[i][0] += 1

    # Remove degenerate segments
    s_result = [[a, b] for a, b, _ in s_updated if a <= b]

    # Determine sign factor for the step
    last = s_labeled[chain[-1]]
    sgn = -1 if is_stop_segment(last,E) else 1

    # Update signs
    E_updated = {}

    for i in chain:
        a, b, t = s_labeled[i]
        if t == 0 and b >= 1:
            E_updated[int(2*(b - 1))] = sgn * E[int(2*b)]
        elif t == 1 and a + b == 0 and b >= 1:
            E_updated[int(2*(b - 1))] = sgn * E[int(2*b)]
        elif a + b == 1 and i not in chain_dual and a != b:
            e = b - 1
            E_updated[int(2*e)] = sgn * (-1) * E.get(int(2*e), 1)

    for a, b in s_result:
        if a + b == 0 and int(2*b) not in E_updated:
            E_updated[int(2*b)] = sgn * E[int(2*b)]

    return s_result, E_updated

@validated(validate_data_good)
def AD_good(m: Multisegment, T: EnhancedParameter, verbose: bool = False) -> LanglandsData:
    """
    Compute the Aubert–Zelevinsky dual in the good parity case.

    Args:
        m: multisegment
        T: enhanced Langlands parameter (tuple of [a, ±1])
        verbose: if True, prints detailed steps

    Returns:
        Langlands data (m', T') after duality
    """
    #Symmetrisation of the datas
    s, E = symmetrisation_good(m, T)

    if verbose:
        _print_step("Good Case", 0, message="Symmetrisation", remaining=tuple(s), signs=format_sign_dict(E))

    dual_segments = []
    dual_signs = {}

    step = 1

    while s:
        #Labelling of the segments
        s_labeled = labeled_segments(s)

        #Construction of the chain of segments and the dual chain
        chain, chain_dual = construct_dual_chains(s_labeled, E)

        #Computing the first segment of the dual
        new_dual, new_signs = compute_first_dual_segment(s_labeled, E, chain)

        #Updating the multisegment and the signs according to the chain and the dual chain
        s, E = update_segments_and_signs(s_labeled, E, chain, chain_dual)

        if verbose:
            _print_step("Good Case", step, labeled_seg=tuple(s_labeled), dual=tuple(new_dual), sign_dual=format_sign_dict(new_signs), remaining=tuple(s), signs=format_sign_dict(E))
            step += 1

        dual_segments.extend(new_dual)
        dual_signs.update(new_signs)

    return langlands_data_good(dual_segments, dual_signs)

# =========================================================
# Common interface
# =========================================================

def AD(
    case: str = None,
    *args,
    verbose: bool = False
) -> Union[Tuple[LanglandsData, LanglandsData], LanglandsData, Multisegment]:
    """
    Unified interface for the Aubert–Zelevinsky involution algorithm.
    
    Handles the following cases:
        - "gl": general linear group (GL(n)), using Moeglin–Waldspurger algorithm
        - "ugly": classical group with ρ ugly
        - "bad": classical group with ρ bad
        - "good": classical group with ρ good

    If case is None, it is inferred from the arguments.

    Args:
        case: one of "gl", "ugly", "bad", "good"
        *args: the Langlands data
        verbose: if True, prints internal processing steps

    Returns:
        dual Langlands data
    """
    if isinstance(case, (list, tuple)):
        args = (case,) + args
        case = None

    def looks_like_enhanced(T):
        return all(isinstance(x, (list, tuple)) and len(x) == 2 for x in T)
    
    if case is not None:
        case = case.lower()

    if case is None:
        # Try to infer case based on argument pattern
        if len(args) == 1:
            case = "gl"
        elif len(args) == 2:
            m, T = args
            if T == ():
                raise ValueError("T is empty: cannot infer between 'good' and 'bad'. Please specify the case explicitly.")
            case = "good" if looks_like_enhanced(T) else "bad"
        elif len(args) == 3:
            case = "ugly"
        else:
            raise ValueError("Cannot infer case from arguments. Please provide 'case' explicitly.")
        
        if verbose:
            print(f"Inferred case: {case}")
            print()
    
    # Dispatch to actual implementation
    if case in {"gln", "gl", "gl_n"}:
        if len(args) != 1:
            raise TypeError("AD('gl', m) attendu.")
        return AD_gl(args[0], verbose=verbose)

    elif case == "ugly":
        if len(args) != 3:
            raise TypeError("AD('ugly', m, m_contr, T) attendu.")
        return AD_ugly(args[0], args[1], args[2], verbose=verbose)

    elif case == "bad":
        if len(args) != 2:
            raise TypeError("AD('bad', m, T) attendu.")
        return AD_bad(args[0], args[1], verbose=verbose)

    elif case == "good":
        if len(args) != 2:
            raise TypeError("AD('good', m, T) attendu.")
        return AD_good(args[0], args[1], verbose=verbose)

    else:
        raise ValueError(f"Unknown case '{case}'. Expected one of: 'gl', 'ugly', 'bad', 'good'.")