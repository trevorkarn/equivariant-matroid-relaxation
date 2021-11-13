Sym = SymmetricFunctions(QQ)

def equivariant_steiner_hyperplane(d, k, n):
    r"""
    Find a hyperplane of the Steiner system `S(d, k, n)`.
    
    The Steiner system `S(d, k, n)` can be thought of as a matroid.
    Endowing it with an action of the Mathieu group `M_n` yields
    an equivariant matroid. The :class:`Sage implementation of `M_n` 
    <sage.groups.perm_gps.permgroup_named.MathieuGroup>` is as a 
    :class:`sage.groups.perm_gps.permgroup_named.PermutationGroup_unique`.
    In order to find a representative hyperplane, we must associate
    the permutation indices in `M_n` to a hyperplane stabilized by the
    correct subgroup.
        
    .. WARNING::
    
        This code does not check if `S(d, k, n)` is a valid Steiner system.
        Passing invalid parameters may lead to mathematically incorrect 
        results.
    
    INPUT:
    
    - ``k`` -- the block size of the Steiner system
    - ``d`` -- the size of the subset contained in a unique hyperplane
    - ``n`` -- a positive integer in `\{ 9, 10, 11, 12, 21, 22, 23, 24 \}`
    
    OUTPUT:
    
    - ``H`` -- a list representing the indices of a representative hyperplane
    
    EXAMPLES:
    
    Find a representative hyperplane for the Steiner system `S(8, 5, 24)`::
    
        sage: equivariant_steiner_hyperplane(8, 5, 24)
        [1, 2, 3, 4, 5, 8, 11, 13]
        
    """
    
    G = MathieuGroup(n)

    from itertools import combinations
    
    # The existence of a hyperplane containing the subset `\{1,\ldots, k\}`
    # is guaranteed.
    guaranteed_hyperplane_subset = list(range(1, k+1))
    
    # Get the size of the orbit of ``H`` which we desire.
    # This should be the number of hyperplanes of `S(k, d, n)` for G to
    # act transitively. 
    orbit_size = binomial(n, k)/binomial(d, k)
    
    # Get the size of the stabilizer of the hyperplane we are looking for
    # by the orbit stabilizer theorem
    stab_size = G.order()/orbit_size

    # The sets `A` such that `[k] \cup A` could be a hyperplane of the
    # Steiner system
    hyperplane_difference_candidates = list(combinations(range(k+1, n+1), d-k))

    for complement in hyperplane_difference_candidates:
        # define the candidate hyperplane
        H = guaranteed_hyperplane_subset + list(complement)
        
        # compute the stabilizer of the hyperplane
        stab = G.stabilizer(H, "OnSets")
        
        # if the stabilizer is the right size, we have found a
        # representative hyperplane
        if len(stab) == stab_size: return H
        
def restrict_symmetric_irreducible(G, mu):
    r"""
    Return the character of `G`.
    
    The permutation group `G` can be naturally included into a bigger
    symmetric group, say `S_N`, which acts by permuting the same 
    underlying set. This function computes the character of the restriction
    of the irreducible `S_N` representation indexed by `mu` to `G`.
    
    INPUT:
    
    - ``G`` -- an instance of a :class:`sage.groups.perm_gps.permgroup.PermutationGroup`
    - ``mu`` -- a partition of `N` - the size of ``G.domain()``
    
    OUTPUT:
    
    - ``chi`` -- a :func:`sage.groups.class_function.ClassFunction` corresponding
      to the restriction of the `S_N` irreducible representation indexed by `mu`
      
    EXAMPLES:
    
        sage: G = MathieuGroup(11)
        sage: restrict_irreducible(G, [7,2,2]).values()
        [385, 9, 1, -1, -1, -2, 0, 0, 0, 0]
        
    TESTS::
        
        sage: G = MathieuGroup(11)
        sage: restrict_irreducible(G, [7])
        Traceback (most recent call last):
        ...
        ValueError: [7] is not a partition of len(G.domain())=11
        
    """
    if len(G.domain()) != Partition(mu).size():
        raise ValueError(f"{mu} is not a partition of {len(G.domain())=}")
    
    # Gap does not provide conjugacy classes in a consistent order. Thus, we
    # create an index function in order to keep track of the conjugacy classes.
    indexing_cf_vals = range(len(G.conjugacy_classes()))
    indexing_cf = ClassFunction(G, indexing_cf_vals)

    # We will apply the Murnagan-Nakayama rule to quickly compute character values
    from sage.combinat.sf.sfa import zee
    s = Sym.schur()
    p = Sym.powersum()

    # Initialize a lookup table. Keys will be cycle types, values will be the character value.
    chi = {}

    for t in [c[0].cycle_type() for c in G.conjugacy_classes()]:
        try:
            # Use the expression of the Murnagan-Nakayama rule
            # in terms of symmetric functions
            chi[t] = p(s[mu]).monomial_coefficients()[t]*zee(t)
        except KeyError:
            # If the character value is 0, the cycle type ``t``
            # will not appear in the expansion of 
            # ``.monomial_coefficients()`` and a KeyError will
            # be raised.
            chi[t] = 0
    
    # Give a list of the character values, sorted into the order which Gap expects
    sorted_chi_vals = sorted(zip([int(indexing_cf(g.representative())) for g in G.conjugacy_classes()],
                                 [chi[c[0].cycle_type()] for c in G.conjugacy_classes()]), key = lambda x: x[0])

    return ClassFunction(G,[b for a,b in sorted_chi_vals])
    
def restrict_symmetric_repn(G, irreps):
    r"""
    Return the restriction of a symmetric group representation to ``G``.
    
    The symmetric group is semisimple, so it decomposes into a list of irreducibles.
    This function allows the user to get the restriction of a symmetric group representation
    to a permutation group ``G`` by passing a list of :class:`sage.combinat.partition.Partition``s.
    The list ``irreps`` is interpreted as the partitions indexing the decomposition of the symmetric
    group into irreducible representations.
    
    INPUT:
    
    - ``G`` -- a :class:`sage.groups.perm_gps.permgroup.PermutationGroup`
    - ``irreps`` -- a list of partitions corresponding to the decomposition
      of a symmetric group representation into irreducibles
    
    OUTPUT:
    
    - ``chi`` -- a :func:`sage.groups.class_function.ClassFunction` on ``G``
      corresponding to the restriction of the symmetric group representation
    
    EXAMPLES:
    
        sage: G = MathieuGroup(11)
        sage: restrict_symmetric_repn(G, [[7,2,2]]).values()
        [385, 9, 1, -1, -1, -2, 0, 0, 0, 0]
        
        sage: restrict_symmetric_repn(SymmetricGroup(3), [[2,1], [1,1,1]]).decompose()
        ((1, Character of Symmetric group of order 3! as a permutation group),
         (1, Character of Symmetric group of order 3! as a permutation group))
    """
    
    # initialize the 0 class function
    chi = ClassFunction(G, [0]*len(G.conjugacy_classes()))
    
    for mu in irreps:
        chi += restrict_symmetric_irreducible(G, mu)
        
    return chi
    
def slack_polynomial_coefficient(k, h, i):
    r"""
    Return the coefficient of `t^i` in the equivariant polynomial `p_{k,h}^{S_h}`.
    
    INPUT:
    
    - ``k`` -- an integer
    - ``h`` -- an integer
    - ``i`` -- an integer
    
    OUTPUT:
    
    - ``chi`` -- the character coresponding to coefficient of `t^i` in the polynomial.
    """
    
    skew_partition_outer = [h - 2 * i + 1] + [k - 2 * i + 1] * i
    skew_partition_inner = [h - 2 * i] + [k - 2 * i - 1] * (i - 1)
    
    s = Sym.schur()
    
    skew_schur = s[skew_partition_outer].skew_by(s[skew_partition_inner])
    
    chi = None
    
    for mu, c in skew_schur.monomial_coefficients().items():
        if not chi:
            chi = SymmetricGroupRepresentation(mu, implementation='specht').to_character()*c
        else:
            chi += SymmetricGroupRepresentation(mu, implementation='specht').to_character()*c
        
    return chi
    
def ind_res_slack_polynomial_coeff(k,h,i,V,W,H):
    r"""
    Compute the induction from ``V`` to ``W`` of the restriction of `p^{S_h}_{k,h}`.
    
    Let `S_h` act on the hyperplane `H`.
    

    """
    Sh = SymmetricGroup(domain=H)
    
    chi = slack_polynomial_coefficient(k, h, i)
    
    # Gap does not provide conjugacy classes in a consistent order. Thus, we
    # create an index function in order to keep track of the conjugacy classes.
    stab_index_cf = ClassFunction(V, range(len(V.conjugacy_classes())))

    sorted_res_chi_values = sorted(zip([int(stab_index_char(g.representative())) for g in H.conjugacy_classes()],
                             [chi(Sd.prod([Sd(c) for c in w[0].cycles() if c in Sd])) for w in H.conjugacy_classes()]),
                             key = lambda x: x[0])
                             
    res_p = ClassFunction(stab, [b for a,b, in sorted_res_chi_values])
    
    return res_p.induct(W)