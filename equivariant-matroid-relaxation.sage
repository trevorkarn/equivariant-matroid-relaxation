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