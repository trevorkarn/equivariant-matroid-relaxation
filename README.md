# Equivariant matroid relaxation

Ferroni, Nasr, and Vecchi show a remarkable connection between the Kazhdan-Lusztig polynomial of a matroid and the operation of matroid relaxation. Forthcoming work shows this extends to an equivariant version. This code computes the decomposition into irreducible Mathieu representations of the polynomial corresponding to various Steiner systems.

For example, we can compute the equivariant Kazhdan-Lusztig polynomial of the largest known Steiner system *S(5, 8, 24)*:

	sage: load('equivariant-matroid-relaxation.sage')
	sage: steiner_system_KL_coeff(d, k, n, 1).values()
	[735, 15, 0, 0, 0, 15, 0, 0, 63, 7, 3, 0, -1, 7, 3, 0, 0, 0, 0, 0, 0, 0, -2, 1, -1, -1]
	sage: steiner_system_KL_coeff(d,k,n,2).values()
	[4830, 6, -5, 1, 1, 54, 0, 0, 110, 2, 6, 0, 0, 6, 2, 0, 0, -2, -2, 0, 0, -1, 1, 0, 0, 0]

We can also decompose it:

	sage: steiner_system_KL_coeff(d,k,n,2).decompose()
	((1,
	  Character of Mathieu group of degree 24 and order 244823040 as a permutation group),
	 (1,
	  Character of Mathieu group of degree 24 and order 244823040 as a permutation group),
	 (1,
	  Character of Mathieu group of degree 24 and order 244823040 as a permutation group))

For more information about `ClassFunction`s in SageMath, [see the SageMath documentation.](https://doc.sagemath.org/html/en/reference/groups/sage/groups/class_function.html)

**WARNING**

We use SageMath functions which are wrappers around GAP functions. GAP does not provide conjugacy classes
in any standard order, nor is it consistent from function-call to function-call. (For more, see [GAP's documentation](https://www.gap-system.org/Manuals/doc/ref/chap39.html#X7D474F8F87E4E5D9))). Thus, values obtained via this code will
agree with the forthcoming paper containing these results as a set, but will not neccessarily be listed in the same
order.
