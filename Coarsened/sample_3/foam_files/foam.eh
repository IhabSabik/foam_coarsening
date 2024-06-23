torus

// we use symbolic periods for all foams so that we can change them
// with applymat (see foam.ed).  The .fe file _must_ define these nine
// parameters before including this file.
periods
p1x p1y p1z
p2x p2y p2z
p3x p3y p3z

// to be applied to all edges of valence 2
quantity dihooke energy modulus 0 method dihedral_hooke

// quantities for stresses
#define nlen sqrt(x4^2+x5^2+x6^2)
quantity Sxxq info_only global_method facet_general_integral
    scalar_integrand:  (nlen-x4^2/nlen)
quantity Syyq info_only global_method facet_general_integral
    scalar_integrand:  (nlen-x5^2/nlen)
quantity Szzq info_only global_method facet_general_integral
    scalar_integrand:  (nlen-x6^2/nlen)
quantity Syzq info_only global_method facet_general_integral
    scalar_integrand: -x5*x6/nlen
quantity Szxq info_only global_method facet_general_integral
    scalar_integrand: -x6*x4/nlen
quantity Sxyq info_only global_method facet_general_integral
    scalar_integrand: -x4*x5/nlen

// quantities for anisotropy tensor Q for bubble shape
quantity Qxxq info_only global_method facet_general_integral
    scalar_integrand:  (nlen/3-x4^2/nlen)
quantity Qyyq info_only global_method facet_general_integral
    scalar_integrand:  (nlen/3-x5^2/nlen)
quantity Qzzq info_only global_method facet_general_integral
    scalar_integrand:  (nlen/3-x6^2/nlen)
quantity Qyzq info_only global_method facet_general_integral
    scalar_integrand: -x5*x6/nlen
quantity Qzxq info_only global_method facet_general_integral
    scalar_integrand: -x6*x4/nlen
quantity Qxyq info_only global_method facet_general_integral
    scalar_integrand: -x4*x5/nlen

// quantities for anisotropy tensor A for edge structure
#define tlen sqrt(x4^2+x5^2+x6^2)
quantity Axxq info_only global_method edge_general_integral
    scalar_integrand:  (x4^2/tlen-tlen/3)
quantity Ayyq info_only global_method edge_general_integral
    scalar_integrand:  (x5^2/tlen-tlen/3)
quantity Azzq info_only global_method edge_general_integral
    scalar_integrand:  (x6^2/tlen-tlen/3)
quantity Ayzq info_only global_method edge_general_integral
    scalar_integrand: x5*x6/tlen
quantity Azxq info_only global_method edge_general_integral
    scalar_integrand: x6*x4/tlen
quantity Axyq info_only global_method edge_general_integral
    scalar_integrand: x4*x5/tlen

define body attribute old_target real
define body attribute radius real
