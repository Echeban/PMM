function [E1,E2,G12,v12,v23,G23] = PMM(EA,ET,GA,vA,vT,EM,vM,Vf)
    // Periodic Microstructure Micromechanics (PMM) transversely isotropic fiber. 
    // (c) Ever J. Barbero (1994-2018)
    // Equations taken from (please cite in your work): 
    //   Barbero, E.J. and Luciano, R., "Micromechanical formulas for the
    //   relaxation tensor of linear viscoelastic composites with transversely
    //   isotropic fibers", International Journal of Solids and Structures,
    //   32(13):1859--1872, 1995.
    //
    // fiber coefficients
    f =ET/EA
    Delta = (1-2*vA^2*f-vT^2-2*vA^2*vT*f)/(EA*ET^2) // (A.3)
    C_11 = (1-vT^2)/(ET^2*Delta)
    C_22 = (1-vA^2*f)/(EA*ET*Delta)
    C_33 = C_22
    C_12 = (vA*f+vA*vT*f)/(ET^2*Delta)
    C_13 = C_12
    C_23 = (vT + vA^2*f)/(EA*ET*Delta)
    C_44 = ET/(2*(1+vT))
    C_55 = GA
    C_66 = C_55
    // multiply by EA*ET^2 to avoid rounding errors
//    Delta = (1-2*vA^2*f-vT^2-2*vA^2*vT*f)               
//    C_11 = (1-vT^2)*EA/(Delta)
//    C_22 = (1-vA^2*f)*ET/(Delta)
//    C_33 = C_22
//    C_12 = (vA*f+vA*vT*f)*EA/(Delta)
//    C_13 = C_12
//    C_23 = (vT + vA^2*f)*ET/(Delta)
//    C_44 = ET/(2*(1+vT))
//    C_55 = GA
//    C_66 = C_55
    // matrix coefficients
    lam_m = (EM*vM)/((1+vM)*(1-2*vM))
    mu_m = EM/2/(1+vM)//isotropic matrix
    // geometry
    S3 = 0.49247 - 0.47603*Vf - 0.02748*Vf^2
    S6 = 0.36844 - 0.14944*Vf - 0.27152*Vf^2
    S7 = 0.12346 - 0.32035*Vf + 0.23517*Vf^2
    // a_values
    a1 = 4*mu_m^2 - 2*mu_m*C_33 + 6*lam_m*mu_m - 2*C_11*mu_m - 2*mu_m*C_23 + C_23*C_11 + 4*lam_m*C_12 - 2*C_12^2 - lam_m*C_33 -2*C_11*lam_m + C_11*C_33 - lam_m*C_23
    a2 = 8*mu_m^3- 8*mu_m^2*C_33 + 12*mu_m^2*lam_m -4*mu_m^2*C_11 - 2*mu_m*C_23^2 + 4*mu_m*lam_m*C_23 + 4*mu_m*C_11*C_33 - 8*mu_m*lam_m*C_33 - 4*mu_m*C_12^2 + 2*mu_m*C_33^2 -4*mu_m*C_11*lam_m + 8*mu_m*lam_m*C_12 + 2*lam_m*C_11*C_33 + 4*C_12*C_23*lam_m - 4*C_12*C_33*lam_m - 2*lam_m*C_11*C_23 - 2*C_23*C_12^2 + C_23^2*C_11 + 2*C_33*C_12^2 - C_11*C_33^2 + lam_m*C_33^2 - lam_m*C_23^2 
    a3 = ((4*mu_m^2 + 4*lam_m*mu_m - 2*C_11*mu_m - 2*mu_m*C_33 - C_11*lam_m - lam_m*C_33 - C_12^2)/a2)  + ((C_11*C_33 + 2*lam_m*C_12)/a2) - ((S3-((S6)/(2-2*vM)))/mu_m) 
    a4 = -1*((-2*mu_m*C_23 + 2*lam_m*mu_m - lam_m*C_23 - C_11*lam_m - C_12^2 + 2*lam_m*C_12 + C_11*C_23)/a2) + (S7)/(mu_m*(2-2*vM))
    // C_values
    C11t = lam_m + 2*mu_m - Vf*(-a4^2 + a3^2)*inv(-1*(((2*mu_m + 2*lam_m - C_33 - C_23)*(a4^2-a3^2))/a1) + ((2*(a4-a3)*(lam_m-C_12)^2)/a1^2)) 
    C12t = lam_m + Vf*(((lam_m-C_12)*(a4-a3))/a1)*inv(1*(((2*mu_m + 2*lam_m - C_33 - C_23)*(a3^2-a4^2))/a1) + ((2*(a4-a3)*(lam_m-C_12)^2)/a1^2))
    C22t = lam_m + 2*mu_m - Vf*(((2*mu_m + 2*lam_m - C_33 - C_23)*a3/a1) -((lam_m-C_12)^2/a1^2))*inv(1*(((2*mu_m + 2*lam_m - C_33 - C_23)*(a3^2-a4^2))/a1) + ((2*(a4-a3)*(lam_m-C_12)^2)/a1^2))
    C23t = lam_m + Vf*(((2*mu_m + 2*lam_m - C_33 - C_23)*a4/a1) -((lam_m-C_12)^2/a1^2))*inv(1*(((2*mu_m + 2*lam_m - C_33 - C_23)*(a3^2-a4^2))/a1) + ((2*(a4-a3)*(lam_m-C_12)^2)/a1^2))
    C44t = mu_m - Vf*inv((2/(2*mu_m - C_22 + C_23)) - inv(mu_m)*(2*S3 - (4*S7/(2-2*vM))))
    C66t = mu_m - Vf*inv(inv(mu_m-C_66) - S3/mu_m)
    // C_total
    C11 = C11t
    C12 = C12t
    C13 = C12t
    C22 = (3/4)*C22t + (1/4)*C23t +(1/4)*C44t // corrected 
    C33 = C22
    C23 = (1/4)*C22t + (3/4)*C23t -(1/4)*C44t // corrected
    C55 = C66t
    C66 = C66t
    C44 = (C22-C23)/2 // do not use C44t here, typo in (A.6) App. A.2, ISBN:978-1-138-19680-3 
    // stiffness
    C = [C11 C12 C13 0 0 0; C12 C22 C23 0 0 0; C13 C23 C33 0 0 0; 0 0 0 C44 0 0; 0 0 0 0 C55 0; 0 0 0 0 0 C66]
    // compliance
    S = inv(C);
    // Elastic properties
    E1 = 1/S(1,1) 
    E2 = 1/S(2,2)
    v12 = -S(2,1)/S(1,1)
    v23 = -S(3,2)/S(2,2) 
    G12 = 1/S(6,6)
    G23 = 1/S(4,4)  //redundant, see (A.7) App. A.2, ISBN:978-1-138-19680-3 
//    DD = 2/(2*mu_m - C_22 + C_23) - (2*S3*(2-2*vM)-4*S7)/mu_m/(2-2*vM)
//    DD = 2/(2*mu_m - C_22 + C_23) - (2*S3 - (4*S7/(2-2*vM)))/mu_m //same
//    G23 = mu_m - Vf/DD
endfunction
